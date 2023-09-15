% Generate specifications and related data for the 1/f project
%
% Cristina Gil Avila, TUM, 15.9.2023

clear all, close all;

%% Generate all the specifications
freq_range = {'2-40', '1-100', '30-45'};
taper = {'hanning', 'dpss'};
zero_padding =  {'0', '5', '10'};
average_psd = {'no', 'yes'};
fooof_knee =  {'no', 'yes'};

[a,b,c,d,e] = ndgrid(freq_range,taper,zero_padding,average_psd,fooof_knee);
nSpec = length(a(:));
specs_id = num2cell(1:nSpec);
specs = [specs_id',a(:),b(:),c(:),d(:),e(:)];
s = cell2table(specs,'VariableNames',{'spec_id','freq_range','taper','zero_padding','average_psd','fooof_knee'});


%% Initialization of paths and other settings

% Import paths
fid = fopen('paths.json');
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
paths = jsondecode(str);
clear fid raw str

% Add functions from the discover-eeg to read preprocessed data
addpath(genpath(paths.discover_eeg));
addpath(paths.fieldtrip);
ft_defaults;
% Add matlab FOOOF functions
addpath('fooof_matlab');

% Load preprocessing parameters
fparams = '../data/blinded/derivatives_v2023_08_18/params.json';
fid = fopen(fparams);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
params = jsondecode(str);
clear fid raw str

% Define the results path
params.RawDataPath = '../data/blinded';
params.PreprocessedDataPath = '../data/blinded/derivatives_v2023_08_18';
params.AtlasPath = 'schaefer_parcellations/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_1mm.Centroid_RAS.csv';
params.SourcePath = '../results/source';
path_power_pfc = '../results/power/PFC';
path_power_s1 = '../results/power/S1';
path_power_vis = '../results/power/VIS';

if ~exist(path_power_pfc)
    mkdir(path_power_pfc)
end
if ~exist(path_power_s1)
    mkdir(path_power_s1)
end
if ~exist(path_power_vis)
    mkdir(path_power_vis)
end


% Parcellation based on Schaefer atlas and source model
atlas_table = readtable(params.AtlasPath);
parcellation.pos = [atlas_table.R, atlas_table.A, atlas_table.S];
parcellation.unit = 'mm';
parcellation.coordsys = 'mni';
parcellation.ROI = atlas_table.ROILabel;
parcellation.ROIlabel = atlas_table.ROIName;
PFC_mask = contains(atlas_table.ROIName,{'Cinga','PFCm','PFCmp'});

% Source model: centroid positons from Schaefer atlas
cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = [atlas_table.R, atlas_table.A, atlas_table.S];
cfg.unit = 'mm';
cfg.headmodel = params.HeadModelPath;
sourcemodel_atlas = ft_prepare_sourcemodel(cfg);
sourcemodel_atlas.coordsys = 'mni';

% Load all subject ids
participants = readtable(fullfile(params.RawDataPath,'participants_rand.tsv'),'Filetype','text');
participant_id = participants.participant_id;
task = 'closed';
nSubj = height(participants);
%%
% Loop over specifications
for iSpec=1:nSpec
       
    % S1. Frequency range
    freq_range = strsplit(s.freq_range{iSpec},'-');
    freq_range = cellfun(@str2num,freq_range);
    params.FreqBand.fullSpectrum = freq_range;  
    
    % Loop over subjects
    for iSubj=1:nSubj

        bidsID = participant_id{iSubj};
        bidsID = [bidsID '_task-' task];

        % ----- Load EEG preprocessed data -----
        data = load_preprocessed_data(params,bidsID);

        % ----- Compute source reconstruction ------
        source = compute_spatial_filter(params,bidsID,'fullSpectrum');

        % ----- Extract virtual channel data -----
        cfg = [];
        cfg.parcellation = 'ROI';
        vdata_trials = ft_virtualchannel(cfg,data,source,parcellation);

        % ----- Estimate power spectra at the source level -----
        cfg = [];
        cfg.foilim = [1 100];
        cfg.method = 'mtmfft';
        % S2. Taper
        cfg.taper = s.taper{iSpec};   
        switch s.taper{iSpec}
            case 'dpss'
                cfg.tapsmofrq = 1;
        end
        % S3. Padding
        switch s.zero_padding{iSpec}
            case '0'
                % Do not pad, don't do anything here
            otherwise
                cfg.pad = str2double(s.zero_padding{iSpec}); 
                cfg.padtype = 'zero';   
        end       
        cfg.output = 'pow';
        cfg.keeptrials ='no';
        power = ft_freqanalysis(cfg, vdata_trials);

        % ----- Extract power at the PFC ------
        cfg = [];
        cfg.channel = find(PFC_mask);
        % S4. Average PSD over channels
        cfg.avgoverchan = s.average_PSD{iSpec};
        power_PFC = ft_selectdata(cfg,power);
          
        % ----- Model power spectrum with FOOOF ------
        % Create a fooof object with default settings
        % S5. Knee parameter
        switch('fooof_knee')
            case 'no'
                fm = fooof('freq_range',freq_range);
            case 'yes'
                fm = fooof('freq_range',freq_range,'aperiodic_mode','knee');
        end
        % Add original freq and power spectrum to the fooof model
        fm.orig_freq = power_PFC.freq;
        fm.orig_pow = power_PFC.powspctrm;
        
        % Select a subset of frequencies defined in freq_range and transform power to log space
        f = and(freq >= fm.freq_range(1),freq<=fm.freq_range(2));
        fm.freq = freq_orig(f);
        fm.pow = log10(pow_orig(f));

        % Fit fooof
        fm = fit_fooof(fm);
    end


    
end