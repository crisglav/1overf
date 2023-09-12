% Estimate power spectrum at mPFC from preprocessed data
clear all;
close all;

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

% Define the frequency band of interest
params.FreqBand.fullSpectrum = [2 40];

%% Parcellation based on Schaefer atlas and source model
atlas_table = readtable(params.AtlasPath);
parcellation.pos = [atlas_table.R, atlas_table.A, atlas_table.S];
parcellation.unit = 'mm';
parcellation.coordsys = 'mni';
parcellation.ROI = atlas_table.ROILabel;
parcellation.ROIlabel = atlas_table.ROIName;

PFC_mask = contains(atlas_table.ROIName,{'Cinga','PFCm','PFCmp'});
% parcellation.PFC = PFC_mask;

S1_mask = contains(atlas_table.ROIName,{'PostC'});
% parcellation.S1 = S1_mask;

VIS_mask = contains(atlas_table.ROIName,{'ExStr_','Striate','StriCal'});
% parcellation.VIS = VIS_mask;

% Source model: centroid positons from Schaefer atlas
cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = [atlas_table.R, atlas_table.A, atlas_table.S];
cfg.unit = 'mm';
cfg.headmodel = params.HeadModelPath;
sourcemodel_atlas = ft_prepare_sourcemodel(cfg);
sourcemodel_atlas.coordsys = 'mni';

%% Load all subject ids
participants = readtable(fullfile(params.RawDataPath,'participants_rand.tsv'),'Filetype','text');
task = 'closed';
n = height(participants);

%% 
% parpool(10);
participant_id = participants.participant_id;
parfor iSubj=1:n
    
    bidsID = participant_id{iSubj};
    bidsID = [bidsID '_task-' task];

    % Load EEG preprocessed data
    data = load_preprocessed_data(params,bidsID);

    % Compute source reconstruction
    source = compute_spatial_filter(params,bidsID,'fullSpectrum');
    
    % Extract virtual channel data
    cfg = [];
    cfg.parcellation = 'ROI';
    vdata_trials = ft_virtualchannel(cfg,data,source,parcellation);

    % Estimate power spectra at the source level
    cfg = [];
    cfg.foilim = [1 100];
    cfg.method = 'mtmfft';
    cfg.taper = 'dpss';
    cfg.tapsmofrq = 1;
    cfg.pad = 5;
    cfg.padtype = 'zero';
    cfg.output = 'pow';
    cfg.keeptrials ='no';
    power = ft_freqanalysis(cfg, vdata_trials);

    % Extract power at the PFC
    cfg = [];
    cfg.channel = find(PFC_mask);
    power_PFC = ft_selectdata(cfg,power);
    freq = power_PFC.freq;
    pow = power_PFC.powspctrm;
    parsave(fullfile(path_power_pfc, [bidsID '_PFC.mat']),freq,pow);
%     save(fullfile(path_power_pfc, [bidsID '_PFC.mat']),'freq','pow');

    % Extract power at the S1
    cfg = [];
    cfg.channel = find(S1_mask);
    power_S1 = ft_selectdata(cfg,power);
    freq = power_S1.freq;
    pow = power_S1.powspctrm;
    parsave(fullfile(path_power_s1, [bidsID '_S1.mat']),freq,pow);
%     save(fullfile(path_power_s1, [bidsID '_S1.mat']),'freq','pow');

    % Extract power at the VIS
    cfg = [];
    cfg.channel = find(VIS_mask);
    power_VIS = ft_selectdata(cfg,power);
    freq = power_VIS.freq;
    pow = power_VIS.powspctrm;
    parsave(fullfile(path_power_vis, [bidsID '_VIS.mat']),freq,pow);
%     save(fullfile(path_power_vis, [bidsID '_VIS.mat']),'freq','pow');
    
end

