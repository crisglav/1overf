% Generate specifications and related data for the 1/f project
%
% Cristina Gil Avila, TUM, 15.9.2023

clear all, close all;

%% Generate all the specifications
taper = {'hanning', 'dpss'};
zero_padding =  {'0', '5'}; % Optional include 10
average_psd = {'no', 'yes'};
fooof_range = {'2-40', '1-100', '40-60'};
fooof_knee =  {'no', 'yes'};

[a,b,c,d,e] = ndgrid(taper,zero_padding,average_psd,fooof_range,fooof_knee);
nSpec = length(a(:));
specs_id = num2cell(1:nSpec);
specs = [specs_id',a(:),b(:),c(:),d(:),e(:)];
s = cell2table(specs,'VariableNames',{'spec_id','taper','zero_padding','average_psd','fooof_range','fooof_knee'});

%% Settings
% Add fieldtrip
addpath('/rechenmagd4/toolboxes_and_functions/fieldtrip');
ft_defaults;

% Load preprocessing parameters
fparams = '../data/blinded/derivatives_v2023_08_18/params.json';
params = load_params(fparams);

% Define the atlas for source localization, create ROIs mask and precompute
% the source model
params.AtlasPath = 'schaefer_parcellations/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_1mm.Centroid_RAS.csv';
params = create_parcellation(params);

% Define the frequency band in which to compute the spatial filter
params.FreqBand.fullSpectrum = [0.5 100.5];

% Define the number of cores for parallelization
params.Ncores = 2;
if(isempty(gcp('nocreate')))
    parObj = parpool(params.Ncores);
end

% Define the results path
params.VdataPath = '../results/vdata';
if ~exist(params.VdataPath)
    mkdir(params.VdataPath)
end
params.PowerPFCPath = '../results/power/PFC';
if ~exist(params.PowerPFCPath)
    mkdir(params.PowerPFCPath)
end
params.FOOOFPath = '../results/fooof_matlab/PFC';
if ~exist(params.FOOOFPath)
    mkdir(params.FOOOFPath)
end

%% Initialization of paths and other settings


% Load all subject ids
participants = readtable(fullfile(params.RawDataPath,'participants_rand.tsv'),'Filetype','text');
participant_id = participants.participant_id;
task = 'closed';
nSubj = height(participants);

%% Loop over all specifications
for iSpec=1:nSpec
      
    
    % Loop over subjects
    for iSubj=1:nSubj

        bidsID = participant_id{iSubj};
        bidsID = [bidsID '_task-' task];
        
        % ----- Load virtual channel data  ------
        load(fullfile(params.VdataPath,[bidsID '_vdata.mat']));

        % ----- Estimate power spectra at the source level -----
        cfg = [];
        cfg.foilim = [1 100];
        cfg.method = 'mtmfft';
        % S1. Taper
        cfg.taper = s.taper{iSpec};   
        switch s.taper{iSpec}
            case 'dpss'
                cfg.tapsmofrq = 1;
        end
        % S2. Padding
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

        % ----- Extract power at the PFC -----
        cfg = [];
        cfg.channel = find(params.PFC_mask);
        % S3. Average PSD over channels / exponents
        cfg.avgoverchan = s.average_psd{iSpec};
        power_PFC = ft_selectdata(cfg,power);
          
        % ----- Model power spectrum with FOOOF ------
        % Initialize a fooof object with settings depending on the specification
        % S4. Fooof range
        % S5. Knee parameter
        fooof_range = strsplit(s.fooof_range{iSpec},'-');
        fooof_range = cellfun(@str2num,fooof_range);
        switch (s.fooof_knee{iSpec})
            case 'no'
                fm = fooof('freq_range',fooof_range);
            case 'yes'
                fm = fooof('freq_range',fooof_range,'aperiodic_mode','knee');
        end
        % Add data in freq_range to the fooof model
        fm = fm.add_data(power_PFC.freq,power_PFC.powspctrm,fooof_range); 
        
        % FIX = fit data per channel when several channels are input

        % Fit fooof
        fm = fit(fm);
        
        % S3. Average PSD / exponents
        switch(s.average_psd)
            case 'no'
                % Average exponents
                
        end
    end
    
    % Extract X for all participants

end