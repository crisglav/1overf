% Estimate the aperiodic component in the PFC of all recordings
% 
% Cristina Gil Avila, TUM, 09.10.2023

clear all,
close all;
%% Settings
% Define the number of cores for parallelization
% params.Ncores = 2;
% if(isempty(gcp('nocreate')))
%     parObj = parpool(params.Ncores);
% end

% Add fieldtrip and analysis functions
addpath('/rechenmagd4/toolboxes_and_functions/fieldtrip');
ft_defaults;
addpath('analysis_functions');
addpath('fooof_matlab');

% Load preprocessing parameters
fparams = '../data/blinded/derivatives_v2023_08_18/params.json';
params = load_params(fparams);

% Define the atlas for source localization, create ROIs mask and precompute
% the source model
params.AtlasPath = '../toolboxes/schaefer_parcellations/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_1mm.Centroid_RAS.csv';
params = create_parcellation(params);

% Define the frequency band in which to compute the spatial filter
params.FreqBand.fullSpectrum = [0.5 100.5];

% Create output folders
params.SourcePath = '../results/features/source';
if ~exist(params.SourcePath)
    mkdir(params.SourcePath)
end
params.VdataPath = '../results/features/vdata';
if ~exist(params.VdataPath)
    mkdir(params.VdataPath)
end
params.PowerPath = '../results/features/power';
if ~exist(params.PowerPath)
    mkdir(params.PowerPath)
end
params.FOOOFPath = '../results/features/fooof_matlab/PFC';
if ~exist(params.FOOOFPath)
    mkdir(params.FOOOFPath)
end

save('../results/features/params.mat','params');
%%
participants = readtable(fullfile(params.RawDataPath,'participants_rand.tsv'),'Filetype','text');
participant_id = participants.participant_id;
n = height(participants);

for iSubj=1:n
    
    bidsID = participant_id{iSubj};
    bidsID = [bidsID '_task-closed'];
    try
        % Load EEG preprocessed data
        data = load_preprocessed_data(params,bidsID);
    catch
        continue
    end
    % Compute source reconstruction
    source = compute_spatial_filter(params,data,'fullSpectrum');
    parsave(fullfile(params.SourcePath,[bidsID '_source.mat']),source);

    % Extract virtual channel data
    cfg = [];
    cfg.parcellation = 'ROI';
    vdata = ft_virtualchannel(cfg,data,source,params.parcellation);
    parsave(fullfile(params.VdataPath, [bidsID '_vdata.mat']),vdata);

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
    power = ft_freqanalysis(cfg, vdata);
    parsave(fullfile(params.PowerPath, [bidsID '_power.mat']),power);
    
    % Extract power at the PFC and average across ROI virtual channels
    cfg = [];
    cfg.channel = find(params.PFC_mask);
    power_PFC = ft_selectdata(cfg,power);
    freq = power_PFC.freq;
    pow = power_PFC.powspctrm;  
    avgpow = mean(pow,1);
    
    % Fit a fooof object in the 2 - 40 Hz freq range 
    fm = fooof('freq_range',[2,40],'aperiodic_mode','fixed');   
    fm = add_data(fm,freq,avgpow);
    fm = fit(fm);
    fname = [bidsID '_PFC_fooofm.mat'];
    parsave(fullfile(params.FOOOFPath,fname),fm);
    
end
