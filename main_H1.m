% Scritpt to test the firt hypothesis
% 
% Cristina Gil Avila, TUM, 09.10.2023

clear all,
close all;
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
% params.SourcePath = '../results/source';
% if ~exist(params.SourcePath)
%     mkdir(params.SourcePath)
% end
params.PowerPFCPath = '../results/power/PFC';
if ~exist(params.PowerPFCPath)
    mkdir(params.PowerPFCPath)
end
params.FOOOFPath = '../results/fooof_matlab/PFC';
if ~exist(params.FOOOFPath)
    mkdir(params.FOOOFPath)
end

%%
% Estimate power spectrum
estimate_power(params);


% FOOOF 
compute_fooof_matlab(params);

