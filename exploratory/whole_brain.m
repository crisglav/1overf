% Estimate the aperiodic component 100 ROIs of all recordings
% 
% Pallotti Flaminia, Cristina Gil, TUM, 30.10.2023

clear all,
close all;

%% Settings

% Add fieldtrip and analysis functions
addpath('/rechenmagd4/toolboxes_and_functions/fieldtrip');
ft_defaults;
addpath('../analysis_functions');
addpath('../fooof_matlab');
addpath('/rechenmagd4/toolboxes_and_functions/plotting_functions');

% Load preprocessing parameters
load('../../results/features/params.mat');

% Create output folders
expPath = '../../results/features/fooof_matlab/wholebrain';
if ~exist(expPath)
    mkdir(expPath)
end
figures_path = '../../results/figures';

%% Load data
participants = readtable(fullfile(params.RawDataPath,'participants_clean.tsv'),'Filetype','text');
participant_id = participants.participant_id;
nSubj = height(participants);
nROI = length(params.parcellation.ROI);

try
    load(fullfile(expPath, 'patients_exponent_values.mat'));
    load(fullfile(expPath, 'healthy_exponent_values.mat'));
    load(fullfile(expPath, 'patients_offset_values.mat'));
    load(fullfile(expPath, 'healthy_offset_values.mat'));

catch
    % Initialize arrays to store exp values and offset values for patients and healthy controls
    patients_exponents = NaN(nSubj, nROI);
    healthy_exponents = NaN(nSubj, nROI);
    patients_offsets = NaN(nSubj, nROI);
    healthy_offsets = NaN(nSubj, nROI);

    for iSubj=1:nSubj

        bidsID = participant_id{iSubj};
        bidsID = [bidsID '_task-closed'];

        % Get info about this subject (patient or hc)
        group = participants.group(iSubj);

        % Loading power data from Cristina (updated version)
        power_struct = load(['../../results/features/power/' bidsID '_power.mat']);
        power = power_struct.power;
        freq = power.freq;


        exponent = nan(1,nROI);
        offset = nan(1,nROI);

        % Looping over each ROI and processing data for each ROI:
        for iROI = 1:nROI

            % Extract power at this ROI
            pow = power.powspctrm(iROI,:);

            % Fit a fooof object in the 2 - 40 Hz freq range for this ROI
            fm = fooof('freq_range',[2,40],'peak_width_limits',[1 12],'aperiodic_mode','fixed');
            fm = add_data(fm,freq,pow);
            fm = fit(fm);

            % Extracting information about the aperiodic component
            exponent(iROI) = fm.aperiodic_params(2); % accessing the second value (exponent) from the aperiodic_params field of fm
            offset(iROI) = fm.aperiodic_params(1);
        end

        % Assigning the exponent and offset for this ROI to the respective array according to the group
        if strcmp(group, 'pa')
            patients_exponents(iSubj, :) = exponent;
            patients_offsets(iSubj, :) = offset;

        elseif strcmp(group, 'hc')
            healthy_exponents(iSubj, :) = exponent;
            healthy_offsets(iSubj, :) = offset;
        end

    end

    % Handling NaN values and saving the results:
    pa_idx = categorical(participants.group) == 'pa';
    hc_idx = categorical(participants.group) == 'hc';

    patients_exponent_values = patients_exponents(pa_idx, :);
    save(fullfile(expPath, 'patients_exponent_values.mat'),'patients_exponent_values');
    healthy_exponent_values = healthy_exponents(hc_idx, :);
    save(fullfile(expPath, 'healthy_exponent_values.mat'),'healthy_exponent_values');
    patients_offset_values = patients_offsets(pa_idx, :);
    save(fullfile(expPath, 'patients_offset_values.mat'),'patients_offset_values');
    healthy_offset_values = healthy_offsets(hc_idx, :);
    save(fullfile(expPath, 'healthy_offset_values.mat'),'healthy_offset_values');
end

%% Plots

% Calculating mean exponent for each ROI, within each group:
patient_mean_exponents = mean(patients_exponent_values);
healthy_mean_exponents = mean(healthy_exponent_values);
exponents_difference = patient_mean_exponents - healthy_mean_exponents;

% Calculating mean offset for each ROI, within each group:
patient_mean_offset = mean(patients_offset_values);
healthy_mean_offset = mean(healthy_offset_values);
offsets_difference = patient_mean_offset - healthy_mean_offset;

 
% Figure 1: Exponents mean values Patients vs. Healthy Controls 
% Figure 2: Difference exponents over 100 ROIs between patients and healthy subjects

[exponents_PvsHC_fig, exponents_difference_fig] = wholebrain_plot(params,patient_mean_exponents, healthy_mean_exponents, exponents_difference); 
% sgtitle(exponents_PvsHC_fig,'Average aperiodic exponents');
title(exponents_difference_fig.Children(end), 'Aperiodic exponents PA-HC');
saveas(exponents_PvsHC_fig,fullfile(figures_path,'wholebrain_exponent_PA_HC.svg'));
saveas(exponents_difference_fig,fullfile(figures_path,'wholebrain_exponent_diff.svg'));

% Figure 3: Offsets mean values Patients vs. Healthy Controls 
% Figure 4: Difference offsets over 100 ROIs between patients and healthy subjects

[offsets_PvsHC_fig, offsets_difference_fig] = wholebrain_plot(params,patient_mean_offset, healthy_mean_offset, offsets_difference); 
% sgtitle(offsets_PvsHC_fig, 'Average aperiodic offsets');
title(exponents_difference_fig.Children(end), 'Aperiodic offsets PA-HC');
saveas(offsets_PvsHC_fig,fullfile(figures_path,'wholebrain_offset_PA_HC.svg'));
saveas(offsets_difference_fig,fullfile(figures_path,'wholebrain_offset_diff.svg'));