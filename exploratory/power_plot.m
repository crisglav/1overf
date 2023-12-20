% Script to plot the power spectrum from subjects with negative offset values 
% checking if negative values of the offset actually make sense 

% Flaminia Pallotti, TUM, 24.11.2023
clear all,
close all;
%%
% Testing for sub-197 (healthy control) at ROI 15
pow_id = 196;
iROI = 15;
power_path = '../../results/features/power';
power_files = dir(fullfile(power_path,'*.mat'));
file = power_files(pow_id);
load(fullfile(file.folder,file.name));

%%
% Extract power at this ROI and average across virtual channels
cfg = [];
cfg.channel = iROI;
power_subj = ft_selectdata(cfg,power);
freq = power_subj.freq;
pow = power_subj.powspctrm;  
avgpow = mean(pow,1);

%% Fit fooof model
% Initialize fooof object with some default settings
fm = fooof('freq_range',[2 40],'aperiodic_mode','knee');
% Add data to the model in the freq range 2 - 40 Hz (pre-specified)
fm = add_data(fm,freq,avgpow);
fm = fit(fm);
fm.get_results()
fm.plot('fig_save',true,'file_name',file.name(1))
title('Power Spectrum sub-197, ROI 15 - Logarithmic scale')
fm.plot_noLog('fig_save',true,'file_name',file.name(1))
title('Power Spectrum sub-197, ROI 15')