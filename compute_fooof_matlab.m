% Script testing the reimplementation of fooof in Matlab
%
% Original python code by Tom Donoghue et. al in https://github.com/fooof-tools/fooof
%
% Cristina Gil Avila, TUM, 01.09.2023
close all
clear all

% Add matlab FOOOF functions
addpath('fooof_matlab');

% FOOOF settings
self.freq_range        = [2 40];
self.peak_width_limits = [0.5 12]; % [1 8];
self.max_n_peaks       = inf; % 3;    
self.min_peak_height   = 0; % 0.15; 
self.peak_threshold    = 2;
% Fixed settings FOOOF settings
self.ap_percentile_thresh = 0.025;

% Input files
power_path = '../results/power/PFC/';
power_path = '../results/power/whole_brain/';


% Output files
fooof_path = '../results/fooof_matlab/whole_brain';
%% Load power data and fooof it

power_files = dir(fullfile(power_path,'*.mat'));

for iFile =1:length(power_files)
    file = power_files(iFile);
    load(fullfile(file.folder,file.name));

    % Average across channels
    avgpow = mean(pow,1);
    pow_orig = avgpow;
    freq_orig = freq;
    
    % Select a subset of frequencies defined in freq_range and transform
    % power to log space
    f = and(freq >= self.freq_range(1),freq<=self.freq_range(2));
    self.freq = freq_orig(f);
    self.pow = log10(pow_orig(f));

    % Fit fooof
    self = fit_fooof(self);
    
    % Save model data
    fname = [file.name(1:end-4) '_fooofm.mat'];
    save(fullfile(fooof_path,fname),'self');
    
end

