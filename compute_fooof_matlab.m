% Script testing the reimplementation of fooof in Matlab
%
% Original python code by Tom Donoghue et. al in https://github.com/fooof-tools/fooof
%
% Cristina Gil Avila, TUM, 01.09.2023
close all
clear all

% Add matlab FOOOF functions
addpath('fooof_matlab');

% Input files
power_path = '../results/power/PFC/';

% Output files
fooof_path = '../results/fooof_matlab/PFC';

%% Load power data and fooof it
power_files = dir(fullfile(power_path,'*.mat'));

for iFile =1:length(power_files)
    file = power_files(iFile);
    load(fullfile(file.folder,file.name));

    % Average across channels
    avgpow = mean(pow,1);
    
    % Initialize a fooof object 
    fm = fooof();
    
    % Add the data in the frequency range 2 - 40 Hz
    fm = add_data(fm,freq,avgpow,[2,40]);
    
    % Fit the model
    fm = fit(fm);
       
    % Save model data
    fname = [file.name(1:end-4) '_fooofm.mat'];
    save(fullfile(fooof_path,fname),'fm');
    
end

