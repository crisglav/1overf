% This scripts compares two different preprocessing strategies.
% 
% A) DISCOVER-EEG pipeline v1.0.0
% This corresponds to
% /rechenmagd3/Experiments/2023_1overf/data/blinded/code_discover_eeg
% branch master
%
% B) DISCOVER-EEG pipeline with a highpass filter before cleanline
% This corresponds to
% /rechenmagd3/Experiments/2023_1overf/data/blinded/code_discover_eeg on
% branch 1overf
% This option was considered to completely eliminate line noise, as this
% could affect fooof estimation in the range 40 - 60 Hz.
% 
% Cristina Gil Avila. TUM. 09.10.2023

clear all,
close all,
addpath /rechenmagd4/toolboxes_and_functions/fieldtrip
ft_defaults
addpath(genpath('/rechenmagd3/Experiments/2023_1overf/data/blinded/code_discover_eeg/external_functions')); % std shade

% Preprocessed data of option A: '../data/blinded/derivatives_v2023_08_18/params.json'
% Preprocessed data of option B: '../data/blinded/derivatives_v2023_08_18/params.json'

% Paths to the power densities
path_a = '../../results/power/PFC_original';
path_b = '../../results/power/PFC_v2023_09_28';

% load all the files
files_a = dir(fullfile(path_a,'*.mat'));
power_a = cell(1,length(files_a));

for i = 1:length(files_a)
    power = load(fullfile(files_a(i).folder,files_a(i).name));
    % Average across channels
    power_a{i} = mean(power.pow,1);
end
a = reshape(cell2mat(power_a),[length(power_a{1}),length(files_a)]);
a = a';

files_b = dir(fullfile(path_b,'*.mat'));
power_b = cell(1,length(files_b));

for i = 1:length(files_b)
    power = load(fullfile(files_a(i).folder,files_b(i).name));
    % Average across channels
    power_b{i} = mean(power.pow,1);
end
b = reshape(cell2mat(power_b),[length(power_b{1}),length(files_b)]);
b = b';

freq = power.freq;

figure;
stdshade(log(a),0.1,[0 0.447 0.741],freq);
hold on
stdshade(log(b),0.1,[0.455 0.674 0.188],freq);

figure;
plot(power.freq,log(b));

