% Subgroup analysis
%
% Cristina Gil, TUM, 22.01.2024

% Load data
clear all,
close all;
%% Settings

% Add fieldtrip and analysis functions
addpath('/rechenmagd4/toolboxes_and_functions/fieldtrip');
ft_defaults;
% addpath('../analysis_functions');
% Load preprocessing parameters
load('../../results/features/params.mat');

%% Load the data
% Load the tsv in which the participant labels are randomized to be
% completely agnostic to the groups.
participants = readtable(fullfile(params.RawDataPath,'participants_clean.tsv'),'Filetype','text');
participants.diagnosis = categorical(participants.diagnosis);
participants.group = categorical(participants.group);

healthy = participants(participants.group == 'hc',:);
patients = participants(participants.group == 'pa',:);
summary(patients)

% The two biggest groups
cwp = participants(participants.diagnosis == 'CWP',:);
cbp = participants(participants.diagnosis == 'CBP',:);