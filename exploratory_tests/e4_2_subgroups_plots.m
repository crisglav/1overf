Compare aperiodic exponents between healthy participants and the two
biggests subgroups of patients (Chronic back pain patients = CBP) and
chronic widespread pain patients (CWP)
%
% Cristina Gil, Flaminia Palloti, TUM, 27.02.2024

% Load data from aperiodic exponents in the mPFC

clear all, close all,

datapath = '../../results/statistics/exp_PFC_real.csv';
data = readtable(datapath);
data.group = categorical(data.group);
data.diagnosis = categorical(data.diagnosis);

% Separate participants by group and diagnosis
hc_mask = data.group == 'hc';
cbp_mask = data.diagnosis == 'CBP';
cwp_mask = data.diagnosis == 'CWP';