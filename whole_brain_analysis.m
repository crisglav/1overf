% Exploratory analysis of the aperiodic exponent at the whole brain level.
%
% Cristina Gil Avila, Technical University of Munich, 15.09.2023

clear all, close all

% Add statisical functions and plotting functions
run('bayes_factor/installBayesFactor.m')
addpath('raincloudplots');
% Paths
results_path = '../results/jasp';
fooof_path = '../results/fooof/';
participants_path = '../data/blinded/';

% Analsysis type:
analysis = 'real'; % 'blinded' / 'real'

switch analysis
    case 'blinded'
        participants_file = 'participants_rand.tsv';
    case 'real'
        participants_file = 'participants_clean.tsv';
end
% Participants.tsv file
participants = readtable(fullfile(participants_path,participants_file),'Filetype','Text');
participants.group = categorical(participants.group);

% Order the participants.tsv in descending order by bidsID
id = cellfun(@(x) str2double(x(5:7)),participants.participant_id,'UniformOutput',false);
id = cell2mat(id);
[~,ix] = sort(id);
participants_sorted = participants(ix,:);

% Separate patients and healthy participants groups
pa_mask = participants_sorted.group == 'pa';
hc_mask = participants_sorted.group == 'hc';

switch analysis
    case 'blinded'
        patients = readtable(fullfile(participants_path,'patients_rand.tsv'),'Filetype','Text');
    case 'real'
        patients = participants(participants.group == 'pa',:);
        % Replace the missing values of patients without avg_pain by their current
        % pain (curr_pain)
        m = ismissing(patients.avg_pain);
        patients.avg_pain(m) = cellfun(@str2double, patients.curr_pain(m));
end
% Regions of interest
rois = {'PFC','S1','VIS'};
%% Load fooof data

for iRoi=1:length(rois)
    
    roi = rois{iRoi};
    roi_path = fullfile(fooof_path,roi);
    
    roi_files = dir(fullfile(roi_path,'*.mat'));
    n = length(roi_files);

    er = nan(1,n);
    r2 = nan(1,n);
    ex = nan(1,n);

    % Load fooof parameters
    for iSubj=1:n
        load(fullfile(roi_files(iSubj).folder,roi_files(iSubj).name));

        er(iSubj) = error;
        r2(iSubj) = r_squared;
        ex(iSubj) = aperiodic_params(2);

        clear error r_squared aperiodic_params gaussian_params peak_params 
    end
    
    errors.(roi) = er;
    r_squareds.(roi) = r2;
    exponents.(roi) = ex;
    
end

% % Plot model errors and r_squared
% figure;
% tlc = tiledlayout(1,2);
% nexttile
% 
% scatter(1:n,er,'filled');
% xlabel('Recording');
% ylabel('error');
% title('Model error');
% 
% nexttile
% scatter(1:n,r2,'filled');
% xlabel('Recording');
% ylabel('R^2');
% title('R^2');

%% Analysis

%% Hypothesis 1: differencens in PFC aperiodic exponents between patients and healthy participants
% Plot aperiodic exponents
cb = lines(2);
f1 = figure();
ax = gca;
h1 = raincloud_plot(exponents.PFC(hc_mask), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', 0.15,'box_col_match', 0, ...
     'cloud_edge_col', cb(1,:),'plot_top_to_bottom',1);
h2 = raincloud_plot(exponents.PFC(pa_mask), 'box_on', 1, 'color', cb(2,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,...
     'cloud_edge_col', cb(2,:),'plot_top_to_bottom',1);
 
legend([h1{1} h2{1}], {'HC', 'PA'});
title('Aperiodic exponents')
set(ax,'XLim', [0 2], 'YLim', [-1.5 2.5]);
set(ax,'YTickLabel',[]);
xlabel('Exponents');
box off

% Hypothesis 1: Bayesian t-test on the aperiodic exponents
[bf10,p] = bf.ttest2(exponents.PFC(hc_mask),exponents.PFC(pa_mask));

str = {sprintf('BF = %0.3f',bf10),sprintf('p = %0.3f',p)};
annotation('textbox',[0.15 0.6 0.3 0.3],'String',str,'FitBoxToText','on');

% Check age and gender covariates ANCOVA in JASP
participants_sorted.exp_PFC = exponents.PFC';
writetable(participants_sorted, fullfile(results_path,['exp_PFC_' analysis '.csv']));

% Save figure
exportgraphics(f1,fullfile(results_path,'hypothesis1.jpg'));

