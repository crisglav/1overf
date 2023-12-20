%% Whole brain - statistical analysis
clear all,
close all;

%%
% Add statisical functions and plotting functions
run('../../toolboxes/bayes_factor/installBayesFactor.m')
addpath('....//toolboxes/raincloudplots');
addpath('fooof_matlab');
addpath('/rechenmagd4/toolboxes_and_functions/fieldtrip');
ft_defaults;

% Paths
participants_path = '../../data/blinded/';

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

% Loading data 
load('../../results/features/fooof_matlab/wholebrain/patients_exponent_values.mat');
load('../../results/features/fooof_matlab/wholebrain/healthy_exponent_values.mat');
load('../../results/features/fooof_matlab/wholebrain/patients_offset_values.mat');
load('../../results/features/fooof_matlab/wholebrain/healthy_offset_values.mat');

%% Hypothesis 1: changes in the aperiodic exponent in any ROI between patients and healthy participants
n_ROIs = 100;
alpha = 0.05;

% Adjusting alpha level (Bonferroni correction), taking into account the multiple comparison problem 
adjusted_alpha = alpha / n_ROIs;

% Inizialize arrays to store results
bf10_results_exp = zeros(1, n_ROIs);
p_values_exp = zeros(1, n_ROIs);
bf10_results_off = zeros(1, n_ROIs);
p_values_off = zeros(1, n_ROIs);

% Iterating over 100 ROIs, performing Bayesian Independent Samples T-test

for iROI = 1:100
    
    [bf10_exp, p_exp] = bf.ttest2(patients_exponent_values(:, iROI), healthy_exponent_values(:, iROI));
    [bf10_off, p_off] = bf.ttest2(patients_offset_values(:, iROI), healthy_offset_values(:, iROI));
    
    % Storing the results for the current ROI:
    bf10_results_exp(iROI)= bf10_exp;
    p_values_exp(iROI) = p_exp; % p_values contains one p_value for each ROI
    
    bf10_results_off(iROI)= bf10_off;
    p_values_off(iROI) = p_off; % p_values contains one p_value for each ROI

end

% Set significant ROIs with BF > 10 or BF < 0.1

significant_ROIs_BF_exp = find(bf10_results_exp > 10 | bf10_results_exp < 0.1); % exponents
significant_ROIs_BF_off = find(bf10_results_off > 10 | bf10_results_off < 0.1); % offsets

%% False Discovery Rate (FDR)
% Traditional Bonferroni method is too conservative
% With FDR you identify as many features as possible while incurring a relatively low proportion of false positives

fdr_exp = mafdr(p_values_exp); % calculating adjusted pFDR values
significant_ROIs_FDR_exp = find (fdr_exp < alpha);

fdr_off = mafdr(p_values_off); % calculating adjusted pFDR values
significant_ROIs_FDR_off = find (fdr_off < alpha);

%% Cluster analysis and cluster based permutation testing (FieldTrip)

% Data preparation for ft_freqstatistics (FieldTrip structure):

data_patient = struct();
data_patient.dimord = 'subj_chan_freq';
data_patient.label = cellstr(num2str((1:100)'));
data_patient.trialinfo = ones(size(patients_exponent_values, 1), 1); % assigning group label
data_patient.freq = 1:0.5:100;
data_patient.exp = patients_exponent_values; % data about the exponents
% data_patient.offsets = patients_offset_values; % data about the offsets

data_healthy = struct();
data_healthy.dimord = 'subj_chan_freq';
data_healthy.label = cellstr(num2str((1:100)')); % labeling each ROI
data_healthy.trialinfo = 2*ones(size(healthy_exponent_values, 1), 1);
data_healthy.freq = 1:0.5:100;
data_healthy.exp = healthy_exponent_values; % data about the exponents 


%% Neighbours configuration (ft_prepare_neighbours)

% Load atlas positions
params.AtlasPath = '../../toolboxes/schaefer_parcellations/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_1mm.Centroid_RAS.csv';
atlas = readtable(params.AtlasPath);
roi_coordinates = [atlas.R, atlas.A, atlas.S];
% roi_name = atlas.ROIName;

% % Defining a layout with ft_prepare_layout
% layout = struct();
% layout.label = cellstr(num2str((1:100)'));
% layout.pos = [atlas.R, atlas.A, atlas.S];
% 
cfg = [];
cfg.layout = 'eeg1010.mat';
cfg.feedback = 'yes';
[layout, cfg] = ft_prepare_layout(cfg);
save('layout.mat', layout);

cfg=[];
cfg.method = 'triangulation';
cfg.layout = 'layout.mat';
neighbours_struct = ft_prepare_neighbours(cfg);

% %% Creating a grid
% grid_spacing = 30; %mm
% grid_dim = ceil((max(roi_coordinates) - min(roi_coordinates)) / grid_spacing); 
% 
% [x, y, z] = meshgrid(...
%     linspace(min(roi_coordinates(:, 1)), max(roi_coordinates(:, 1)), grid_dim(1)), ...
%     linspace(min(roi_coordinates(:, 2)), max(roi_coordinates(:, 2)), grid_dim(2)), ...
%     linspace(min(roi_coordinates(:, 3)), max(roi_coordinates(:, 3)), grid_dim(3)));
% 
% % Reshaping the grid 
% grid_points = [x(:), y(:), z(:)];
% 
% % Visualize the grid 
% figure();
% scatter3(roi_coordinates(:, 1), roi_coordinates(:, 2), roi_coordinates(:, 3), 'filled', 'MarkerEdgeColor', 'k');
% hold on; 
% scatter3(grid_points(:, 1), grid_points(:, 2), grid_points(:, 3), 50, 'r', 'filled');
% grid on; 
% 
% % Find the closest grid point for each ROI
% nearest_indices = zeros(size(roi_coordinates, 1), size(grid_points, 1));
% n_ROI = 100;
% neighbours = struct('label', cell(1, n_ROI), 'neighblabel', cell(1, n_ROI));
% 
% for i=1:size(roi_coordinates, 1)
%     distances = sqrt(sum((grid_points - roi_coordinates(i, :)).^2, 2));
%     [~, sorted_indices] = sort(distances);
%     nearest_indices(i, :) = sorted_indices';
%    
%     for j=2:length(sorted_indices)
%         closest_index = sorted_indices(j);
%         closest_rois = find(nearest_indices (:, j) == closest_index);
%         closest_rois = setdiff(closest_rois, i);
% 
%         if ~isempty(closest_rois)  
%             break;
%         end
%     
%     end
%     
%     neighbours(i).label = num2str(i);
%     neighbours(i).neighblabel = cellfun(@num2str, num2cell(closest_rois), 'UniformOutput', false)';
%     
%     for j = 1: length(closest_rois)
%         idx = closest_rois (j);
% 
%         if ~ismember(num2str(i), neighbours(idx).neighblabel)
%             neighbours(idx).neighblabel{end+1} = num2str(i);
%         end
%     end
% end
% 
% save('neighbours_template.mat', 'neighbours')
% cfg=[];
% cfg.method = 'template';
% cfg.template = 'neighbours_template.mat';
% neighbours_struct = ft_prepare_neighbours(cfg);
%% 
% Define the FieldTrip configuration for ft_freqstatistics 
cfg = [];
cfg.channel = 'all';
cfg.avgoverchan = 'no';
cfg.frequency = 'all';
cfg.avgoverfreq = 'yes';
cfg.parameter = 'exp';
cfg.method = 'ft_statistics_montecarlo';
cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.cluster.statistics = 'maxsum';
cfg.clusterthreshold = 'nonparametric_common';
cfg.minnbchan = 2;
cfg.tail = 0;
cfg.alpha = 0.05;
cfg.correcttail = 'alpha';
cfg.computerprob = 'yes';
cfg.numrandomization = 1000;
cfg.neighbours = neighbours_struct; % previously created

design_matrix = zeros(1, size(patients_exponent_values, 1) + size(healthy_exponent_values,1));
design_matrix(1, 1:size(patients_exponent_values, 1)) = 1;
design_matrix(1, (size(patients_exponent_values, 1) + 1): end) = 2;
cfg.design = design_matrix;
cfg.ivar = 1; % 1 independent variable (group assignment)

% Performing statistical analysis:
stat2 = ft_freqstatistics (cfg, data_patient, data_healthy);

% %%
% % Plotting the results 
% cfg = [];
% cfg.alpha = stat2.cfg.alpha;
% cfg.parameter = 'stat';
% cfg.zlim = [-3 3];
% cfg.elec = elec;
% ft_clusterplot(cfg, stat2);
% 
% cfg=[];
% cfg.elec = elec;
% cfg.zlim = [1.5 3];
% cfg.xlim = [8 15];
% cfg.parameter = 'powspctrm_b';
% cfg.markersymbol = '.';
% cfg.comment = 'no';
% cfg.colormap = 'jet';
% cfg.colorbar = 'no';
% 
% figure('position', [680 240 1039 420]);
% subplot(1,2,1); ft_topoplotER(cfg, data_patient.exponents); colorbar; title('Patients - Exponents')
% subplot(1,2,2); ft_topoplotER(cfg, data_healthy.exponents); colorbar; title('Healthy Controls - Exponents')

%% GLOBAL APERIODIC EXPONENT
% Average aperiodic exponents across ROIS
global_pa = mean(patients_exponent_values,2);
global_hc = mean(healthy_exponent_values,2);

c = cell(264,1);
c(1:149) = {'pa'};
c(150:end) = {'hc'};

t = table([global_pa;global_hc],c,'VariableNames',{'Global_exp','Group'});
writetable(t,'../../results/statistics/global_exponent.csv');

% Raincloud plots
cb = lines(2);
f1 = figure();
ax = gca;
h1 = raincloud_plot(global_hc, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', 0.15,'box_col_match', 0, ...
     'cloud_edge_col', cb(1,:),'plot_top_to_bottom',1);
h2 = raincloud_plot(global_pa, 'box_on', 1, 'color', cb(2,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,...
     'cloud_edge_col', cb(2,:),'plot_top_to_bottom',1);
 
legend([h1{1} h2{1}], {'HC', 'PA'});
title('Aperiodic exponents')
set(ax,'XLim', [0 2], 'YLim', [-1.5 2.5]);
set(ax,'YTickLabel',[]);
xlabel('Exponents');
box off
exportgraphics(f1,fullfile(figures_path,'hypothesis1_globalexp.jpg'));
