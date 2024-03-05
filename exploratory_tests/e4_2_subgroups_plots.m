% Compare aperiodic exponents between healthy participants and the two
% biggests subgroups of patients (Chronic back pain patients = CBP) and
% chronic widespread pain patients (CWP). Generate raincloudplots
%
% Cristina Gil, Flaminia Palloti, TUM, 27.02.2024


clear all, close all,

% Add raincloudplots
addpath('../../toolboxes/raincloudplots');

% Path for the figures
figures_path = '../../results/figures';

% Load prepare data from csv, were we matched the HC in age, gender and eeg setting to the patient subgroup.
% First you have to run e4_1_subgroups_statistics.R
datapath = '../../results/statistics';
cwp = readtable(fullfile(datapath,'e4_cwp.csv'));
cbp = readtable(fullfile(datapath,'e4_cbp.csv'));
cwp.group = categorical(cwp.group);
cbp.group = categorical(cbp.group);

% Separate residuals per group
cwp_hc_mask = cwp.group == 'hc';
cwp_pa_mask = cwp.group == 'pa';
apexp_cwp_hc = cwp.res_apexp(cwp_hc_mask);
apexp_cwp_pa = cwp.res_apexp(cwp_pa_mask);

cbp_hc_mask = cbp.group == 'hc';
cbp_pa_mask = cbp.group == 'pa';
apexp_cbp_hc = cbp.res_apexp(cbp_hc_mask);
apexp_cbp_pa = cbp.res_apexp(cbp_pa_mask);


%% Hypothesis 1: Do aperiodic exponents in the mPFC differ between subgroups of patients and healthy participants?

% Plot aperiodic exponents for patients and healthy participants
f1 = figure('Units','centimeters','Position',[25 25 22 20]);
tiledlayout(2,2);
ax=nexttile;
h1 = raincloud_plot_vertical(apexp_cbp_hc,'box_dodge',1, 'color', '#69A5C4', 'alpha', 1, 'bxcl', [.2 .2 .2], 'line_width', 1.5, 'box_dodge_amount', 0.5,'dot_dodge_amount',0.5,'wdth',0.3);
h2 = raincloud_plot_vertical(apexp_cbp_pa, 'box_dodge', 1, 'color', "#D95D4A", 'alpha',1, 'bxcl', [.2 .2 .2], 'line_width', 1.5, 'box_dodge_amount', 1,'dot_dodge_amount',1,'wdth',0.3);
xlim([-3, 1.5])
ylim([-0.8 0.8])
legend([h1{1} h2{1}], {'matched HC', 'CBP'},'Location','southwest');
% title('CBP patients')
set(ax,'XTick',[],'XTickLabel',[]);
ylabel('Aperiodic exponents residuals');
box off
grid on

ax = nexttile(3);
% xlim([-3, 1.5])
h1 = raincloud_plot_vertical(apexp_cwp_hc,'box_dodge',1, 'color', '#69A5C4', 'alpha', 1, 'bxcl', [.2 .2 .2], 'line_width', 1.5, 'box_dodge_amount', 0.5,'dot_dodge_amount',0.5,'wdth',0.3);
h2 = raincloud_plot_vertical(apexp_cwp_pa, 'box_dodge', 1, 'color', "#D9B64A", 'alpha',1, 'bxcl', [.2 .2 .2], 'line_width', 1.5, 'box_dodge_amount', 1,'dot_dodge_amount',1,'wdth',0.3);
xlim([-3, 1.5])
ylim([-0.8 0.8])
legend([h1{1} h2{1}], {'matched HC', 'CWP'},'Location','southwest');
% title('CWP patients')
set(ax,'XTick',[],'XTickLabel',[]);
ylabel('Aperiodic exponents residuals');
box off
grid on
% saveas(f1,fullfile(figures_path,'e4_subgroups_h1.svg'));



%% Hypothesis 2: Do aperiodic exponents in the mPFC correlate with pain intensity in subgroups of patients?
% Plot aperiodic exponents for patients and healthy participants
% f2 = figure('Units','centimeters','Position',[25 25 22 9]);
% tlc = tiledlayout(1,2);
ax=nexttile;
model_cbp = fitlm(cbp.res_pain(cbp_pa_mask), cbp.res_apexp(cbp_pa_mask));
h = plot(model_cbp);
h(1).Marker = 'o';
h(1).MarkerFaceColor = '#D95D4A';
h(1).MarkerEdgeColor = 'none';
h(2).Color = [0.5 0.5 0.5];
h(2).LineWidth = 1.5;
h(3).Color = [0.5 0.5 0.5];
h(3).LineWidth = 1.5;
h(4).Color = [0.5 0.5 0.5];
h(4).LineWidth = 1.5;
xlabel('Pain residuals');
ylabel('Aperiodic exponent residuals');
title('');
box off;
legend off;
grid on

model_cwp = fitlm(cwp.res_pain(cwp_pa_mask), cwp.res_apexp(cwp_pa_mask));
nexttile;
h = plot(model_cwp);
h(1).Marker = 'o';
h(1).MarkerFaceColor = '#D9B64A';
h(1).MarkerEdgeColor = 'none';
h(2).Color = [0.5 0.5 0.5];
h(2).LineWidth = 1.5;
h(3).Color = [0.5 0.5 0.5];
h(3).LineWidth = 1.5;
h(4).Color = [0.5 0.5 0.5];
h(4).LineWidth = 1.5;
xlabel('Pain residuals');
ylabel('Aperiodic exponent residuals');
title('');
box off;
legend off;
ylim([- 0.6 0.6])
grid on

% saveas(f2,fullfile(figures_path,'e4_subgroups_h2.svg'));
saveas(f1,fullfile(figures_path,'e4_subgroups.svg'))