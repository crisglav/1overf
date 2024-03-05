% Compare offset between patients with chronic pain and healthy
% participants (visualization only)
%
% Cristina Gil Avila, TUM, 26.02.2024

clear all, close all

% Add plotting functions
addpath('../../toolboxes/raincloudplots');
addpath('../fooof_matlab'); % This is needed to read fooof objects

% Paths
participants_path = '../../data/blinded/';
fooof_path = '../../results/features/fooof_matlab/';
stats_path = '../../results/statistics';
figures_path = '../../results/figures';
if ~exist(stats_path,'dir')
    mkdir(stats_path)
end
if ~exist(figures_path,'dir')
    mkdir(figures_path)
end

%% Load demographic data 
participants_file = 'participants_clean.tsv';
% Participants.tsv file
participants = readtable(fullfile(participants_path,participants_file),'Filetype','Text');
participants.group = categorical(participants.group);

% Order the participants table in descending order by bidsID
participants_id = cellfun(@(x) str2double(x(5:7)),participants.participant_id,'UniformOutput',false);
participants_id = cell2mat(participants_id);
[participants_id_sorted,ix] = sort(participants_id);
participants_sorted = participants(ix,:);

% Separate patients and healthy participants groups
pa_mask = participants_sorted.group == 'pa';
hc_mask = participants_sorted.group == 'hc';

% Load patients data
patients = readtable(fullfile(participants_path,'patients_clean.tsv'),'Filetype','Text');

% Order the patients table in descending order by bidsID
patients_id = cellfun(@(x) str2double(x(5:7)),patients.participant_id,'UniformOutput',false);
patients_id = cell2mat(patients_id);
[~,ix] = sort(patients_id);
patients_sorted = patients(ix,:);

%% Load fooof data
fooof_files = dir(fullfile(fooof_path,'PFC','*.mat'));
n = length(fooof_files);

offset = nan(1,n);

% Load fooof parameters
for iSubj=1:n
    load(fullfile(fooof_files(iSubj).folder,fooof_files(iSubj).name));
    offset(iSubj) = fm.aperiodic_params(1);
end


%% Hypothesis 1: Do aperiodic offsets in the mPFC differ between patients and healthy participants?

% For visualization here only, regress out age from the aperiodic offsets
age = participants.age;
model_offset = fitlm(age, offset);
residuals_offset = model_offset.Residuals.Raw;

% Separate residuals per group
res_offset_hc = residuals_offset(hc_mask);
res_offset_pa = residuals_offset(pa_mask);

% Plot aperiodic offsets for patients and healthy participants
f1 = figure('Units','centimeters','Position',[25 25 22 9]);
tiledlayout(1,2);
ax = nexttile;
h1 = raincloud_plot_vertical(res_offset_hc,'box_dodge',1, 'color', '#69A5C4', 'alpha', 1, 'bxcl', [.2 .2 .2], 'line_width', 1.5, 'box_dodge_amount', 0.5,'dot_dodge_amount',0.5,'wdth',0.3);
h2 = raincloud_plot_vertical(res_offset_pa, 'box_dodge', 1, 'color', "#D98F4A", 'alpha',1, 'bxcl', [.2 .2 .2], 'line_width', 1.5, 'box_dodge_amount', 1,'dot_dodge_amount',1,'wdth',0.3);
xlim([-3, 1.5])
legend([h1{1} h2{1}], {'HC', 'PA'},'Location','southwest');
title('')
set(ax,'XTick',[],'XTickLabel',[]);
ylabel('Aperiodic offsets residuals');
box off
grid on

% Save aperiodic offsets in a .csv file
participants_sorted.offset_PFC = offset';
participants_sorted.group_binary = pa_mask;
writetable(participants_sorted, fullfile(stats_path,'e2_offset_h1.csv'));

% Save figure
% saveas(f1,fullfile(figures_path,'e2_offset_h1.svg'));

%% Hypothesis 2: Do aperiodic offsets in the mPFC correlate with pain ratings in patients?
% Get avg pain ratings and age and discard the patient from which we don't have a pain rating
pain_mask = ~isnan(patients_sorted.avg_pain);
avg_pain = patients_sorted.avg_pain(pain_mask);
age_pa = patients_sorted.age(pain_mask);
offset_pa_original = offset(pa_mask);
offset_pa = offset_pa_original(pain_mask);

% For visualization here only, regress out age from the pain ratings
model_pain = fitlm(age_pa, avg_pain);
residuals_pain = model_pain.Residuals.Raw;

% For visualization here only, regress out age from the aperiodic offsets
model_offset = fitlm(age_pa, offset_pa);
residuals_offset = model_offset.Residuals.Raw;

% Scatter plot of residuals
model_res = fitlm(model_pain.Residuals.Raw, model_offset.Residuals.Raw);
% f2 = figure;
nexttile;
h = plot(model_res);
h(1).Marker = 'o';
h(1).MarkerFaceColor = '#D98F4A';
h(1).MarkerEdgeColor = 'none';
h(2).Color = [0.5 0.5 0.5];
h(2).LineWidth = 1.5;
h(3).Color = [0.5 0.5 0.5];
h(3).LineWidth = 1.5;
h(4).Color = [0.5 0.5 0.5];
h(4).LineWidth = 1.5;
xlabel('Pain residuals');
ylabel('Aperiodic offset residuals');
title('');
box off;
legend off;
grid on

% Save aperiodic offsets in a csv file
patients_sorted.offset_PFC = offset_pa_original';
writetable(patients_sorted, fullfile(stats_path,'e2_offset_h2.csv'));

% Save figure
% saveas(f2,fullfile(figures_path,'e2_offset_h2.svg'));
saveas(f1,fullfile(figures_path,'e2_offset.svg'));
