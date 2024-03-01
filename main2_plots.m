% Prepare the data for statistical analyisis and generate plots for the main analysis
% Statistical analysis are performed in R in main3_statistics.R
%
% Cristina Gil Avila, TUM, 12.02.2024

clear all, close all

% Analsysis type:
analysis = 'real'; % 'blinded' / 'real'

% Add plotting functions
addpath('../toolboxes/raincloudplots');
addpath('fooof_matlab'); % This is needed to read fooof objects

% Paths
participants_path = '../data/blinded/';
fooof_path = '../results/features/fooof_matlab/';
stats_path = '../results/statistics';
figures_path = '../results/figures';
if ~exist(stats_path,'dir')
    mkdir(stats_path)
end
if ~exist(figures_path,'dir')
    mkdir(figures_path)
end

%% Load demographic data (randomized or not)
switch analysis
    case 'blinded'
        participants_file = 'participants_rand.tsv';
    case 'real'
        participants_file = 'participants_clean.tsv';
end
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

% Load patients data with randomized or not pain ratings
switch analysis
    case 'blinded'
        patients = readtable(fullfile(participants_path,'patients_rand.tsv'),'Filetype','Text');
    case 'real'
        patients = readtable(fullfile(participants_path,'patients_clean.tsv'),'Filetype','Text');
end
% Order the patients table in descending order by bidsID
patients_id = cellfun(@(x) str2double(x(5:7)),patients.participant_id,'UniformOutput',false);
patients_id = cell2mat(patients_id);
[~,ix] = sort(patients_id);
patients_sorted = patients(ix,:);

%% Load fooof data
fooof_files = dir(fullfile(fooof_path,'PFC','*.mat'));
n = length(fooof_files);

errors = nan(1,n);
r_squareds = nan(1,n);
apexp = nan(1,n);

% Load fooof parameters
for iSubj=1:n
    load(fullfile(fooof_files(iSubj).folder,fooof_files(iSubj).name));
    
    errors(iSubj) = fm.error(1); % mae
    r_squareds(iSubj) = fm.r_squared;
    apexp(iSubj) = fm.aperiodic_params(end);
end

% % Plot fooof model errors 
% figure;
% tlc = tiledlayout(1,2);
% nexttile
% scatter(1:n,errors,'filled');
% xlabel('Recording');
% ylabel('error');
% title('Model error');
% nexttile
% scatter(1:n,r_squareds,'filled');
% xlabel('Recording');
% ylabel('R^2');
% title('R^2');

%% Hypothesis 1: Do aperiodic exponents in the mPFC differ between patients and healthy participants?

% For visualization here only, regress out age from the aperiodic exponents
age = participants_sorted.age;
model_apexp = fitlm(age, apexp);
residuals_apexp = model_apexp.Residuals.Raw;

% Separate residuals per group
res_apexp_hc = residuals_apexp(hc_mask);
res_apexp_pa = residuals_apexp(pa_mask);

% Plot aperiodic exponents for patients and healthy participants
f1 = figure('Units','centimeters','Position',[25 25 11 9]);
ax = gca;
h1 = raincloud_plot_vertical(res_apexp_hc,'box_dodge',1, 'color', '#69A5C4', 'alpha', 1, 'bxcl', [.2 .2 .2], 'line_width', 1.5, 'box_dodge_amount', 0.5,'dot_dodge_amount',0.5,'wdth',0.3);
h2 = raincloud_plot_vertical(res_apexp_pa, 'box_dodge', 1, 'color', "#D98F4A", 'alpha',1, 'bxcl', [.2 .2 .2], 'line_width', 1.5, 'box_dodge_amount', 1,'dot_dodge_amount',1,'wdth',0.3);
 
legend([h1{1} h2{1}], {'HC', 'PA'},'Location','southwest');
% title('H1')
set(ax,'XTick',[],'XTickLabel',[]);
ylim([-0.8, 0.8]);
xlim([-3, 1.5]);
ylabel('Aperiodic exponent residuals');
box off
grid on
% Save aperiodic exponents in a .csv file
participants_sorted.exp_PFC = apexp';
participants_sorted.group_binary = pa_mask;
writetable(participants_sorted, fullfile(stats_path,['exp_PFC_' analysis '.csv']));

% Save figure
saveas(f1,fullfile(figures_path,['hypothesis1_' analysis '.svg']));

%% Hypothesis 2: Do aperiodic exponents in the mPFC correlate with pain ratings in patients?
% To select the subset of aperiodic exponents belonging to the real
% patients we need to extract the mask of real patients.  
pa_mask_real = ismember(participants_id_sorted,patients_id);
apexp_pa_original = apexp(pa_mask_real);
% Check that you did it well
% all(ismember(participants_sorted(pa_mask_real,:).participant_id,patients_sorted.participant_id))

% Create a column for age in the patients_sorted file based on the
% participants_sorted table
patients_sorted.age = participants_sorted.age(pa_mask_real);

% Get avg pain ratings and age and discard the patient from which we don't have a pain rating
pain_mask = ~isnan(patients_sorted.avg_pain);
avg_pain = patients_sorted.avg_pain(pain_mask);
age_pa = patients_sorted.age(pain_mask);
apexp_pa = apexp_pa_original(pain_mask);

% For visualization here only, regress out age from the pain ratings
model_pain = fitlm(age_pa, avg_pain);
residuals_pain = model_pain.Residuals.Raw;

% For visualization here only, regress out age from the aperiodic exponents
model_apexp = fitlm(age_pa, apexp_pa);
residuals_apexp = model_apexp.Residuals.Raw;

% figure, plot(model_pain);
% figure, plot(model_apexp);

% Scatter plot of residuals
model_res = fitlm(model_pain.Residuals.Raw, model_apexp.Residuals.Raw);
f2 = figure('Units','centimeters','Position',[25 25 11 9]);
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
ylabel('Aperiodic exponent residuals');
ylim([-0.6 0.6]);
% title('H2');
box off;
legend off;
grid on;
% Save figure
saveas(f2,fullfile(figures_path,['hypothesis2_' analysis '.svg']));

% Save aperiodic exponents in a csv file
patients_sorted.exp_PFC = apexp_pa_original';
writetable(patients_sorted, fullfile(stats_path,['patients_PFC_' analysis '.csv']));



