% Prepare the data for statistical analyisis and generate plots for the main analysis
%
% Cristina Gil Avila, TUM, 12.02.2024

clear all, close all

% Add statisical functions and plotting functions
% run('../toolboxes/bayes_factor/installBayesFactor.m')
addpath('../toolboxes/raincloudplots');
addpath('fooof_matlab');

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
rois = {'PFC'};
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

        er(iSubj) = fm.error(1); % mae
        r2(iSubj) = fm.r_squared;
        ex(iSubj) = fm.aperiodic_params(end);
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

%% Analyses

%% Hypothesis 1: Do aperiodic exponents in the mPFC differ between patients and healthy participants?
% Plot aperiodic exponents
cb = lines(2);
f1 = figure();
ax = gca;
h1 = raincloud_plot_vertical(exponents.PFC(hc_mask),'box_dodge',1, 'color', cb(1,:), 'alpha', 0.5, 'box_dodge_amount', 0.5,'dot_dodge_amount',0.5);
h2 = raincloud_plot_vertical(exponents.PFC(pa_mask), 'box_dodge', 1, 'color', cb(2,:), 'alpha', 0.5, 'box_dodge_amount', 1,'dot_dodge_amount',1);
 
legend([h1{1} h2{1}], {'HC', 'PA'});
title('Aperiodic exponents')
set(ax,'XTick',[],'XTickLabel',[]);
ylabel('Exponents');
box off

% Check age covariates with ANCOVA in JASP
participants_sorted.exp_PFC = exponents.PFC';
writetable(participants_sorted, fullfile(stats_path,['exp_PFC_' analysis '.csv']));

% Save figure
saveas(f1,fullfile(figures_path,'hypothesis1.svg'));

%% Hypothesis 2: Are changes in aperiodic exponents spatially specific for the mPFC?

% Ploting
% rearrange data
for i = 1:length(rois)
    roi = rois{i};
    data{i,1} = exponents.(roi)(hc_mask);
    data{i,2} = exponents.(roi)(pa_mask);
end

colors = repmat(lines(2), [1 1 3]);
f2 = figure();
ax = gca;
h   = rm_raincloud(data, cb,1, 'ks', [],1);
% h = rm_raincloud_cg(data,'colours',colors,'bandwidth',[],'dist_plots',0.1,...
%     'line_width',1,'raindrop_size',10,'opacity',0.4);
box off
% ax.Color = 'none';
set(ax,'YTickLabel',flip(rois));
xlabel('Exponents');
title('Aperiodic exponents - ROIs')

% Do statistics in JASP if there is evidence for H1
% participants_sorted.exp_PFC = exponents.PFC';
% participants_sorted.exp_S1 = exponents.S1';
% participants_sorted.exp_VIS = exponents.VIS';
% writetable(participants_sorted, fullfile(results_path,['exp_allROIs_' analysis '.csv']));

% Save figure
% exportgraphics(f2,fullfile(figures_path,'hypothesis2.jpg'));

%% Hypothesis 3: Do aperiodic exponents correlate with pain in patients?

% Sort the patients table
id = cellfun(@(x) str2double(x(5:7)),patients.participant_id,'UniformOutput',false);
id = cell2mat(id);
[id_sorted,ix] = sort(id);
patients_sorted = patients(ix,:);

% Load PFC fooof files
pfc_files = dir(fullfile(fooof_path,'PFC','*.mat'));
% Extracf bids Id from files names
fooofid = cellfun(@(x) str2double(x(5:7)),{pfc_files.name},'UniformOutput',false);
fooofid = cell2mat(fooofid);

% Real patients mask
mask_pa = ismember(fooofid,id);
% Check that you did it well
% all(ismember(participants_sorted(mask,:).participant_id,patients_sorted.participant_id));

% Extract aperiodic exponents for the patients
exp_pa = exponents.PFC(mask_pa);
patients_sorted.exp_PFC = exp_pa';

% Discard the patient from which we don't have a pain rating
mask = ~isnan(patients_sorted.avg_pain);
avg_pain = patients_sorted.avg_pain(mask);
exp_pa = exp_pa(mask)';
age = patients_sorted.age(mask);

% For correlating avg_pain and exp_pa we want to regress out age to account
% for this covariate.
% Fit a linear model to the data and obtain the residulas
model_pain = fitlm(age,avg_pain);
model_apexp = fitlm(age,exp_pa);
figure, plot(model_pain);
figure, plot(model_apexp);
residuals_pain = model_pain.Residuals.Raw;
residuals_apexp = model_apexp.Residuals.Raw;

% Scatter plot of residuals
model_res = fitlm(model_pain.Residuals.Raw, model_apexp.Residuals.Raw);
f3_0 = figure;
h = plot(model_res);
h(1).Marker = 'o';
h(1).MarkerFaceColor = cb(2,:);
h(1).MarkerEdgeColor = 'none';
h(2).Color = [0.7 0.7 0.7];
h(2).LineWidth = 1.5;
h(3).Color = [0.7 0.7 0.7];
h(3).LineWidth = 1.5;
h(4).Color = [0.7 0.7 0.7];
h(4).LineWidth = 1.5;
xlabel('Pain residuals');
ylabel('Aperiodic exponent residuals');
title('Correlation aperidoic exponents - pain ratings');
box off;

% save residuals to the table
patients_sorted.residuals_pain(mask) = model_pain.Residuals.Raw;
patients_sorted.residuals_apexp(mask) = model_apexp.Residuals.Raw;

% Scatter plot
f3 = figure();
model = fitlm(avg_pain,exp_pa);
h = plot(model);
h(1).Marker = 'o';
h(1).MarkerFaceColor = cb(2,:);
h(1).MarkerEdgeColor = 'none';
h(2).Color = [0.7 0.7 0.7];
h(2).LineWidth = 1.5;
h(3).Color = [0.7 0.7 0.7];
h(3).LineWidth = 1.5;
h(4).Color = [0.7 0.7 0.7];
h(4).LineWidth = 1.5;
ylabel('Aperiodic exponent');
xlabel('Pain rating');
xlim([1 10])
title('Correlation aperidoic exponents - pain ratings');
box off;

% Save to JASP
writetable(patients_sorted, fullfile(stats_path,['patients_PFC_' analysis '.csv']));

% Save figure
saveas(f3,fullfile(figures_path,'hypothesis3.svg'));
saveas(f3_0,fullfile(figures_path,'hypothesis3_residuals.svg'));

