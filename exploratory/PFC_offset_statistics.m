% Statistical analyis of the aperiodic offsets in the PFC
%
% Cristina Gil, TUM, 2023
clear all, close all

% Add statisical functions and plotting functions
run('../../toolboxes/bayes_factor/installBayesFactor.m')
addpath('../../toolboxes/raincloudplots');
addpath('../fooof_matlab');

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

%% Load fooof data
roi_path = fullfile(fooof_path,'PFC');
roi_files = dir(fullfile(roi_path,'*.mat'));
n = length(roi_files);
offset = nan(1,n);

% Load fooof parameters
for iSubj=1:n
    load(fullfile(roi_files(iSubj).folder,roi_files(iSubj).name));
    offset(iSubj) = fm.aperiodic_params(1);
end


%% Hypothesis 1: differences at PFC in aperiodic offset between patients and healthy participants
% Plot offsets 
cb = lines(2);
f1 = figure();
ax = gca;
h1 = raincloud_plot(offset(hc_mask), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', 0.15,'box_col_match', 0, ...
     'cloud_edge_col', cb(1,:),'plot_top_to_bottom',1);
h2 = raincloud_plot(offset(pa_mask), 'box_on', 1, 'color', cb(2,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,...
     'cloud_edge_col', cb(2,:),'plot_top_to_bottom',1);
 
legend([h1{1} h2{1}], {'HC', 'PA'});
title('Aperiodic offsets in PFC')
set(ax,'XLim', [-0.5 2.5], 'YLim', [-1.5 2.5]);
set(ax,'YTickLabel',[]);
xlabel('Offset');
box off

% Check age and gender covariates ANCOVA in JASP
participants_sorted.offset_PFC = offset';
writetable(participants_sorted, fullfile(stats_path,['offset_PFC_' analysis '.csv']));

% Save figure
exportgraphics(f1,fullfile(figures_path,'hypothesis1_offset.jpg'));

%% Hypothesis 3: Offsets correlate with pain levels in patients
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
mask = ismember(fooofid,id);
% Check that you did it well
% all(ismember(participants_sorted(mask,:).participant_id,patients_sorted.participant_id));

% Extract fooof values for the patients
offset_pa = offset(mask);

% Plotting
f3 = figure();
ax = gca;
scatter(offset_pa,patients.avg_pain,[],cb(2,:),'filled');
xlim([-0.2, 2])
ylim([0,10])
xlabel('Offset');
ylabel('Pain rating');
title('Correlation aperiodic offset - pain ratings in PFC');

% Save to do the statistical analysis in JASP
patients_sorted.offset_PFC = offset_pa';
writetable(patients_sorted, fullfile(stats_path,['offset_patients_PFC_' analysis '.csv']));

% Save figure
exportgraphics(f3,fullfile(figures_path,'hypothesis3_offset.jpg'));