% Analyze fooof fitting from matlab implementation
clear all, close all

% Add statisical functions and plotting functions
run('bayes_factor/installBayesFactor.m')
addpath('raincloudplots');
% Paths
results_path = '../results/jasp';
fooof_path = '../results/fooof_matlab/';
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

        er(iSubj) = fm.error_mae;
        r2(iSubj) = fm.r2;
        ex(iSubj) = fm.aperiodic_params(2);

        clear error r_squared aperiodic_params gaussian_params peak_params 
    end
    
    errors.(roi) = er;
    r_squareds.(roi) = r2;
    exponents.(roi) = ex;
    
end

% Plot model errors and r_squared
figure;
tlc = tiledlayout(1,2);
nexttile

scatter(1:n,er,'filled');
xlabel('Recording');
ylabel('error');
title('Model error');

nexttile
scatter(1:n,r2,'filled');
xlabel('Recording');
ylabel('R^2');
title('R^2');

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
% writetable(participants_sorted, fullfile(results_path,['exp_PFC_' analysis '.csv']));

% Save figure
% exportgraphics(f1,fullfile(results_path,'hypothesis1.jpg'));

%% Hypothesis 2: Changes in mPFC exponent are spatially specific for the PFC

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

% Do statistics in JASP
participants_sorted.exp_PFC = exponents.PFC';
participants_sorted.exp_S1 = exponents.S1';
participants_sorted.exp_VIS = exponents.VIS';
% writetable(participants_sorted, fullfile(results_path,['exp_allROIs_' analysis '.csv']));

% Save figure
% exportgraphics(f2,fullfile(results_path,'hypothesis2.jpg'));

%% Hypothesis 3: Exponents correlate with pain in patients

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
all(ismember(participants_sorted(mask,:).participant_id,patients_sorted.participant_id));

% Extract fooof values for the patients
exp_pa = exponents.PFC(mask);
patients_sorted.exp_PFC = exp_pa';

% Plotting
f3 = figure();
ax = gca;
scatter(exp_pa,patients_sorted.avg_pain,[],cb(2,:),'filled');
xlabel('Exponent');
ylabel('Pain rating');
title('Correlation aperidoic exponents - pain ratings');

% Correlate exponents with the pain ratings
avg_pain = patients_sorted.avg_pain(~isnan(patients_sorted.avg_pain));
exp_pa = exp_pa(~isnan(patients_sorted.avg_pain))';
[bf10,r,p] = bf.corr(exp_pa,avg_pain);

str = {sprintf('BF = %0.3f',bf10),sprintf('r = %0.3f',r),sprintf('p = %0.3f',p)};
annotation('textbox',[0.15 0.6 0.3 0.3],'String',str,'FitBoxToText','on','BackgroundColor','w');

% Save to JASP
% writetable(patients_sorted, fullfile(results_path,['patients_PFC_' analysis '.csv']));

% Save figure
% exportgraphics(f3,fullfile(results_path,'hypothesis3.jpg'));

