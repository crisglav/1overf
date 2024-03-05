% Estimate aperiodic exponents and offsets at 100 different brain ROIs
%
% Cristina Gil, Flaminia Palloti, TUM, 26.02.2024

%% Settings

clear all,
close all;

% Add fieldtrip and analysis functions
addpath('/rechenmagd4/toolboxes_and_functions/fieldtrip');
ft_defaults;
addpath('../fooof_matlab');
addpath('../../toolboxes/matplotlib');
addpath('../../toolboxes');

% Load parameters
load('../../results/features/params.mat','params');

% Create output folders
out_path = '../../results/features/fooof_matlab/whole_brain';
if ~exist(out_path,'dir')
    mkdir(out_path)
end
figures_path = '../../results/figures';
power_path = '../../results/features/power';

% For plotting
surf = ft_read_headshape('surface_white_both.mat');
pos = params.sourcemodel_atlas.pos;

% Define the number of cores for parallelization
% params.Ncores = 25;
% if(isempty(gcp('nocreate')))
%     parObj = parpool(params.Ncores);
% end

%% Demographics data
% Participants file
participants = readtable(fullfile(params.RawDataPath,'participants_clean.tsv'),'Filetype','text');
participants.group = categorical(participants.group);

% Order the participants table in descending order by bidsID
participants_id = cellfun(@(x) str2double(x(5:7)),participants.participant_id,'UniformOutput',false);
participants_id = cell2mat(participants_id);
[participants_id_sorted,ix] = sort(participants_id);
participants_sorted = participants(ix,:);
age = participants_sorted.age;

% Separate patients and healthy participants groups
pa_mask = participants_sorted.group == 'pa';
hc_mask = participants_sorted.group == 'hc';

% Load patients data
patients = readtable(fullfile(params.RawDataPath,'patients_clean.tsv'),'Filetype','Text');

% Order the patients table in descending order by bidsID
patients_id = cellfun(@(x) str2double(x(5:7)),patients.participant_id,'UniformOutput',false);
patients_id = cell2mat(patients_id);
[~,ix] = sort(patients_id);
patients_sorted = patients(ix,:);

nRoi = 100;
nSubj = height(participants);
nPa = height(patients);
%% Fooof data
% Check if the data exists and load it
if exist(fullfile(out_path, 'e3_whole_brain.mat'),'file')
    load(fullfile(out_path, 'e3_whole_brain.mat'));
else

    % Power files
    power_files = dir(fullfile(power_path, '*.mat'));
    
    % Preallocate variables for aperiodic exponents and offsets
    apexp = nan(nSubj,nRoi);
    offset = nan(nSubj,nRoi);
    
    parfor iSubj=1:nSubj
        
        fileID = power_files(iSubj).name;
        
        % Loading power data
        power_struct = load(fullfile(power_path,fileID));
        freq = power_struct.power.freq;
        pow = power_struct.power.powspctrm;
        
        % Fit a fooofGroup object in the 2 - 40 Hz freq range, which fits a
        % power spectra of several channels
        fm = fooofGroup('freq_range',[2,40],'peak_width_limits',[1 12],'aperiodic_mode','fixed');
        fm = add_data(fm,freq,pow);
        fm = fit(fm);
        
        % Extract aperiodic exponents and offsets for the 100 ROIS
        apexp(iSubj,:) = cellfun(@(x) x.aperiodic_params(end), fm.group_results);
        offset(iSubj,:) = cellfun(@(x) x.aperiodic_params(1), fm.group_results);
    end
    % Save into disk
    save(fullfile(out_path, 'e3_whole_brain.mat'),'apexp','offset','pa_mask','hc_mask','age');
end


%% Hypothesis 1: Do aperiodic exponents in 100 rois differ between patients and healthy participants?

% Regress out age from the aperiodic exponents
apexp_res = nan(size(apexp));
for iRoi=1:nRoi
    model_apexp = fitlm(age,apexp(:,iRoi));
    apexp_res(:,iRoi) = model_apexp.Residuals.Raw;
end

% Separate residuals per groups
apexp_res_pa = apexp_res(pa_mask,:);
apexp_res_hc = apexp_res(hc_mask,:);

% Statistics
% As, as far as we know there is no possiblity to correct for multiple comparisions
% with bayesian statistics. Therfore we performed 100 frequentist t-tests bewtween groups
% (one test per ROI) and corrected the p-values with the false discovery
% rate.
pval_h1 = nan(nRoi,1);
for iRoi = 1:nRoi
    [~,pval_h1(iRoi)] = ttest2(apexp_res_hc(:,iRoi),apexp_res_pa(:,iRoi),'tail','both');
end
[~, ~, padj_h1] = fdr(pval_h1);
any(padj_h1 < 0.05);

%% Plot
f1 = figure('Units','centimeters','Position',[0 0 30 10]);
tiledlayout(1,3);

% Create a color scale for data from the patients and healthy, where the
% color at the bottom of the scale is min(min(data_pa),min(data_hc)) and
% the color at the top of the scale is max(max(data_pa),max(data_hc))
data_hc = mean(apexp_res_hc);
data_pa = mean(apexp_res_pa);
cmin = min(min(data_pa),min(data_hc));
cmax = max(max(data_pa),max(data_hc));

ax = nexttile;
ax = wholebrain_plot(ax,data_hc,cmin,cmax,surf,pos);
title('Healthy');

ax = nexttile;
ax = wholebrain_plot(ax,data_pa,cmin,cmax,surf,pos);
title('Patients');

data_diff = mean(apexp_res_hc)-mean(apexp_res_pa);
% cmin = min(data_diff);
% cmax = max(data_diff);
cmin = -max(abs(data_diff));
cmax = max(abs(data_diff));

ax = nexttile;
ax = wholebrain_plot(ax,data_diff,cmin,cmax,surf,pos);
title('Healthy - Patients');

suptitle('Age-corrected aperiodic exponent')
saveas(f1,fullfile(figures_path,'e3_whole_brain_h1_plasma.png'));

%% Hypothesis 2: Do aperiodic exponents in 100 rois correlate with pain intensity in patients?
% Get avg pain ratings and age and discard the patient from which we don't have a pain rating
pain_mask = ~isnan(patients_sorted.avg_pain);
avg_pain = patients_sorted.avg_pain(pain_mask);
age_pa = patients_sorted.age(pain_mask);
apexp_pa = apexp(pa_mask,:);
apexp_pa = apexp_pa(pain_mask,:);

% Regress out age from the pain ratings
model_pain = fitlm(age_pa, avg_pain);
pain_res = model_pain.Residuals.Raw;

% Regress out age from aperiodic exponents
apexp_res_pa = nan(size(apexp_pa));
for iRoi=1:nRoi
    model_apexp_pa = fitlm(age_pa,apexp_pa(:,iRoi));
    apexp_res_pa(:,iRoi) = model_apexp_pa.Residuals.Raw;
end

% Statistics
% Correlate aperiodic exponents residuals with pain residuals with
% Pearson's correlation for each ROI and correct all the p-values with fdr.
[rho, pval_h2] = corr(pain_res,apexp_res_pa);
[~, ~, padj_h2] = fdr(pval_h2);
any(padj_h2 < 0.05);

%% Plot
f2 = figure('Units','centimeters','Position',[0 0 10 10]);

% Color-code the correlation
cmin = -max(abs(rho));
cmax = max(abs(rho));

ax = gca;
ax = wholebrain_plot(ax,rho,cmin,cmax,surf,pos);
title('Correlation');
saveas(f2,fullfile(figures_path,'e3_whole_brain_h2_plasma.png'));
%%
function ax = wholebrain_plot (ax,data,cmin,cmax,surf,pos)
try
    colors = plasma;
catch
    colors = parula(256);
end
% Assign data to colors
index = fix((data-cmin)/(cmax-cmin)*256)+1;
rgb = squeeze(ind2rgb(index,colors));

% Plot walnut brain
ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv','facealpha',0.2,'facecolor','brain');

% Overlay ROIS with the data color-coded
ft_plot_mesh(pos, 'vertexsize',20, 'vertexcolor',rgb);

% Set colormap and color limits for all subplots
set(ax, 'Colormap', colors,'CLim', [cmin cmax])

% Add colorbar
colorbar(ax(end),'eastoutside');
end