% Estimate aperiodic exponents and offsets at 100 different brain ROIs
%
% Cristina Gil, Flaminia Palloti, TUM, 26.02.2024

%% Settings

clear all,
close all;

% Add fieldtrip and analysis functions
addpath('/rechenmagd4/toolboxes_and_functions/fieldtrip');
ft_defaults;
addpath('../analysis_functions');
addpath('../fooof_matlab');
addpath('../../toolboxes/raincloudplots');
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

% Define the number of cores for parallelization
% params.Ncores = 25;
% if(isempty(gcp('nocreate')))
%     parObj = parpool(params.Ncores);
% end
%% Load data and fit fooof

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

nRoi = 100;
nSubj = height(participants);

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

% Regress out age from aperiodic exponents and offsets (for visualization)
apexp_res = nan(size(apexp));
offset_res = nan(size(offset));
for iRoi=1:nRoi
    model_apexp = fitlm(age,apexp(:,iRoi));
    apexp_res(:,iRoi) = model_apexp.Residuals.Raw;
    model_offset = fitlm(age,offset(:,iRoi));
    offset_res(:,iRoi) = model_offset.Residuals.Raw;
end

% Separate data into groups
apexp_pa = apexp(pa_mask,:);
apexp_hc = apexp(hc_mask,:);
offset_pa = offset(pa_mask,:);
offset_hc = offset(hc_mask,:);

apexp_res_pa = apexp_res(pa_mask,:);
apexp_res_hc = apexp_res(hc_mask,:);
offset_res_pa = offset_res(pa_mask,:);
offset_res_hc = offset_res(hc_mask,:);

%% Statistics
% As, as far as we know there is no possiblity to correct for multiple comparisions
% with bayesian statistics, we performed 100 frequentist t-tests bewtween groups
% (one test per ROI) and corrected the p-values with the false discovery
% rate.
p_value_apexp = nan(nRoi,1);
p_value_apexp_res = nan(nRoi,1);
for iRoi = 1:nRoi
    [~,p_value_apexp(iRoi)] = ttest2(apexp_hc(:,iRoi),apexp_pa(:,iRoi),'tail','both');
    [~,p_value_apexp_res(iRoi)] = ttest2(apexp_res_hc(:,iRoi),apexp_res_pa(:,iRoi),'tail','both');
end
[pthr, pcor, padj] = fdr(p_value_apexp_res);
[pthr, pcor, padj] = fdr(p_value_apexp);



%% Plots
surf = ft_read_headshape('surface_white_both.mat');
pos = params.sourcemodel_atlas.pos;

% % Figure 1: Average aperiodic exponents, difference between groups
% f_apexp = wholebrain_plot(mean(apexp_pa),mean(apexp_hc),pos,surf);
% suptitle('Aperiodic exponent')
% saveas(f_apexp,fullfile(figures_path,'e3_whole_brain_apexp.fig'));
% saveas(f_apexp,fullfile(figures_path,'e3_whole_brain_apexp.svg'));


% % Figure 2: average aperiodic offsets, difference between groups
% f_offset = wholebrain_plot(mean(offset_pa), mean(offset_hc),pos,surf);
% suptitle('Aperiodic offset')
% saveas(f_offset,fullfile(figures_path,'e3_whole_brain_offset.fig'));
% saveas(f_offset,fullfile(figures_path,'e3_whole_brain_offset.svg'));

%%
% Figure 3: Age-corrected aperiodic exponents
f_apexp_res = wholebrain_plot(mean(apexp_res_pa),mean(apexp_res_hc),pos,surf);
suptitle('Age-corrected aperiodic exponent')
% saveas(f_apexp,fullfile(figures_path,'e3_whole_brain_apexp_res.fig'));
saveas(f_apexp_res,fullfile(figures_path,'e3_whole_brain_apexp_res.svg'));


% % Figure 4: Age-corrected aperiodic offsets
% f_offset_res = wholebrain_plot(mean(offset_res_pa),mean(offset_res_hc),pos,surf);
% suptitle('Age-corrected aperiodic offset')
% % saveas(f_offset_res,fullfile(figures_path,'e3_whole_brain_offset_res.fig'));
% saveas(f_offset_res,fullfile(figures_path,'e3_whole_brain_offset_res.svg'));
% 
% %% Plot rois
% index = fix((1:100-1)/(99)*256)+1;
% rgb = squeeze(ind2rgb(1:100,plasma));
% figure;
% ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv','facealpha',0.2);
% ft_plot_mesh(pos(53,:), 'vertexsize',20);
%%
function [main_figure] = wholebrain_plot (patients,healthy,pos,surf)

% Difference between groups
difference = healthy - patients;

% Color limits
cmax_p = max(patients);
cmin_p = min(patients);
cmax_hc = max(healthy);
cmin_hc = min(healthy);
cmax_diff = max(difference);
cmin_diff = min(difference);

% Create a figure with tiled layout
main_figure = figure('Units','centimeters','Position',[0 0 30 10]);
tiledlayout(1,3);
try
    colors = viridis;
catch
    colors = parula(256);
end
% Healthy
index = fix((healthy-cmin_hc)/(cmax_hc-cmin_hc)*256)+1;
rgb = squeeze(ind2rgb(index,colors));
ax_hc = nexttile;
ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv','facealpha',0.2);
ft_plot_mesh(pos, 'vertexsize',20, 'vertexcolor',rgb);
title('Healthy')

% Patients
index = fix((patients-cmin_p)/(cmax_p-cmin_p)*256)+1;
rgb = squeeze(ind2rgb(index,colors));
ax_p = nexttile;
ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv','facealpha',0.2);
ft_plot_mesh(pos, 'vertexsize',20, 'vertexcolor',rgb);
title('Patients')

% Difference
index = fix((difference - cmin_diff)/(cmax_diff-cmin_diff)*256)+1;
rgb = squeeze(ind2rgb(index,colors));
ax_diff = nexttile;
ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv','facealpha',0.2);
ft_plot_mesh(pos, 'vertexsize',20, 'vertexcolor',rgb);
title('Healthy - Patients')

%%
% Set colormap and color limits for all subplots
set(ax_hc, 'Colormap', colors, 'CLim', [min(cmin_p,cmin_hc) max(cmax_p,cmax_hc)])
set(ax_p, 'Colormap', colors, 'CLim', [min(cmin_p,cmin_hc) max(cmax_p,cmax_hc)])
set(ax_diff, 'Colormap', colors, 'CLim', [cmin_diff cmax_diff]);

% assign color bar to one tile
colorbar(ax_hc(end),'eastoutside');
colorbar(ax_p(end),'eastoutside');
colorbar(ax_diff(end),'eastoutside');
end