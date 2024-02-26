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
params.Ncores = 25;
if(isempty(gcp('nocreate')))
    parObj = parpool(params.Ncores);
end
%% Load data and fit fooof

% Participants file
participants = readtable(fullfile(params.RawDataPath,'participants_clean.tsv'),'Filetype','text');
participants.group = categorical(participants.group);

% Order the participants table in descending order by bidsID
participants_id = cellfun(@(x) str2double(x(5:7)),participants.participant_id,'UniformOutput',false);
participants_id = cell2mat(participants_id);
[participants_id_sorted,ix] = sort(participants_id);
participants_sorted = participants(ix,:);

% Separate patients and healthy participants groups
pa_mask = participants_sorted.group == 'pa';
hc_mask = participants_sorted.group == 'hc';

% Power files
power_files = dir(fullfile(power_path, '*.mat'));
nSubj = length(power_files);
nRoi = 100;

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
save(fullfile(out_path, 'e3_whole_brain.mat'),'apexp','offset','pa_mask','hc_mask');


% Separate data into groups
apexp_pa = apexp(pa_mask,:);
apexp_hc = apexp(hc_mask,:);
offset_pa = apexp(pa_mask,:);
offset_hc = apexp(hc_mask,:);
    

%% Plots
surf = ft_read_headshape('surface_white_both.mat');
pos = params.sourcemodel_atlas.pos;

% Calculating mean exponent for each ROI, within each group:
apexp_pa_avg = mean(apexp_pa);
apexp_hc_avg = mean(apexp_hc);

% Calculating mean offset for each ROI, within each group:
offset_pa_avg = mean(offset_pa);
offset_hc_avg = mean(offset_hc);

% Figure 1: Exponents mean values Healthy controls vs. Patients. and
% difference between groups
f_apexp = wholebrain_plot(apexp_pa_avg,apexp_hc_avg,pos,surf);
suptitle('Aperiodic Exponent - Whole Brain Analysis')
saveas(f_apexp,fullfile(figures_path,'wholebrain_exponent_PA_HC.fig'));


% Figure 2: Offsets mean values Healthy controls vs. Patients. and
% difference between groups
f_offset = wholebrain_plot(offset_pa_avg, offset_hc_avg,pos,surf);
suptitle('Aperiodic Offset - Whole Brain Analysis')
saveas(f_offset,fullfile(figures_path,'wholebrain_offset_PA_HC.fig'));


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
tcl = tiledlayout(1,3);
try
    colors = plasma;
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