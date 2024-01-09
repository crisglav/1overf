%% Plotting whole brain (top view), highlighting 100 ROIs and color-coding according to the exponent or the offset

function [main_figure, difference_figure] = wholebrain_plot (params,patients, healthy, difference)

% Load surface structure
surf = ft_read_headshape('surface_white_both.mat');

pos = params.sourcemodel_atlas.pos;

cmax_p = max(patients);
cmin_p = min(patients);
cmax_hc = max(healthy);
cmin_hc = min(healthy);

% Creating a figure with tiled layout
main_figure = figure('Units','centimeters','Position',[0 0 20 10]);

index = fix((patients-cmin_p)/(cmax_p-cmin_p)*256)+1;
rgb = squeeze(ind2rgb(index,parula(256)));
tcl = tiledlayout(1,2);
ax_p = nexttile;
ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv','facealpha',0.2);
ft_plot_mesh(pos, 'vertexsize',20, 'vertexcolor',rgb);
title('Patients')

index = fix((healthy-cmin_hc)/(cmax_hc-cmin_hc)*256)+1;
rgb = squeeze(ind2rgb(index,parula(256)));

ax_hc = nexttile;
ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv','facealpha',0.2);
ft_plot_mesh(pos, 'vertexsize',20, 'vertexcolor',rgb);
title('Healthy')

% Set colormap and color limits for all subplots
set(ax_p, 'Colormap', parula, 'CLim', [min(cmin_p,cmin_hc) max(cmax_p,cmax_hc)])
set(ax_hc, 'Colormap', parula, 'CLim', [min(cmin_p,cmin_hc) max(cmax_p,cmax_hc)])
% assign color bar to one tile
colorbar(ax_p(end),'eastoutside');
colorbar(ax_hc(end),'eastoutside');

%% Plotting the difference between groups 
difference_figure = figure('Units','centimeters','Position',[0 0 10 10]);

cmax_diff = max(difference);
cmin_diff = min(difference);
lim = max(abs(cmin_diff), abs(cmax_diff));


index = fix((difference - cmin_diff)/(cmax_diff-cmin_diff)*256)+1;
rgb = squeeze(ind2rgb(index,parula(256)));

ax_diff = nexttile;
ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv','facealpha',0.2);
ft_plot_mesh(pos, 'vertexsize',20, 'vertexcolor',rgb);

% Set colormap and color limits
set(ax_diff, 'Colormap', parula , 'CLim', [cmin_diff cmax_diff]);
% assign color bar to one tile
colorbar(ax_diff(end),'eastoutside');
