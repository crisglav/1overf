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
cb = colorbar(ax(end),'eastoutside');
end