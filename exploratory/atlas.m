clear all; close all;
fieldtrip_path = '/rechenmagd4/toolboxes_and_functions/fieldtrip';
addpath(fieldtrip_path);
ft_defaults

% afni = ft_read_atlas(fullfile(fieldtrip_path,'template/atlas/afni/TTatlas+tlrc.HEAD'));


%% Source model: centroid positons from Schaefer atlas
clear h;
% Schaefer atlas
nnetworks = '17'; % 7 or 17
nsources = '100'; % 100, 200, ..., 1000
schaefer_path = ['/rechenmagd3/Experiments/2023_1overf/toolboxes/schaefer_parcellations/Schaefer2018_' nsources 'Parcels_' nnetworks 'Networks_order_FSLMNI152_1mm.Centroid_RAS.csv'];
schaefer = readtable(schaefer_path);

% Volume conduction model
volpath = fullfile(fieldtrip_path,'template','headmodel','standard_bem.mat');
load(volpath,'vol');

% Source model
cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = [schaefer.R, schaefer.A, schaefer.S];
cfg.unit = 'mm';
cfg.headmodel = volpath;
sourcemodel_atlas = ft_prepare_sourcemodel(cfg);
sourcemodel_atlas.coordsys = 'mni';


% ROI definition
rois = {'PFC','S1','VIS'};

% PFC definition
PFC_mask = contains(schaefer.ROIName,{'Cinga','PFCm','PFCmp'});
sch.PFC = schaefer(PFC_mask,:);
aux = cellfun(@(x) strsplit(x,'_'),sch.PFC.ROIName,'UniformOutput',false);
labels.PFC = unique(cellfun(@(x) x{4},aux,'UniformOutput',false));
pos.PFC = sourcemodel_atlas.pos(PFC_mask,:); % Points by PFC label

% Somatosensory cortex definition
S1_mask = contains(schaefer.ROIName,{'PostC'});
sch.S1 = schaefer(S1_mask,:);
aux = cellfun(@(x) strsplit(x,'_'),sch.S1.ROIName,'UniformOutput',false);
labels.S1 = unique(cellfun(@(x) x{4},aux,'UniformOutput',false));
pos.S1 = sourcemodel_atlas.pos(S1_mask,:);

% Visual cortex definition
VIS_mask = contains(schaefer.ROIName,{'ExStr_','Striate','StriCal'});
sch.VIS = schaefer(VIS_mask,:);
aux = cellfun(@(x) strsplit(x,'_'),sch.VIS.ROIName,'UniformOutput',false);
labels.VIS = unique(cellfun(@(x) x{4},aux,'UniformOutput',false));
pos.VIS = sourcemodel_atlas.pos(VIS_mask,:);



f = figure('Position',[0 0 1000 4000]);
t = tiledlayout(3,2,'TileSpacing','none');
for r=1:3
    roi = rois{r};
    
    c = lines(length(labels.(roi)));
    
    tx = nexttile;
    ft_plot_mesh(vol.bnd(3),'facealpha',0.05,'facecolor',[0.1 0.1 0.1],'edgecolor',[0.5 0.5 0.5],'edgealpha',0.5); % brain
    for i=1:length(labels.(roi))
    %     hold on
        idx  = find(cellfun(@(x) contains(x,['_' labels.(roi){i} '_']), sch.(roi).ROIName));
        ft_plot_mesh(pos.(roi)(idx,:), 'vertexsize',15, 'vertexcolor',c(i,:));   
    end
    view(90, 0);

    nexttile
    ft_plot_mesh(vol.bnd(3),'facealpha',0.05,'facecolor',[0.1 0.1 0.1],'edgecolor',[0.5 0.5 0.5],'edgealpha',0.5); % brain
    for i=1:length(labels.(roi))
        hold on
        idx  = find(cellfun(@(x) contains(x,['_' labels.(roi){i} '_']), sch.(roi).ROIName));
        ft_plot_mesh(pos.(roi)(idx,:), 'vertexsize',15, 'vertexcolor',c(i,:));
        h(i) = scatter(NaN,NaN,10,c(i,:),'filled'); % dummy for legend
    end
    legend(h,labels.(roi),'color','none','Location','EastOutside');
    view(0, 90);
end

% title (t,['Schaefer atlas ' nnetworks ' networks, ' nsources ' sources, mPFC = ' num2str(height(PFC)) ])
% saveas(f,'Schaefer atlas ROIS.jpeg');

%%
figure;
ft_plot_mesh(vol.bnd(3),'facealpha',0.05,'facecolor',[0.1 0.1 0.1],'edgecolor',[0.5 0.5 0.5],'edgealpha',0.5); % brain
roi = 'PFC';
for i=1:length(labels.(roi))
    %     hold on
    idx  = find(cellfun(@(x) contains(x,['_' labels.(roi){i} '_']), sch.(roi).ROIName));
    ft_plot_mesh(pos.PFC(idx,:), 'vertexsize',15, 'vertexcolor','k');
end
view(0, 90);



%% Points by network
PFC_networks = unique(cellfun(@(x) x{3},aux,'UniformOutput',false));
pos_PFC = sourcemodel_atlas.pos(PFC_mask,:);
c = lines(length(PFC_networks));
f = figure('Position',[651 613 972 684]);
ft_plot_mesh(vol.bnd(3),'facealpha',0.1,'facecolor',[0.1 0.1 0.1],'edgecolor',[1 1 1],'edgealpha',0.5); % brain
for i=1:length(PFC_networks)
    hold on
    idx  = find(cellfun(@(x) contains(x,['_' PFC_networks{i} '_']), PFC.ROIName));
    ft_plot_mesh(pos_PFC(idx,:), 'vertexsize',15, 'vertexcolor',c(i,:));
    h(i) = scatter(NaN,NaN,10,c(i,:),'filled'); % dummy for legend
end
legend(h,PFC_networks,'color','none');
view(0,90);
title (['Schaefer atlas ' nnetworks ' networks, ' nsources ' sources, PFC = ' num2str(height(PFC)) ])
% saveas(f,fullfile(params.preprocessed_data_path,'EEG_features','Schaefer atlas.svg'));
% close(f);
