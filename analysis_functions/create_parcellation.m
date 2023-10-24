function params = create_parcellation(params)
% Create parcellation based on atlas and precompute the source model
atlas_table = readtable(params.AtlasPath);
parcellation.pos = [atlas_table.R, atlas_table.A, atlas_table.S];
parcellation.unit = 'mm';
parcellation.coordsys = 'mni';
parcellation.ROI = atlas_table.ROILabel;
parcellation.ROIlabel = atlas_table.ROIName;
params.parcellation = parcellation;

params.PFC_mask = contains(parcellation.ROIlabel,{'Cinga','PFCm','PFCmp'});
params.S1_mask = contains(parcellation.ROIlabel,{'PostC'});
params.VIS_mask = contains(parcellation.ROIlabel,{'ExStr_','Striate','StriCal'});

%% Source model 
% Source model: centroid positons from Schaefer atlas
atlas400 = readtable(params.AtlasPath);
cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = [atlas400.R, atlas400.A, atlas400.S];
cfg.unit = 'mm';
cfg.headmodel = params.HeadModelPath;
sourcemodel_atlas = ft_prepare_sourcemodel(cfg);
sourcemodel_atlas.coordsys = 'mni';
params.sourcemodel_atlas = sourcemodel_atlas;

end

