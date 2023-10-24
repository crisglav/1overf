function source = compute_spatial_filter(params,data,freqBand)
% Similar as in discover-EEG. Adjusted to make the computation faster

%% Bandpass the data in the relevant frequency band
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = params.FreqBand.(freqBand);
data = ft_preprocessing(cfg, data);

%% Source_model is precomputed in create_parcellation.m

%% Compute the covariance matrix from the data
% First normalize time axis of the data (otherwise it cracks).
% Here we loose the temporal order of the epochs
temptime = data.time{1};
[data.time{:}] = deal(temptime);

% Compute the average covaraciance matrix from the sensor data
cfg = [];
cfg.covariance = 'yes';
cfg.keeptrials = 'no';
cfg.removemean = 'yes';
tlock = ft_timelockanalysis(cfg,data);

%%  Computation of the spatial filter
% Forward model (leadfield)
cfg = [];
cfg.sourcemodel = params.sourcemodel_atlas;
cfg.headmodel = params.HeadModelPath;
cfg.normalize = 'yes';
lf = ft_prepare_leadfield(cfg, data);

% Spatial filter
cfg = [];
cfg.method = 'lcmv';
cfg.keeptrials = 'yes';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.lambda = '5%';
cfg.lcmv.fixedori = 'yes';
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.weightnorm = 'arraygain';
% cfg.lcmv.weightnorm = 'nai';
cfg.sourcemodel = lf;
source = ft_sourceanalysis(cfg, tlock);

% save(fullfile(params.SourcePath,[bidsID '_source_' freqBand '.mat']),'source');
end