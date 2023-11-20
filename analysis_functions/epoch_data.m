function data = epoch_data(params,data)
dataset = data.cfg.dataset;

% Cut the data into epochs
cfg = [];
cfg.length = params.EpochLength;
cfg.overlap = params.EpochOverlap;
data = ft_redefinetrial(cfg,data);

% Remove the epochs that contain a discontinuity (introduced in the
% preprocessing when removing bad time segments)
events = ft_read_event(dataset);
sample = [events.sample];
duration = [events.duration];
boundary = ceil(sample(~isnan(duration))); % ceil because the resampling in EEGLab makes some latencies at half timepoints
% Double check with plot_bad_timesegements_singlestudy
%     duration = duration(~isnan(duration));
%     dur = [0 duration(1:end-1)];
%     boundary + cumsum(dur); % this is the same as the first row of badsegs (the last element is not present but we don't care about the last segment removed if it is at the end of the chunk)
badIntervalMatrix = [boundary', boundary']; % first column start artifact, second column end artifact

cfg = [];
cfg.artfctdef.reject = 'complete';
cfg.artfctdef.clean_raw_data.artifact = badIntervalMatrix;
data = ft_rejectartifact(cfg,data);

end