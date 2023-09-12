% Estimate power spectrum at mPFC from preprocessed data
clear all;
close all;

% Import paths
fid = fopen('paths.json');
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
paths = jsondecode(str);
clear fid raw str

% Add functions from the discover-eeg to read preprocessed data
addpath(genpath(paths.discover_eeg));
addpath(paths.fieldtrip);
ft_defaults;

% Load preprocessing parameters
fparams = '../data/blinded/derivatives_v2023_08_18/params.json';
fid = fopen(fparams);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
params = jsondecode(str);
clear fid raw str

% Define the results path
params.RawDataPath = '../data/blinded';
params.PreprocessedDataPath = '../data/blinded/derivatives_v2023_08_18';
params.SourcePath = '../results/source';
params.AtlasPath = 'schaefer_parcellations/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_1mm.Centroid_RAS.csv';
params.FiguresPath = '../results/figures';
%% Parcellation based on Schaefer atlas
atlas_table = readtable(params.AtlasPath);
parcellation.pos = [atlas_table.R, atlas_table.A, atlas_table.S];
parcellation.unit = 'mm';
parcellation.coordsys = 'mni';
parcellation.ROI = atlas_table.ROILabel;
parcellation.ROIlabel = atlas_table.ROIName;

PFC_mask = contains(atlas_table.ROIName,{'Cinga','PFCm','PFCmp'});
parcellation.PFC = PFC_mask;

S1_mask = contains(atlas_table.ROIName,{'PostC'});
parcellation.S1 = S1_mask;

VIS_mask = contains(atlas_table.ROIName,{'ExStr_','Striate','StriCal'});
parcellation.VIS = VIS_mask;
%% Load all subject ids
% participants = readtable(fullfile(params.RawDataPath,'participants_rand.tsv'),'Filetype','text');
% 
% n = height(participants);
% task = 'closed';
% 
% % Define the frequency band of interest
% params.FreqBand.fullSpectrum = [2 40];
% 
% % for iSubj=1:n
% iSubj = 1;
% bidsID = participants.participant_id{iSubj};
% bidsID = [bidsID '_task-' task];

% Compute source reconstruction - CHANGE TO THIS FUNCTION IN THE END
% source = compute_spatial_filter(params,bidsID,'fullSpectrum');
%% % Compute_spatial_filter function
% Hard-coded
participants = {'sub-001_task-closed','sub-002_task-closed','sub-003_task-closed','sub-004_task-closed','sub-005_task-closed'};
params.FreqBand.fullSpectrum = [2 40];

for iSubj=1:length(participants)
% iSubj = 1;
bidsID = participants{iSubj};

% Load EEG data
data = load_preprocessed_data(params,bidsID);

% Source model: centroid positons from Schaefer atlas
atlas400 = readtable(params.AtlasPath);
cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = [atlas400.R, atlas400.A, atlas400.S];
cfg.unit = 'mm';
cfg.headmodel = params.HeadModelPath;
sourcemodel_atlas = ft_prepare_sourcemodel(cfg);
sourcemodel_atlas.coordsys = 'mni';

% Bandpass the data in the relevant frequency band
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = params.FreqBand.('fullSpectrum');
datafilt = ft_preprocessing(cfg, data);

% Compute the average covaraciance matrix from the sensor data
cfg = [];
cfg.covariance = 'yes';
cfg.keeptrials = 'no';
cfg.removemean = 'yes';
tlock = ft_timelockanalysis(cfg,datafilt);

%  Computation of the spatial filter
% Forward model (leadfield)
cfg = [];
cfg.sourcemodel = sourcemodel_atlas;
cfg.headmodel = params.HeadModelPath;
cfg.normalize = 'yes';
lf = ft_prepare_leadfield(cfg, datafilt);

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

% % Parcellate the source structure based on the regions of interest (not needed, done later)
% cfg = [];
% cfg.parcellation = 'ROI';
% sourcep = ft_sourceparcellate(cfg,source,parcellation);

% % Plot source power
% surf = ft_read_headshape('surface_white_both.mat');
% tmpcfg = [];
% tmpcfg.method = 'nearest';
% tmpcfg.parameter = 'pow';
% sourceInterp = ft_sourceinterpolate(tmpcfg, source, surf);
%
% cfg = [];
% cfg.funparameter = 'pow';
% cfg.method = 'surface';
% cfg.funcolormap = 'jet';
% cfg.atlas = atlas_ft;
% ft_sourceplot(cfg,sourceInterp,surf);

%% Extract virtual channel data
cfg = [];
cfg.parcellation = 'ROI';
vdata_trials = ft_virtualchannel(cfg,data,source,parcellation);

% % ======== Testing =========
% % source.avg.mom (instantaneous dipole moment each source position, aka the virtual time series)
% % Create a data structure with the virtual data so that it can be given as input to ft_freqanalysis
% vdata_filt_s = [];
% vdata_filt_s.time = source.time;
% vdata_filt_s.label = parcellation.ROIlabel;
% vdata_filt_s.avg = cell2mat(source.avg.mom);
% 
% % The following is exactly the same as getting source.avg.mom when applied to tlock, but it is in
% % a ft structure that can be inputed to ft_freqanalysis
% cfg = [];
% cfg.parcellation = 'ROI';
% vdata_filt = ft_virtualchannel(cfg,tlock,source,parcellation);
% 
% % Apply the spatial filter to the unfiltered data. Option A and B yield the same result
% % Option A) Average over channels the sensor data (unfiltered), then apply spatial filter to the average
% % NOTE: We should not average over repetitions here, but when computing the power spectrum!!!!
% cfg = [];
% cfg.avgoverrpt = 'yes';
% [data_avg] = ft_selectdata(cfg, data);
% cfg = [];
% cfg.parcellation = 'ROI';
% vdataA = ft_virtualchannel(cfg,data_avg,source,parcellation);
% % Option B) Apply spatial filter to the sensor data (unfiltered), then average the source localized data
% cfg = [];
% cfg.parcellation = 'ROI';
% vdata_trials = ft_virtualchannel(cfg,data,source,parcellation);
% cfg = [];
% cfg.avgoverrpt = 'yes';
% [vdataB] = ft_selectdata(cfg, vdata_trials);
% 
% 
% % Check that the virtual data looks ok
% figure;
% plot(vdataA.time{:},vdataA.trial{1}(1,:));
% hold on
% plot(vdataB.time{:},vdataB.trial{1}(1,:));
% hold on
% plot(vdata_filt.time,vdata_filt.avg(1,:));
% hold on
% plot(vdata_filt_s.time,vdata_filt_s.avg(1,:));
% legend({'option A', 'option B','virtual filtered data','mom'})
%% Estimate power spectra

% Extract the power spectra for all the virtualdata trials, and average over trials
cfg = [];
cfg.foilim = [1 100];
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.tapsmofrq = 1;
cfg.pad = 5;
cfg.padtype = 'zero';
cfg.output = 'pow';
cfg.keeptrials ='no';
power = ft_freqanalysis(cfg, vdata_trials);

% ========== Testing =======
% This is not okay
% This is doing the power spectrum of the average across all time series.
% cfg = [];
% cfg.foilim = [1 100];
% cfg.method = 'mtmfft';
% cfg.taper = 'dpss';
% cfg.tapsmofrq = 1;
% cfg.pad = 5;
% cfg.padtype = 'zero';
% cfg.output = 'pow';
% cfg.keeptrials ='no';
% power_old = ft_freqanalysis(cfg, vdataA);
% 
% figure;
% plot(power.freq,power.powspctrm(1,:));
% hold on
% plot(power_old.freq,power_old.powspctrm(1,:));
% legend({'avg power of each vdata epoch, avgd across epochs','power of the avg vdata, avgd across epochs'})

%% Extract power at the PFC
% Extract power at the PFC and average across channels in the ROI
cfg = [];
cfg.channel = find(parcellation.PFC);
power_PFC = ft_selectdata(cfg,power);

cfg = [];
cfg.avgoverchan = 'yes';
avgpower_PFC = ft_selectdata(cfg,power_PFC);

figure();
plot(power_PFC.freq,power_PFC.powspctrm);
hold on
plot(avgpower_PFC.freq,avgpower_PFC.powspctrm,'k','LineWidth',2);



% Check if this is the same
cfg = [];
cfg.foilim = [1 100];
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.tapsmofrq = 1;
cfg.pad = 5;
cfg.padtype = 'zero';
cfg.output = 'pow';
cfg.keeptrials ='yes';
power2 = ft_freqanalysis(cfg, vdata_trials);

cfg = [];
cfg.channel = find(parcellation.PFC);
cfg.avgoverchan = 'yes';
% cfg.avgoverrpt = 'yes';
avgpower_PFC2 = ft_selectdata(cfg,power2);

cfg = [];
% cfg.channel = find(parcellation.PFC);
% cfg.avgoverchan = 'yes';
cfg.avgoverrpt = 'yes';
avgpower_PFC3 = ft_selectdata(cfg,avgpower_PFC2);

%% Power spectra at the sensor level
cfg = [];
cfg.foilim = [1 100];
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.tapsmofrq = 1;
cfg.pad = 5;
cfg.padtype = 'zero';
cfg.output = 'pow';
cfg.keeptrials ='no';
power_sensor = ft_freqanalysis(cfg, data);

%% Plot power spectra at the sensor and source levels
f = figure('Units','centimeters','Position',[0 0 20 16]);
tlc = tiledlayout(2,2);
nexttile()
% Plot power
plot(power_sensor.freq,power_sensor.powspctrm);
hold on
% Plot avg power across channels
cfg = [];
cfg.avgoverchan = 'yes';
avgpower_sensor = ft_selectdata(cfg,power_sensor);
plot(avgpower_sensor.freq,avgpower_sensor.powspctrm,'k','LineWidth',2);
title('Sensor space');

nexttile();
% Plot power
plot(power.freq,power.powspctrm);
% Plot avg power across channels
hold on,
cfg = [];
cfg.avgoverchan = 'yes';
avgpower = ft_selectdata(cfg,power);
plot(avgpower.freq, avgpower.powspctrm,'k','LineWidth',2);
title('Source space')

nexttile()
% Plot power
plot(power_sensor.freq,log(power_sensor.powspctrm));
hold on
% Plot avg power across channels
plot(avgpower_sensor.freq,log(avgpower_sensor.powspctrm),'k','LineWidth',2);
title('Sensor space log scale');

nexttile();
% Plot power
plot(power.freq,log(power.powspctrm));
% Plot avg power across channels
hold on,
plot(avgpower.freq, log(avgpower.powspctrm),'k','LineWidth',2);
title('Source space log scale')

saveas(f,fullfile(params.FiguresPath, [bidsID '_power.svg']));
end


%% Testing
% % Load EEG data
% data = load_preprocessed_data(params,bidsID);
% % Bandpass the data in the relevant frequency band
% cfg = [];
% cfg.bpfilter = 'yes';
% cfg.bpfreq = params.FreqBand.('fullSpectrum');
% data = ft_preprocessing(cfg, data);
% % Compute the average covaraciance matrix from the sensor data
% cfg = [];
% cfg.covariance = 'yes';
% cfg.keeptrials = 'no';
% cfg.removemean = 'yes';
% tlock = ft_timelockanalysis(cfg,data);
%
% % Reconstruct the virtual time series (apply spatial filter to sensor level
% % data)
% % Here you reconstruct each trial.
% cfg  = [];
% cfg.pos = source.pos;
% virtChan_data = ft_virtualchannel(cfg,data,source);
%
% % Let's try to average first the data (e.g. as in timelock) and then apply
% % the filter
% virtChan_data_avg = ft_virtualchannel(cfg,tlock,source);
%
% % theoretically, source.avg.mom should be the same as virtChan_data
% figure;
% plot(source.time,source.avg.mom{1});
%
% cfg = [];
% cfg.avgoverrpt = 'yes';
% [data2] = ft_selectdata(cfg, virtChan_data);
% hold on
% plot(source.time,data2.trial{1,1}(1,:));
%
% hold on
% plot(source.time,virtChan_data_avg.avg(1,:));
% legend('source.avg.mom','virtChan_data_trial','virtChan_data_avg')
%
% % The difference between this two lines is that the first one is the
% % reconstruction of the average trials
% % the second one is the reconstruction for each trial and then averaged.

%% % %% Aperiodic estimation in matlab
% cfg = [];
% % cfg.channel = find(PFC_mask);
% % cfg.avgoverchannel = 'yes';
% cfg.foilim = [1 100];
% cfg.method = 'mtmfft';
% cfg.taper = 'dpss';
% cfg.tapsmofrq = 1;
% cfg.pad = 5;
% cfg.padtype = 'zero';
% cfg.output = 'fooof_aperiodic';
% cfg.keeptrials ='no';
% power_fooof = ft_freqanalysis(cfg, vdata_trials);
% 
% % Brainstorm function included in fieldtrip
% opts_bst  = getfield(process_fooof('GetDescription'), 'options');
% % Fetch user settings, this is a chunk of code copied over from
% % process_fooof, to bypass the whole database etc handling.
% opt                     = ft_getopt(cfg, 'fooof', []);
% opt.freq_range          = ft_getopt(opt, 'freq_range', [2 40]);
% opt.peak_width_limits   = ft_getopt(opt, 'peak_width_limits', opts_bst.peakwidth.Value{1});
% opt.max_peaks           = ft_getopt(opt, 'max_peaks',         opts_bst.maxpeaks.Value{1});
% opt.min_peak_height     = ft_getopt(opt, 'min_peak_height',   opts_bst.minpeakheight.Value{1}/10); % convert from dB to B
% opt.aperiodic_mode      = ft_getopt(opt, 'aperiodic_mode',    opts_bst.apermode.Value);
% opt.peak_threshold      = ft_getopt(opt, 'peak_threshold',    2);   % 2 std dev: parameter for interface simplification
% opt.return_spectrum     = ft_getopt(opt, 'return_spectrum',   1);   % SPM/FT: set to 1
% opt.border_threshold    = ft_getopt(opt, 'border_threshold',  1);   % 1 std dev: proximity to edge of spectrum, static in Python
% % Matlab-only options
% opt.power_line          = ft_getopt(opt, 'power_line',        '50'); % for some reason it should be a string, if you don't want a notch, use 'inf'. Brainstorm's default is '60'
% opt.peak_type           = ft_getopt(opt, 'peak_type',         opts_bst.peaktype.Value);
% opt.proximity_threshold = ft_getopt(opt, 'proximity_threshold', opts_bst.proxthresh.Value{1});
% opt.guess_weight        = ft_getopt(opt, 'guess_weight',      opts_bst.guessweight.Value);
% opt.thresh_after        = ft_getopt(opt, 'thresh_after',      true);   % Threshold after fitting always selected for Matlab (mirrors the Python FOOOF closest by removing peaks that do not satisfy a user's predetermined conditions)
% 
% % Output options
% opt.sort_type  = opts_bst.sorttype.Value;
% opt.sort_param = opts_bst.sortparam.Value;
% opt.sort_bands = opts_bst.sortbands.Value;
% 
% hasOptimTools = 0;
% if exist('fmincon', 'file')
%     hasOptimTools = 1;
%     disp('Using constrained optimization, Guess Weight ignored.')
% end
%     
% opt = [];
% opt.freq_range = [2 40];
% [fs, fg] = process_fooof('FOOOF_matlab', avgpower_PFC.powspctrm, avgpower_PFC.freq, opt, hasOptimTools);
% 
% 
% % Plot
% figure();
% plot(power_PFC.freq,power_PFC.powspctrm);
% hold on
% plot(avgpower_PFC.freq,avgpower_PFC.powspctrm,'k','LineWidth',2);
% 
% % To mimick the fooof plot in python
% figure;
% freqrange = find(freq==2):find(freq==40);
% plot(freq(freqrange),log(avgpower_PFC.powspctrm(freqrange)),'k');
% xlabel('Frequency');
% ylabel('logPower');


