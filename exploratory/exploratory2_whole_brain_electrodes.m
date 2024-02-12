% Whole brain analyisis in the elctrode space
%
% Cristina Gil, TUM, 20.12.2023

% Load data
clear all,
close all;
%% Settings

% Add fieldtrip and analysis functions
addpath('/rechenmagd4/toolboxes_and_functions/fieldtrip');
ft_defaults;
addpath('../analysis_functions');
addpath('../fooof_matlab');
% Load preprocessing parameters
load('../../results/features/params.mat');

% Output folders
figures_path = '../../results/figures';
expPath = '../../results/features/electrode_space';
if ~exist(expPath)
    mkdir(expPath)
end
% params.Ncores = 15;
% if(isempty(gcp('nocreate')))
%     parObj = parpool(params.Ncores);
% end
%% Load the data
% Load the tsv in which the participant labels are randomized to be
% completely agnostic to the groups.
participants = readtable(fullfile(params.RawDataPath,'participants_clean.tsv'),'Filetype','text');
participant_id = participants.participant_id;
n = height(participants);
nElectrodes = 62;

try
    load(fullfile(expPath, 'exponents.mat'))
catch
    exponents = nan(n,nElectrodes);
    parfor iSubj=1:n
        
        bidsID = participant_id{iSubj};
        bidsID = [bidsID '_task-closed'];
        
        try
            power = load(fullfile(expPath,[bidsID '_power.mat']))
            power = power.power;
        catch
            try
                % Load EEG preprocessed data
                data = load_preprocessed_data(params,bidsID);
            catch
                continue
            end
            
            % Cut the data into 2s epochs
            data = epoch_data(params,data);
            
            % Estimate power in the electrode space
            cfg = [];
            cfg.foilim = [1 100];
            cfg.method = 'mtmfft';
            cfg.taper = 'dpss';
            cfg.tapsmofrq = 1;
            cfg.output = 'pow';
            cfg.keeptrials ='no';
            power = ft_freqanalysis(cfg, data);
            parsave(fullfile(expPath, [bidsID '_power.mat']),power);
        end
        
        % Create a fooof object handling multiple channels, add data to the
        % object and fit the model to each channel
        fm = fooofGroup('freq_range',[2,40],'peak_width_limits',[1 12],'aperiodic_mode','fixed');
        fm = add_data(fm,power.freq,power.powspctrm);
        fm = fit(fm);
        parsave(fullfile(expPath, [bidsID '_fm.mat']),fm);

        % Extract exponents
        exponents(iSubj,:) = cellfun(@(x) x.aperiodic_params(end), fm.group_results);
        
    end
    
    pa_mask = categorical(participants.group) == 'pa';
    hc_mask = categorical(participants.group) == 'hc';
    
    save(fullfile(expPath, 'exponents.mat'),'exponents','pa_mask','hc_mask');
end

%% Topoplots with spatial distribution of aperiodic exponents
pa_exponents = exponents(pa_mask,:);
hc_exponents = exponents(hc_mask,:);

% Use the power strucutre as a basis to plot the topographic distributions
bidsID = [participant_id{1} '_task-closed'];
load(fullfile(expPath,[bidsID '_power.mat']))

top_hc = power;
top_hc.exponent = mean(hc_exponents)';
top_pa = power;
top_pa.exponent = mean(pa_exponents)';
top_diff = power;
top_diff.exponent = mean(hc_exponents)' - mean (pa_exponents)';

cfg = [];
cfg.xlim = [1, 1]; % Fake dimension (frequency) just for the plot
cfg.zlim = [1 1.5];
cfg.parameter = 'exponent';
cfg.comment = 'no';

fhc = figure('Units','centimeters','Position',[0 0 10 10]);
ft_topoplotER(cfg,top_hc);
colorbar;
title('Healthy');
saveas(fhc,fullfile(figures_path,'wholebrain_exponent_electrode_HC.svg'));

fpa = figure('Units','centimeters','Position',[0 0 10 10]);
ft_topoplotER(cfg,top_pa);
colorbar;
title('Patients');
saveas(fpa,fullfile(figures_path,'wholebrain_exponent_electrode_PA.svg'));

fdiff = figure('Units','centimeters','Position',[0 0 10 10]);
cfg.zlim = [];
% cfg.zlim = [-0.06 0.06];
ft_topoplotER(cfg,top_diff);
colorbar;
title('Healthy - Patients');
saveas(fdiff,fullfile(figures_path,'wholebrain_exponent_electrode_HCPA.svg'));


%% Load exponent data from disk
data = load_preprocessed_data(params,'sub-035_task-closed');
label = data.label;
layout = ft_prepare_layout([],data);

% Load all fm files
fm_all = {};
for iSubj = 1:n
    bidsID = participant_id{iSubj};
    bidsID = [bidsID '_task-closed'];
    
    load(fullfile(expPath,[bidsID '_fm.mat']));
    % Put the data in a structure that ft_multiplotER can read
    aux = power;
    aux.freq = fm.freqs;
    aux.powspctrm = fm.power_spectra;
    aux.ap_fit = fm.ap_fit_group;
    aux.fooofed_spectrum = fm.fooofed_spectrum_group;
    aux.exponent = cellfun(@(x) x.aperiodic_params(end), fm.group_results)';
    fm_all{iSubj} = aux;
end
fm_pa = fm_all(pa_mask);
fm_hc = fm_all(hc_mask);


%% Cluster-based permutatuion statistics

% Define neighbours of electrodes (based on Son's code)
cfg = [];
cfg.method = 'distance';
cfg.layout = layout;
cfg.neighbourdist = 0.18; % this value is taken from Son's code (18 mm)
neighbours = ft_prepare_neighbours(cfg, data.elec);

% Cluster-based permutations statistics
cfg = [];
cfg.method = 'montecarlo';
cfg.parameter = 'exponent';

cfg.correctm = 'cluster';
cfg.minnbchan = 2;
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.alpha = 0.05;
cfg.tail = 0;
cfg.correcttail = 'prob';
cfg.clustertail = 0;
cfg.numrandomization = 500;
cfg.statistic = 'indepsamplesT';
cfg.neighbours = neighbours;
cfg.design = [ones(1,length(fm_hc)), ones(1,length(fm_pa)) + 1];
cfg.ivar = 1;
cfg.frequency = [2 2]; % fake freq
cfg.avgoverfreq = 'no';

stat = ft_freqstatistics(cfg, fm_hc{:}, fm_pa{:});


%% Multiplot
cfg = [];
cfg.parameter = {'ap_fit','fooofed_spectrum','powspctrm'};
fm_pa_avg = ft_freqgrandaverage(cfg,fm_pa{:});
fm_hc_avg = ft_freqgrandaverage(cfg,fm_hc{:});
% fm_pa_avg.freq = log10(fm_pa_avg.freq);
% fm_hc_avg.freq = log10(fm_hc_avg.freq);

% % Multiplot
% figure;
% cfg = [];
% cfg.parameter = 'ap_fit';
% cfg.layout = layout;
% ft_multiplotER(cfg,fm_pa_avg,fm_hc_avg);

% Multiplot
figure;
cfg = [];
cfg.ylim = [-2.67 1.08]; % As in the preprocessed data
cfg.parameter = 'powspctrm';
cfg.layout = layout;
ft_multiplotER(cfg,fm_pa_avg,fm_hc_avg);


% % Load all power files
% power_all = {};
% for iSubj = 1:n
%     bidsID = participant_id{iSubj};
%     bidsID = [bidsID '_task-closed'];
%     
%     load(fullfile(expPath,[bidsID '_power.mat']));
%     power.powspctrm = log10(power.powspctrm);
%     power_all{iSubj} = power;
% end
% power_pa = power_all(pa_mask);
% power_hc = power_all(hc_mask);
% 
% cfg = [];
% power_pa_avg = ft_freqgrandaverage(cfg,power_pa{:});
% power_hc_avg = ft_freqgrandaverage(cfg,power_hc{:});
% % % log the power
% % power_pa_avg.logpow = log10(power_pa_avg.powspctrm);
% % power_hc_avg.logpow = log10(power_hc_avg.powspctrm);
% 
% % Multiplot
% figure;
% cfg = [];
% cfg.parameter = 'powspctrm';
% cfg.layout = layout;
% ft_multiplotER(cfg,power_pa_avg,power_hc_avg);