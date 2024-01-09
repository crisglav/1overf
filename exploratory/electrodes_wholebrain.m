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

load(fullfile(expPath, 'exponents.mat'));

%% Load the data
% Load the tsv in which the participant labels are randomized to be
% completely agnostic to the groups.
participants = readtable(fullfile(params.RawDataPath,'participants_clean.tsv'),'Filetype','text');
participant_id = participants.participant_id;
n = height(participants);
nElectrodes = 62;

exponents = nan(n,nElectrodes);
for iSubj=1:n
    
    bidsID = participant_id{iSubj};
    bidsID = [bidsID '_task-closed'];
    
    try
        load(fullfile(expPath,[bidsID '_power.mat']))
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

    % Extract exponents
    exponents(iSubj,:) = cellfun(@(x) x.aperiodic_params(end), fm.group_results);
    
end

pa_mask = categorical(participants.group) == 'pa';
hc_mask = categorical(participants.group) == 'hc';

pa_exponents = exponents(pa_mask,:);
hc_exponents = exponents(hc_mask,:);
save(fullfile(expPath, 'exponents.mat'),'exponents','pa_mask','hc_mask');

%% Topoplots with spatial distribution of aperiodic exponents

% Use the power strucutre as a basis to plot the topographic distributions
top_pa = power;
top_pa.exponent = median(pa_exponents)';
top_hc = power;
top_hc.exponent = median(hc_exponents)';
top_diff = power;
top_diff.exponent = mean(pa_exponents)' - mean (hc_exponents)';

figure;
cfg = [];
cfg.xlim = [1, 1]; % Fake dimension (frequency) just for the plot
cfg.zlim = [1 1.5];
cfg.parameter = 'exponent';
cfg.comment = 'no';
ft_topoplotER(cfg,top_pa);
colorbar;
title('Patients');
saveas(gca,fullfile(figures_path,'exp_patients.jpeg'));

figure;
ft_topoplotER(cfg,top_hc);
colorbar;
title('Healthy');
saveas(gca,fullfile(figures_path,'exp_healthy.jpeg'));

figure;
cfg.zlim = [];
% cfg.zlim = [-0.06 0.06];
% cfg.colormap = colormap('turbo');
ft_topoplotER(cfg,top_diff);
colorbar;
title('PA - HC');
saveas(gca,fullfile(figures_path,'exp_diff.jpeg'));


%% Cluster-based permutatuion statistics (TO DO)


% Load one dataset to get the electrode layout
data = load_preprocessed_data(params,bidsID);

% Topoplots for both groups
layout = ft_prepare_layout([],data);
figure;
ft_plot_layout(layout);

figure;
ft_plot_sens(data.elec,'label','yes')
ft_plot_axes(data.elec)



% Define neighbours of electrodes (based on Son's code)
cfg = [];
cfg.method = 'distance';
cfg.layout = layout;
cfg.neighbourdist = 0.18; % this value is taken from Son's code (18 mm)
neighbours = ft_prepare_neighbours(cfg, data.elec);





% From Son
% computes a single cluster-based permutation test between the power values of 2 groups
% minimal number of channels for a cluster to be formed: 2
% cluster alpha = 0.05
% statistical test alpha = 0.025 (two-sided)
% based on an independent sample T-test
% 
%   input:  - PA: cell array of FieldTrip pow structures for 1st group (patients)
%           - HC: cell array of FieldTrip pow structures for 2nd group (healthy control)
%           - params: structure containing all relevant parameters
%           - startFreq: lower frequency limit
%           - endFreq: upper frequency limit
%           - relAbs: 'absolute' or 'relative', flag defining whether to
%           use the absolute or relative power
%   output: - stat: FieldTrip structure containing the cluster statistics
% 
% Son Ta Dinh, TU Munich, 22.09.2017
function stat = computeCluster(PA, HC, params, startFreq, endFreq, relAbs)
cfg = [];
cfg.method = 'montecarlo';
if strcmp(relAbs, 'absolute')
    cfg.parameter = 'powspctrm';
elseif strcmp(relAbs, 'relative')
    cfg.parameter = 'relPowspctrm';
else
    disp('Default: using absolute power values')
end
cfg.correctm = 'cluster';
cfg.minnbchan = 2;
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.alpha = 0.05/4;
cfg.tail = 0;
cfg.correcttail = 'prob';
cfg.clustertail = 0;
cfg.numrandomization = params.NREPS;
cfg.statistic = 'indepsamplesT';
design = [ones(1,length(PA)), ones(1,length(HC)) + 1];
cfg.neighbours = params.neighbours;
cfg.design = design;
cfg.ivar = 1;
cfg.frequency = [startFreq endFreq];
cfg.avgoverfreq = 'yes';

stat = ft_freqstatistics(cfg, PA{:}, HC{:});
end