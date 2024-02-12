% Generate specifications and related data for the 1/f project
%
% Cristina Gil Avila, TUM, 15.9.2023

clear all, close all;
% Settings
% Add fieldtrip and analysis functions
addpath('/rechenmagd4/toolboxes_and_functions/fieldtrip');
ft_defaults;
addpath('analysis_functions');

% Load parameter files and define paths
load('../results/features/params.mat');
results_path = '../results/sca/power';
if ~exist(results_path,'dir')
    mkdir(results_path)
end

% Epoch lengths
epoch_length =  {'2', '5'};

%% Load subject ids
participants = readtable(fullfile(params.RawDataPath,'participants_clean.tsv'),'Filetype','text');

% Order the participants.tsv in descending order by bidsID
participants.group = categorical(participants.group);
id = cellfun(@(x) str2double(x(5:7)),participants.participant_id,'UniformOutput',false);
id = cell2mat(id);
[~,ix] = sort(id);
participants_sorted = participants(ix,:);
participant_id = participants_sorted.participant_id;
nSubj = height(participants);


%% Loop over subjects (264)
for iSubj=1:nSubj

    bidsID = participant_id{iSubj};
    bidsID = [bidsID '_task-closed'];

    for iEpoch=1:length(epoch_length)
        ep = epoch_length{iEpoch};

        try
            % Load EEG preprocessed data
            data = load_preprocessed_data(params,bidsID);
        catch
            continue
        end

        % ---- Cut the data into epochs and normalize time axis of the data
        switch ep
            case '2'
                params.EpochLength = 2;
                ep_field = 's2';
            case '5'
                params.EpochLength = 5;
                ep_field = 's5';
        end
        data = epoch_data(params,data);
        temptime = data.time{1};
        [data.time{:}] = deal(temptime);

        % ----- Compute source reconstruction ------
        source = compute_spatial_filter(params,data,'fullSpectrum');

        % ----- Extract virtual channel data -----
        cfg = [];
        cfg.parcellation = 'ROI';
        vdata = ft_virtualchannel(cfg,data,source,params.parcellation);

        % ----- Estimate power spectra at the source level in the PFC only -----
        cfg = [];
        cfg.foilim = [1 100];
        cfg.channel = find(params.PFC_mask);
        cfg.method = 'mtmfft';
        cfg.taper = 'dpss';
        cfg.tapsmofrq = 1;
        cfg.output = 'pow';
        cfg.keeptrials ='no';
        power_dpss.(ep_field){iSubj} = ft_freqanalysis(cfg, vdata);


        cfg = [];
        cfg.foilim = [1 100];
        cfg.channel = find(params.PFC_mask);
        cfg.method = 'mtmfft';
        cfg.taper = 'hanning';
        cfg.output = 'pow';
        cfg.keeptrials ='no';
        power_hanning.(ep_field){iSubj} = ft_freqanalysis(cfg, vdata);

    end

end
power_2s_dpss = power_dpss.s2;
power_2s_hanning = power_hanning.s2;
power_5s_dpss = power_dpss.s5;
power_5s_hanning = power_hanning.s5;
save(fullfile(results_path,'power_sca.mat'),'power_2s_dpss','power_2s_hanning','power_5s_dpss','power_5s_hanning');