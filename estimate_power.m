function estimate_power(params)

participants = readtable(fullfile(params.RawDataPath,'participants_rand.tsv'),'Filetype','text');
task = 'closed';
n = height(participants);
 
participant_id = participants.participant_id;
for iSubj=1:n
    
    bidsID = participant_id{iSubj};
    bidsID = [bidsID '_task-' task];

    % Load EEG preprocessed data
    data = load_preprocessed_data(params,bidsID);

    % Compute source reconstruction
    source = compute_spatial_filter(params,data,'fullSpectrum');
    
    % Extract virtual channel data
    cfg = [];
    cfg.parcellation = 'ROI';
    vdata_trials = ft_virtualchannel(cfg,data,source,params.parcellation);
    parsave_vdata(fullfile(params.VdataPath, [bidsID '_vdata.mat']),vdata_trials);

    % Estimate power spectra at the source level
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

    % Extract power at the PFC
    cfg = [];
    cfg.channel = find(params.PFC_mask);
    power_PFC = ft_selectdata(cfg,power);
    freq = power_PFC.freq;
    pow = power_PFC.powspctrm;
    parsave_pow(fullfile(params.PowerPFCPath, [bidsID '_PFC.mat']),freq,pow);
    
end
end

function parsave_vdata(fname,vdata_trials)

    save(fname,'vdata_trials');

end

function parsave_pow(fname,freq,pow)

    save(fname,'freq','pow');

end

