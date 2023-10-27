% Generate specifications and related data for the 1/f project
%
% Cristina Gil Avila, TUM, 15.9.2023

clear all, close all;

%% Generate all the specifications
taper = {'hanning', 'dpss'};
zero_padding =  {'0', '5'}; % Optional include 10
average_psd = {'no', 'yes'};
fooof_range = {'2-40', '1-100', '40-60'};
fooof_knee =  {'no', 'yes'};

[a,b,c,d,e] = ndgrid(taper,zero_padding,average_psd,fooof_range,fooof_knee);
nSpec = length(a(:));
specs_id = num2cell(1:nSpec);
specs = [specs_id',a(:),b(:),c(:),d(:),e(:)];
s = cell2table(specs,'VariableNames',{'spec_id','taper','zero_padding','average_psd','fooof_range','fooof_knee'});
clear a b c d e specs_id taper zero_padding average_psd fooof_range fooof_knee
%% Settings
% Define the number of cores for parallelization
% params.Ncores = 2;
% if(isempty(gcp('nocreate')))
%     parObj = parpool(params.Ncores);
% end

% Add fieldtrip and analysis functions
addpath('C:\Users\Mitarbeiter\fieldtrip');
ft_defaults;
addpath('analysis_functions');
addpath('fooof_matlab');
load('../results/features/params.mat');

%%

% Load all subject ids
participants = readtable(fullfile(params.RawDataPath,'participants_rand.tsv'),'Filetype','text');
participant_id = participants.participant_id;
nSubj = height(participants);

%% Loop over all specifications
% Hardcoded
participant_id = {'sub-002'};
nSubj = length(participant_id);
avg_exp = nan(nSpec,nSubj);
t01 = tic;
for iSpec=1:nSpec
    t1 = tic;
    % Loop over subjects
    for iSubj=1:nSubj

        bidsID = participant_id{iSubj};
        bidsID = [bidsID '_task-closed'];
        
        % ----- Load precomputed virtual channel data ------
        load(fullfile(params.VdataPath,[bidsID '_vdata.mat']));

        % ----- Estimate power spectra at the source level in the PFC only -----
        cfg = [];
        cfg.foilim = [1 100];
        cfg.channel = find(params.PFC_mask);
        cfg.method = 'mtmfft';
        % S1. Taper
        cfg.taper = s.taper{iSpec};   
        switch s.taper{iSpec}
            case 'dpss'
                cfg.tapsmofrq = 1;
        end
        % S2. Padding
        switch s.zero_padding{iSpec}
            case '0'
                % Do not pad, don't do anything here
            otherwise
                cfg.pad = str2double(s.zero_padding{iSpec}); 
                cfg.padtype = 'zero';   
        end       
        cfg.output = 'pow';
        cfg.keeptrials ='no';
        power = ft_freqanalysis(cfg, vdata);
         
        % ----- Model power spectrum with FOOOF ------
        % Initialize a fooof object with settings depending on the specification
        fooof_range = strsplit(s.fooof_range{iSpec},'-');
        fooof_range = cellfun(@str2num,fooof_range);
        
        switch s.average_psd{iSpec}
            % S3. Average PSD over channels
            case 'yes'
                pow = mean(power.powspctrm,1);
                
                % Initialize a fooof object handling one channel
                % S4. Fooof range
                % S5. Knee parameter
                switch (s.fooof_knee{iSpec})
                    case 'no'
                        fm = fooof('freq_range',fooof_range);
                    case 'yes'
                        fm = fooof('freq_range',fooof_range,'aperiodic_mode','knee');
                end
            % S3. Average aperiodic exponents    
            case 'no'
                pow = power.powspctrm;
                
                % Initialize a fooof object handling several channels
                % S4. Fooof range
                % S5. Knee parameter
                switch (s.fooof_knee{iSpec})
                    case 'no'
                        fm = fooofGroup('freq_range',fooof_range);
                    case 'yes'
                        fm = fooofGroup('freq_range',fooof_range,'aperiodic_mode','knee');
                end
        end
        
        % Add data to the fooof model (in the predefined freq-range)
        fm = fm.add_data(power.freq,pow); 
        % Fit fooof
        fm = fit(fm);
        
        % S3. Average PSD / exponents
        switch s.average_psd{iSpec}
            case 'yes'
                avg_exp(iSpec,iSubj) = fm.aperiodic_params(end);
            case 'no'
                % Average exponents
                exps = cellfun(@(x) x.aperiodic_params(end), fm.group_results);
                avg_exp(iSpec,iSubj) = mean(exps);                
        end
        % Plot power spectrum
        plot_fm_specs(fm,specs(iSpec,:),'fig_save',true,'file_name',[bidsID '_spec'],'file_path','C:\Users\Mitarbeiter\1overf\results\figures');
    end
    
    t2 = toc(t1);
    fprintf('Specification %d took %.2f seconds \n',[iSpec t2]);

end
t02 = toc(t01);
fprintf('The analysis took %.2f seconds \n',t02);
