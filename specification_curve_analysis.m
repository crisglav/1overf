% Generate specifications and related data for the 1/f project
%
% Cristina Gil Avila, TUM, 15.9.2023

clear all, close all;
rng('default');
% Settings
% Add fieldtrip and analysis functions
addpath('/rechenmagd4/toolboxes_and_functions/fieldtrip');
ft_defaults;
addpath('analysis_functions');
addpath('fooof_matlab');
run('../toolboxes/bayes_factor/installBayesFactor.m')

% Load parameter files and define paths
load('../results/features/params.mat');
figures_path ='../results/figures/';
results_path = '../results/sca/';
if ~exist(results_path,'dir')
    mkdir(results_path)
end
% Define the number of cores for parallelization
params.Ncores = 20;
if(isempty(gcp('nocreate')))
    parObj = parpool(params.Ncores);
end
%% Generate all the specifications
if(~exist(fullfile(results_path,'specifications.txt'),'file'))
    epoch_length =  {'2', '5'};
    taper = {'dpss','hanning'};
    average_psd = {'yes','no'};
    fooof_range = {'2-40', '1-100', '40-60'};
    fooof_knee =  {'no', 'yes'};

    [a,b,c,d,e] = ndgrid(epoch_length,taper,average_psd,fooof_range,fooof_knee);
    nSpec = length(a(:));
    specs_id = num2cell(1:nSpec);
    specs = [specs_id',a(:),b(:),c(:),d(:),e(:)];
    s = cell2table(specs,'VariableNames',{'spec_id','epoch_length','taper','average_psd','fooof_range','fooof_knee'});
    orig_spec = 1; % Hard coded. This is the specification of the main analysis
    clear a b c d e specs_id taper zero_padding average_psd fooof_range fooof_knee
    writetable(s,fullfile(results_path,'specifications.txt'),'Delimiter',',');
else
    s = readtable(fullfile(results_path,'specifications.txt'),'Format','%d%s%s%s%s%s');
    nSpec = height(s);
end
%% Load subject ids
participants = readtable(fullfile(params.RawDataPath,'participants_clean.tsv'),'Filetype','text');

% Order the participants.tsv in descending order by bidsID
participants.group = categorical(participants.group);
id = cellfun(@(x) str2double(x(5:7)),participants.participant_id,'UniformOutput',false);
id = cell2mat(id);
[~,ix] = sort(id);
participants_sorted = participants(ix,:);

% Separate patients and healthy participants groups
pa_mask = participants_sorted.group == 'pa';
hc_mask = participants_sorted.group == 'hc';

participant_id = participants_sorted.participant_id;
nSubj = height(participants);

% Load pre-computed PSDs
load('../results/sca/power/power_sca.mat');

% list executed randomizations, extract randomizations ids
rfiles = dir(fullfile(results_path,'specs_bf_rand*.txt'));
ids = str2double(extract({rfiles.name},digitsPattern(3)));
start_rand = max(ids)+1;
if isempty(start_rand)
    start_rand = 0;
end
% Open log file
fid = fopen(fullfile(results_path,'logfile.txt'),'a');

%% Loop over all randomizations
nRand = 500;

t001 = tic;
for iRand=start_rand:nRand
    t01 = tic;
    % Preallocate variables
    exp = nan(nSpec,nSubj);
    r_squared = nan(nSpec,nSubj);
    mae = nan(nSpec,nSubj);
    succeed = zeros(nSpec,nSubj);
    
    bayesfactor = nan(nSpec,1);
    p = nan(nSpec,1);
    effect = nan(nSpec,1);
    ci_inf = nan(nSpec,1);
    ci_sup = nan(nSpec,1);

    % Randomization 0 is the original specification curve
    if iRand>0
        ix = randperm(nSubj);
        rand_group = participants_sorted.group(ix);
        participants_sorted.group = rand_group;
        % Separate patients and healthy participants groups
        pa_mask = participants_sorted.group == 'pa';
        hc_mask = participants_sorted.group == 'hc';
    end

    % Loop over all specifications (48)
    for iSpec=1:nSpec
        t1 = tic;

        % Specifications
        epoch_length = s.epoch_length{iSpec};
        taper = s.taper{iSpec};
        average_psd = s.average_psd{iSpec};
        fooof_range_s = s.fooof_range{iSpec};
        fooof_knee = s.fooof_knee{iSpec};

        % S1. Epoch length
        % S2. Taper
        if and(strcmp(epoch_length,'2'),strcmp(taper,'dpss'))
            power_all = power_2s_dpss;
        elseif and(strcmp(epoch_length,'2'),strcmp(taper,'hanning'))
            power_all = power_2s_hanning;
        elseif and(strcmp(epoch_length,'5'),strcmp(taper,'dpss'))
            power_all = power_5s_dpss;
        elseif and(strcmp(epoch_length,'5'),strcmp(taper,'hanning'))
            power_all = power_5s_hanning;
        end
      
        % Loop over subjects (264)
        parfor iSubj=1:nSubj

            power = power_all{iSubj};

            % ----- Model power spectrum with FOOOF ------
            % Initialize a fooof object with settings depending on the specification
            fooof_range = strsplit(fooof_range_s,'-');
            fooof_range = cellfun(@str2num,fooof_range);

            % S3. Average PSD
            switch average_psd
                % S3a. Average PSD over channels
                case 'yes'
                    pow = mean(power.powspctrm,1);

                    % Initialize a fooof object handling one channel
                    % S4. Fooof range
                    % S5. Knee parameter
                    switch (fooof_knee)
                        case 'no'
                            fm = fooof('freq_range',fooof_range,'peak_width_limits',[1 12]);
                        case 'yes'
                            fm = fooof('freq_range',fooof_range,'peak_width_limits',[1 12],'aperiodic_mode','knee');
                    end
                    % S3b. Average aperiodic exponents
                case 'no'
                    pow = power.powspctrm;

                    % Initialize a fooof object handling several channels
                    % S4. Fooof range
                    % S5. Knee parameter
                    switch (fooof_knee)
                        case 'no'
                            fm = fooofGroup('freq_range',fooof_range,'peak_width_limits',[1 12]);
                        case 'yes'
                            fm = fooofGroup('freq_range',fooof_range,'peak_width_limits',[1 12],'aperiodic_mode','knee');
                    end
            end

            % Add data to the fooof model (in the predefined freq-range)
            fm = fm.add_data(power.freq,pow);
            % Fit fooof
            try
                fm = fit(fm);
            catch
                succeed(iSpec,iSubj) = -1;
                continue;
            end

            % S3. Average PSD / exponents
            switch average_psd
                case 'yes'
                    exp(iSpec,iSubj) = fm.aperiodic_params(end);
                case 'no'
                    % Average exponents
                    exps = cellfun(@(x) x.aperiodic_params(end), fm.group_results);
                    exp(iSpec,iSubj) = mean(exps);
            end

            % Extract error model
            if isprop(fm,'group_results')
                r_squared(iSpec,iSubj) = median(cellfun(@(x) x.r_squared, fm.group_results));
                mae(iSpec,iSubj) = median(cellfun(@(x) x.error(1), fm.group_results));
            else
                r_squared(iSpec,iSubj) = fm.r_squared;
                mae(iSpec,iSubj) = fm.error(1);
            end

            % If everything went fine acknowledge it
            succeed(iSpec,iSubj) = 1;

            %         % Plot power spectrum
            %         fig = plot_fm_specs(fm,s(iSpec,:),'loglog',false,'fig_save',false,'file_name',[bidsID '_spec'],'file_path',figures_path);
            %         close;

        end
        % Statistical test: differences in aperiodic component between groups
        [bayesfactor(iSpec),p(iSpec)] = bf.ttest2(exp(iSpec,hc_mask),exp(iSpec,pa_mask));
        % Effect size (hedge's g and confidence interval (works from ML22a onwards)
        d = meanEffectSize(exp(iSpec,hc_mask),exp(iSpec,pa_mask),Effect='cohen',ConfidenceIntervalType='bootstrap',NumBootstraps=1000);%    Works in ML2022a onwards
        effect(iSpec) = d.Effect;
        ci_inf(iSpec) = d.ConfidenceIntervals(1);
        ci_sup(iSpec) = d.ConfidenceIntervals(2);

        t2 = toc(t1);
        fprintf(fid,'Specification %d took %.2f seconds \n',[iSpec t2]);

    end
    t02 = toc(t01);
    fprintf(fid,'The randomization %d took %.2f minutes \n',[iRand t02/60]);

    % Save results
    results = s;
    results.bayes_factor = bayesfactor;
    results.p_value = p;
    results.effect_size = effect;
    results.ci_inf = ci_inf;
    results.ci_sup = ci_sup;
    if iRand ==0
        writetable(results,fullfile(results_path,'specs_bf.txt'),'Delimiter',',');
        save(fullfile(results_path,'specs_ap_exp.mat'),'exp','r_squared','mae','succeed');
    else
        writetable(results,fullfile(results_path,sprintf('specs_bf_rand%.3d.txt',iRand)),'Delimiter',',');
        save(fullfile(results_path,sprintf('specs_ap_exp_rand%.3d.mat',iRand)),'exp','r_squared','mae','succeed');
    end

end
t002 = toc(t001);
fprintf(fid,'Script ended successfully. It took: %.2f hours \n',t002/(60*60));
fclose(fid);
