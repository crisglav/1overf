% Generate specifications and related data for the 1/f project
%
% Cristina Gil Avila, TUM, 15.9.2023

clear all, close all;

% Settings
% Add fieldtrip and analysis functions
addpath('/rechenmagd4/toolboxes_and_functions/fieldtrip');
ft_defaults;
addpath('analysis_functions');
addpath('fooof_matlab');
run('../toolboxes/bayes_factor/installBayesFactor.m')

load('../results/features/params.mat');
figures_path ='../results/figures/';
results_path = '../results/sca/';

% Define the number of cores for parallelization
params.Ncores = 15;
if(isempty(gcp('nocreate')))
    parObj = parpool(params.Ncores);
end
%% Generate all the specifications
if(~exist(fullfile(results_path,'specifications.txt'),'file'))
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
    orig_spec = 8; % Hard coded. This is the specification of the analysis in the pre-registration
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

% % Hardcoded
% participant_id = {'sub-001','sub-002'};
% nSubj = length(participant_id);

% Preallocate variables
exp = nan(nSpec,nSubj);
r_squared = nan(nSpec,nSubj);
mae = nan(nSpec,nSubj);

bayesfactor = nan(nSpec,1);
p = nan(nSpec,1);
effect = nan(nSpec,1);
ci_inf = nan(nSpec,1);
ci_sup = nan(nSpec,1);

% Open log file
fid = fopen(fullfile(results_path,'logfile.txt'),'a');

%% Loop over all specifications
% Take out some variable for improved performance
VdataPath = params.VdataPath;
PFC_mask = params.PFC_mask;

t01 = tic;
for iSpec=1:nSpec
    t1 = tic;
    
    % Specifications
    taper = s.taper{iSpec};
    zero_padding = s.zero_padding{iSpec};
    average_psd = s.average_psd{iSpec};
    fooof_range_s = s.fooof_range{iSpec};
    fooof_knee = s.fooof_knee{iSpec};
    
    % Loop over subjects
    for iSubj=1:nSubj

        bidsID = participant_id{iSubj};
        bidsID = [bidsID '_task-closed'];
        
        % ----- Load precomputed virtual channel data ------
        aux = load(fullfile(VdataPath,[bidsID '_vdata.mat']));

        % ----- Estimate power spectra at the source level in the PFC only -----
        cfg = [];
        cfg.foilim = [1 100];
        cfg.channel = find(PFC_mask);
        cfg.method = 'mtmfft';
        % S1. Taper
        cfg.taper = taper;   
        switch taper
            case 'dpss'
                cfg.tapsmofrq = 1;
        end
        % S2. Padding
        switch zero_padding
            case '0'
                % Do not pad, don't do anything here
            otherwise
                cfg.pad = str2double(zero_padding); 
                cfg.padtype = 'zero';   
        end       
        cfg.output = 'pow';
        cfg.keeptrials ='no';
        power = ft_freqanalysis(cfg, aux.vdata);
         
        % ----- Model power spectrum with FOOOF ------
        % Initialize a fooof object with settings depending on the specification
        fooof_range = strsplit(fooof_range_s,'-');
        fooof_range = cellfun(@str2num,fooof_range);
        
        switch average_psd
            % S3. Average PSD over channels
            case 'yes'
                pow = mean(power.powspctrm,1);
                
                % Initialize a fooof object handling one channel
                % S4. Fooof range
                % S5. Knee parameter
                switch (fooof_knee)
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
                switch (fooof_knee)
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

%         % Plot power spectrum
%         plot_fm_specs(fm,s(iSpec,:),'loglog',false,'fig_save',false,'file_name',[bidsID '_spec'],'file_path',figures_path);
%         close;
        
    end
    % Statistical test: differences in aperiodic componentbetween groups
    [bayesfactor(iSpec),p(iSpec)] = bf.ttest2(exp(iSpec,hc_mask),exp(iSpec,pa_mask));
    % Effect size (hedge's g and confidence interval (works from ML22a onwards)
    d = meanEffectSize(exp(iSpec,hc_mask),exp(iSpec,pa_mask),Effect='cohen',ConfidenceIntervalType='bootstrap',NumBootstraps=1000);%    Works in ML2022a onwards
    effect(iSpec) = d.Effect;
    ci_inf(iSpec) = d.ConfidenceIntervals(1);
    ci_sup(iSpec) = d.ConfidenceIntervals(2);

    % Cohen's d by hand
%     n1 = sum(hc_mask);
%     n2 = sum(pa_mask);
%     m1 = mean(exp(iSpec,hc_mask));
%     m2 = mean(exp(iSpec,pa_mask));
%     s1 = nanvar(exp(iSpec,hc_mask))*(n1-1);
%     s2 = nanvar(exp(iSpec,pa_mask))*(n2-1);
%     pooledsd = sqrt((s1 + s2)/(n1 + n2 - 2));
%     effect(iSpec) = (m1 - m2)/pooledsd;
    
    t2 = toc(t1);
    fprintf(fid,'Specification %d took %.2f seconds \n',[iSpec t2]);
    
end
t02 = toc(t01);
fprintf(fid,'The analysis took %.2f minutes \n',t02/60);

% Save results
results = s;
results.bayes_factor = bayesfactor;
results.p_value = p;
results.effect_size = effect;
results.ci_inf = ci_inf;
results.ci_sup = ci_sup;
writetable(results,fullfile(results_path,'specs_bf.txt'),'Delimiter',',');
save(fullfile(results_path,'specs_ap_exp.m'),'exp','r_squared','mae');
%% Plot all exponents gathered from all the specifications for one subject
c = lines(2);
figure();
scatter(1:nSpec,exp(:,1)',20,c(1,:));
hold on,
scatter(orig_spec,exp(orig_spec,1),20,'k','filled');
hold on,
plot(1:nSpec,r_squared(:,1),'Color',c(1,:))
hold on,
% Subject 2
scatter(1:nSpec,exp(:,2)',20,c(2,:));
hold on,
scatter(orig_spec,exp(orig_spec,2),20,'k','filled');
hold on,
plot(1:nSpec,r_squared(:,2),'Color',c(2,:))
