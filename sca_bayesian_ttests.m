% Load the data from the SCA randomizations, compute bayesian ttests in R
%
% Cristina Gil, Felix Bott, 06.12.2023

% define and create output folder
datapath = '/rechenmagd3/Experiments/2023_1overf/results/sca';
addpath('Rinterface');
run('../toolboxes/bayes_factor/installBayesFactor.m')

% Hardcoded number of randomizations and specifications
nRand = 500;
nSpec = 48;

for iRand=1:nRand
    % load statistics of randomizations
    load(fullfile(datapath,sprintf('specs_ap_exp_rand%.3d.mat',iRand)));
    hc = exp(hc_mask);
    pa = exp(pa_mask);
    
    % Frequentist ttest
    [~,p,ci,stats] = ttest2(hc,pa);
    tstat = stats.tstat;

    % Bayesian ttest in R
    restab = ttestBF2(hc,pa);
    
%     % Bayesian ttest in matlab
%     [b,p] = bf.ttest2(hc,pa);
%     [b,p] = bf.ttest('T',tstat,'N',[length(hc) length(pa)],'tail','both');

    % Effect size (Hedges g)        
    d = meanEffectSize(hc,pa,Effect='cohen',ConfidenceIntervalType='bootstrap',NumBootstraps=1000);%    Works in ML2022a onwards
    effect(iSpec) = d.Effect;
    ci_inf(iSpec) = d.ConfidenceIntervals(1);
    ci_sup(iSpec) = d.ConfidenceIntervals(2);

    results = table();
%     results.bayes_factor = bayesfactor;
    results.p_value = p;
    results.effect_size = effect;
    results.ci_inf = ci_inf;
    results.ci_sup = ci_sup;


end
