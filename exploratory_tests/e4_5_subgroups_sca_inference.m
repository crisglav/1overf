clear all, close all,
addpath('/rechenmagd3/Experiments/2023_1overf/code/')

% define and create output folder
sca_path = '/rechenmagd3/Experiments/2023_1overf/results/sca/';
figures_path = '../../results/figures/sca_inference/';
if ~exist(figures_path,'dir')
    mkdir(figures_path);
end

% Direction of hypothesis has to be specified by the researcher
tail = 'both'; % 'right', 'left'
groups = {'cwp','cbp'};
%% LOAD THE DATA
nSpec = 48;
nRand = 500;

%% HYPOTHESIS 1

%========= CWP ==========
% LOAD ORIGINAL CURVE
results_orig = readtable(fullfile(sca_path,'randomizations_e4_h1','stats_orig.csv'));
[d_orig, ix] = sort(results_orig.d_cwp);
bf_temp = results_orig.BF_cwp;
bf_orig = bf_temp(ix); % Sort according the effect size
pvalues_temp = results_orig.pvalue_cwp;
pvalues_orig = pvalues_temp(ix);

% LOAD RANDOMIZED SPECIFICTIONS 
d_rand = nan(nSpec,nRand);
bf_rand = nan(nSpec,nRand);
pvalues_rand = nan(nSpec,nRand);
for iRand=1:nRand
    % load statistics of randomizations
    results = readtable(fullfile(sca_path,'randomizations_e4_h1',sprintf('stats_rand%.3d.csv',iRand)));
    [d_rand(:,iRand),ix] = sort(results.d_cwp);
    bf_temp = results.BF_cwp;
    bf_rand(:,iRand) = bf_temp(ix);
    pvalues_temp = results.pvalue_cwp;
    pvalues_rand(:,iRand) = pvalues_temp(ix);   
end

% INFERENTIAL TESTS
[p_median,p_share,p_aggregate] = sca_inference_tests(d_orig,bf_orig,pvalues_orig,d_rand,bf_rand,pvalues_rand,tail);
inference.h1.cwp = [p_median,p_share,p_aggregate];

% Save data in csv
lowPC = prctile(d_rand,2.5,2);
highPC = prctile(d_rand,97.5,2);
t = table(d_orig,sort(prctile(d_rand,50,2)),lowPC,highPC,'VariableNames',{'d_orig','d_median','lowPC','highPC'});
writetable(t, fullfile(sca_path,'e4_sca_inference_h1_cwp.csv'));

%========= CBP ==========
% ORIGINAL CURVE
[d_orig, ix] = sort(results_orig.d_cbp);
bf_temp = results_orig.BF_cbp;
bf_orig = bf_temp(ix); % Sort according the effect size
pvalues_temp = results_orig.pvalue_cbp;
pvalues_orig = pvalues_temp(ix);

% LOAD RANDOMIZED SPECIFICTIONS 
d_rand = nan(nSpec,nRand);
bf_rand = nan(nSpec,nRand);
pvalues_rand = nan(nSpec,nRand);
for iRand=1:nRand
    % load statistics of randomizations
    results = readtable(fullfile(sca_path,'randomizations_e4_h1',sprintf('stats_rand%.3d.csv',iRand)));
    [d_rand(:,iRand),ix] = sort(results.d_cbp);
    bf_temp = results.BF_cbp;
    bf_rand(:,iRand) = bf_temp(ix);
    pvalues_temp = results.pvalue_cbp;
    pvalues_rand(:,iRand) = pvalues_temp(ix);   
end

% INFERENTIAL TESTS
[p_median,p_share,p_aggregate] = sca_inference_tests(d_orig,bf_orig,pvalues_orig,d_rand,bf_rand,pvalues_rand,tail);
inference.h1.cbp = [p_median,p_share,p_aggregate];

% Save data in csv
lowPC = prctile(d_rand,2.5,2);
highPC = prctile(d_rand,97.5,2);
t = table(d_orig,sort(prctile(d_rand,50,2)),lowPC,highPC,'VariableNames',{'d_orig','d_median','lowPC','highPC'});
writetable(t, fullfile(sca_path,'e4_sca_inference_h1_cbp.csv'));


%% HYPOTHESIS 2

%========= CWP ==========
% LOAD ORIGINAL CURVE
results_orig = readtable(fullfile(sca_path,'randomizations_e4_h2','stats_orig.csv'));
[d_orig, ix] = sort(results_orig.R_cwp);
bf_temp = results_orig.BF_cwp;
bf_orig = bf_temp(ix); % Sort according the effect size
pvalues_temp = results_orig.pvalue_cwp;
pvalues_orig = pvalues_temp(ix);

% LOAD RANDOMIZED SPECIFICTIONS 
d_rand = nan(nSpec,nRand);
bf_rand = nan(nSpec,nRand);
pvalues_rand = nan(nSpec,nRand);
for iRand=1:nRand
    % load statistics of randomizations
    results = readtable(fullfile(sca_path,'randomizations_e4_h2',sprintf('stats_rand%.3d.csv',iRand)));
    [d_rand(:,iRand),ix] = sort(results.R_cwp);
    bf_temp = results.BF_cwp;
    bf_rand(:,iRand) = bf_temp(ix);
    pvalues_temp = results.pvalue_cwp;
    pvalues_rand(:,iRand) = pvalues_temp(ix);   
end

% INFERENTIAL TESTS
[p_median,p_share,p_aggregate] = sca_inference_tests(d_orig,bf_orig,pvalues_orig,d_rand,bf_rand,pvalues_rand,tail);
inference.h2.cwp = [p_median,p_share,p_aggregate];

% Save data in csv
lowPC = prctile(d_rand,2.5,2);
highPC = prctile(d_rand,97.5,2);
t = table(d_orig,sort(prctile(d_rand,50,2)),lowPC,highPC,'VariableNames',{'R_orig','R_median','lowPC','highPC'});
writetable(t, fullfile(sca_path,'e4_sca_inference_h2_cwp.csv'));

%========= CBP ==========
% LOAD ORIGINAL CURVE
results_orig = readtable(fullfile(sca_path,'randomizations_e4_h2','stats_orig.csv'));
[d_orig, ix] = sort(results_orig.R_cbp);
bf_temp = results_orig.BF_cbp;
bf_orig = bf_temp(ix); % Sort according the effect size
pvalues_temp = results_orig.pvalue_cbp;
pvalues_orig = pvalues_temp(ix);

% LOAD RANDOMIZED SPECIFICTIONS 
d_rand = nan(nSpec,nRand);
bf_rand = nan(nSpec,nRand);
pvalues_rand = nan(nSpec,nRand);
for iRand=1:nRand
    % load statistics of randomizations
    results = readtable(fullfile(sca_path,'randomizations_e4_h2',sprintf('stats_rand%.3d.csv',iRand)));
    [d_rand(:,iRand),ix] = sort(results.R_cbp);
    bf_temp = results.BF_cbp;
    bf_rand(:,iRand) = bf_temp(ix);
    pvalues_temp = results.pvalue_cbp;
    pvalues_rand(:,iRand) = pvalues_temp(ix);   
end

% INFERENTIAL TESTS
[p_median,p_share,p_aggregate] = sca_inference_tests(d_orig,bf_orig,pvalues_orig,d_rand,bf_rand,pvalues_rand,tail);
inference.h2.cbp = [p_median,p_share,p_aggregate];

% Save data in csv
lowPC = prctile(d_rand,2.5,2);
highPC = prctile(d_rand,97.5,2);
t = table(d_orig,sort(prctile(d_rand,50,2)),lowPC,highPC,'VariableNames',{'R_orig','R_median','lowPC','highPC'});
writetable(t, fullfile(sca_path,'e4_sca_inference_h2_cbp.csv'));

save(fullfile(sca_path,'e4_inference.mat'),'inference');
%% PLOT
grey = [0.8 0.8 0.8];
f = figure;
plot(1:nSpec,d_orig, 'Color', '#5400a2', 'LineWidth', 1.5); hold on;
plot(1:nSpec,sort(prctile(d_rand,50,2)),'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(lowPC), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(highPC), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
xlim([1 nSpec]);
% ylim([-0.5,0.5])
yline(0);
title('Null distribution of specification curves');
ylabel('Effect size (R)')
xlabel('Specifications (nr, sorted by effect size)')
set(gca, 'TickDir', 'out');
box off
legend({'Observed data','Median under the null','2.5th and 97.5th under the null'})
