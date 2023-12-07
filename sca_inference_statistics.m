% Inferential statistics for the specifiction curve analysis
%
% Based on Elisabeth's May script
% /rechenmagd4/Experiments/2020_10_alpha_peak_frequency/code/matlab/functions/APF_RQ2_S5_inferentialCurve_statistics.m
% 
% Cristina Gil Avila, 08.11.2023

% define and create output folder
datapath = '/rechenmagd3/Experiments/2023_1overf/results_v4/sca';
figures_folder = '../results_v3/figures/sca_inference/';
if ~exist(figures_folder,'dir')
    mkdir(figures_folder);
end
% Hardcoded number of randomizations and specifications
nRand = 500;
nSpec = 48;

% LOAD ORIGINAL CURVE
results_orig = readtable(fullfile(datapath,'specs_bf.txt'));
[d_orig, ix] = sort(results_orig.effect_size);
bf_temp = results_orig.bayes_factor;
bf_orig = bf_temp(ix); % Sort according the effect size
pvalues_temp = results_orig.p_value;
pvalues_orig = pvalues_temp(ix);

% Predominant direction of effect for the original data
pos = sum(d_orig > 0);
neg = sum(d_orig < 0);
if pos >= neg
    sign_orig = 1;
elseif neg > pos
    sign_orig = -1;
end

% LOAD RANDOMIZED SPECIFICTIONS 
d_rand = nan(nSpec,nRand);
bf_rand = nan(nSpec,nRand);
pvalues_rand = nan(nSpec,nRand);
sign = nan(1,nRand);
for iRand=1:nRand
    % load statistics of randomizations
    results = readtable(fullfile(datapath,sprintf('specs_bf_rand%.3d.txt',iRand)));
    [d_rand(:,iRand),ix] = sort(results.effect_size);
    bf_temp = results.bayes_factor;
    bf_rand(:,iRand) = bf_temp(ix);
    pvalues_temp = results.p_value;
    pvalues_rand(:,iRand) = pvalues_temp(ix);
    
    % Calculate the dominant sign of effects for each randomization
    pos = sum(d_rand(:,iRand) > 0);
    if (pos/nSpec) > 0.5
        sign(iRand) = 1;
    elseif (pos/nSpec) < 0.5
        sign(iRand) = -1;
    else
        sign(iRand) = sign_orig; % if there is no dominant sign of effects, assume direction of original effect
    end
end
% Flip the curves that do not have the dominant sign. The dominant sign is
% based on the original curve.
d_2s = d_rand.*sign*sign_orig;

%% LOAD ELISABETHS DATA
nRand = 500;
nSpec = 72;

% Load original curve
datapath = '/rechenmagd4/Experiments/2020_10_alpha_peak_frequency/results/stats/';
results_orig = readtable(fullfile(datapath,'results_RQ2_session1_all_specs_orig_stats.csv'));
rowsElec = find(results_orig.electrodes == "global");
[d_orig, ix]= sort(results_orig.R(rowsElec,:));
bf_temp = results_orig.BF10(rowsElec,:);
bf_orig = bf_temp(ix);
% Predominant direction of effect for the original data
pos = sum(d_orig > 0);
neg = sum(d_orig < 0);
if pos >= neg
    sign_orig = 1;
elseif neg > pos
    sign_orig = -1;
end

% LOAD RANDOMIZED SPECIFICTIONS
% create some variables for collection of results
d_rand  =  nan(nSpec,nRand);
bf_rand    = nan(nSpec,nRand);

for iRand = 1:nRand

    % load statistics of randomizations
    results = readtable([datapath 'results_RQ2_session1_all_specs_rand' num2str(iRand) '_stats.csv']);
    rowsElec = find(results.electrodes == "global");
    % sort according to effect size
    [d_rand(:,iRand), ix] = sort(results.R(rowsElec,:));
    % again, sort Bayes Factors in same way to keep assignment
    bf_temp = results.BF10(rowsElec,:);
    bf_rand(:,iRand)   = bf_temp(ix);

    % calculate the dominant sign of effects for each randomization
    pos = sum(d_rand(:,iRand) > 0);
    if (pos/nSpec) > 0.5
        sign(iRand) = 1;
    elseif (pos/nSpec) < 0.5
        sign(iRand) = -1;
    else
        sign(iRand) = sign_orig; % if there is no dominant sign of effects, assume direction of original effect
    end
end
d_2s = d_rand.*sign*sign_orig;
%% INFERENTIAL TESTS

% 1. TEST OF MEDIAN EFFECT SIZE
% ===================================================================

% One-sided test
% ==============
% calculate p-value as percentage of randomizations with median
% effects size larger/smaller than original effect size - direction
% depends on sign of the original effect size
median_d_orig = median(d_orig);
median_d_rand = median(d_rand);
if sign_orig == 1
    t = sum(median_d_rand >= median_d_orig);
elseif sign_orig == -1
    t = sum(median_d_rand <= median_d_orig);
end
p_median_1s = t/nRand;

% Two-sided test
% ==============
% calculate p-value as percentage of randomizations with median
% effects size more extreme than original effect size on both ends
% of the distribution of median effect sizes
t = sum(abs(median_d_rand) >= abs(median_d_orig));
p_median_2s = t/nRand;

% Test flipping the curves (2 sided)
% ==============
% This test is conceptually the same as the two-sided test.
median_d_rand_2s = median(d_2s);
if sign_orig == 1
    t = sum(median_d_rand_2s >= median_d_orig);
elseif sign_orig == -1 
    t = sum(median_d_rand_2s <= median_d_orig);
end
p_median_flip = t/nRand;


% 2. TEST OF SHARE OF 'SIGNIFICANT' RESULTS (HERE SIGNIFICANT: BF > 3)
% ====================================================================
% Extract the effects (d) of the 'significant' specifications
significant_d_orig = d_orig(bf_orig > 3);
% Count only the significant specifications that have the same sign as the
% dominant curve effect
if ~isempty(significant_d_orig)
    BF_orig_count = sum((significant_d_orig * sign_orig) > 0);
else
    BF_orig_count = 0;
end
% Extract the effects of 'significant' specifications for all
% randomizations
significant_d_rand = nan(nSpec,nRand); 
significant_d_rand(bf_rand > 3) = d_rand(bf_rand > 3);

% One-sided test
% ==============
% Count the significant specifications that have the same sign
% as the dominant curve effect from the original curve (!)
BF_rand_count = sum((significant_d_rand * sign_orig) > 0,1);
t = sum(BF_rand_count>=BF_orig_count); 
p_share_1s = t/nRand;

% Two-sided test
% ==============
% For each randomization, count the significant specifications
% that have the same sign as the dominant curve effect of that
% randomization (!)
BF_rand_count = sum((significant_d_rand .* sign) > 0,1);
t = sum(BF_rand_count>=BF_orig_count); 
p_share_2s = t/nRand;

% Test flipping the curves (2 sided)
% ==============
% This test is conceptually the same as the two-sided test.
% Calculate the dominant sign based on the flipped curves (edom)
edom = nan(nSpec,nRand);
if sign_orig == 1
    edom(d_2s >= 0) = 1;
    edom(d_2s < 0) = 0;
else
    edom(d_2s >= 0) = 0;
    edom(d_2s < 0) = 1;
end
sig_freq = bf_rand > 3;
sig_dom = and(sig_freq,edom);
sig_dom_freq = sum(sig_dom);
t = sum(sig_dom_freq >= BF_orig_count);
p_share_flip = t/nRand;


% 3. TEST OF AGGREGATED P-VALUES
% =================================================================
% For each randomization, average pvalues following the Stouffer's method.
% First get z-scores for each pvalue and average across specifications
% Pvalues are diveded by two because the original t-tests between patients
% and healthy participants were two-sided.
z_temp = norminv(pvalues_orig./2);
z_orig = sum(z_temp)/sqrt(nSpec);

z_temp = norminv(pvalues_rand./2);
z_rand = sum(z_temp,1)./sqrt(nSpec);

% One-sided test
% ==============
% calculate p-value as percentage of randomizations with average z value
% larger/smaller than original average z value - direction
% depends on sign of the original z value
if z_orig > 0
    t = sum(z_rand >= z_orig);
elseif z_orig < 0
    t = sum(z_rand <= z_orig);
else
    t = nan;
end
p_share_1s = t/nRand;

% Two-sided test
% ==============
% calculate p-value as percentage of randomizations with average z value 
% more extreme than original average z value on both ends of the distribution of z
% values rand
t = sum(abs(z_rand) >= abs(z_orig));
p_share_2s = t/nRand;

% Test flipping the curves (2 sided)
% ==============
pvalues_2s = pvalues_rand.*sign;
z_temp = norminv(pvalues_2s./2);
z_rand_2s = sum(z_temp,1)./sqrt(nSpec);

if sign_orig == 1
    t = sum(z_rand_2s >= z_orig);
elseif sign_orig == -1 
    t = sum(z_rand_2s <= z_orig);
end
p_share_flip = t/nRand;



%% PLOT MEDIAN CURVES
figure;
grey = [0.8 0.8 0.8];

tiledlayout(3,3)

% ALL THE RANDOMIZED CURVES
nexttile
plot(1:nSpec,d_rand), hold on;
plot(1:nSpec,d_orig, 'Color', 'k', 'LineWidth', 2); hold on;
scatter(24,d_orig(24),'black','filled'); hold on;
ylim([-0.5,0.5])
yline(0);
ylabel('d');
title('All randomizations')

nexttile
plot(1:nSpec,d_2s), hold on;
plot(1:nSpec,d_orig, 'Color', 'k', 'LineWidth', 2); hold on;
scatter(24,d_orig(24),'black','filled'); hold on;
ylim([-0.5,0.5])
yline(0);
title('All randomizations, flipped sign')

nexttile
plot(1:nSpec,sort(d_2s)), hold on;
plot(1:nSpec,d_orig, 'Color','k', 'LineWidth', 2); hold on;
scatter(24,d_orig(24),'black','filled'); hold on;
ylim([-0.5,0.5])
yline(0);
title('All randomizations, flipped sign, sorted')

% MEDIAN TESTS
% Plot the median ttest Cristina's verstion
nexttile
histogram(median(d_rand),'BinWidth',0.02);
hold on
x1 = xline(median(d_orig),'k','LineWidth',2);
hold on
x2 = xline(prctile(median(d_rand),2.5),'r');
hold on
xline(prctile(median(d_rand),97.5),'r')
xlabel('Median d')
legend([x1, x2],'observed d','tcrit 2.5%')
text(min(xlim), 17, ['p = ' num2str(round(p_median_cris,2))]);
box off;
title('Median test');

% Plot the median ttest original version
nexttile;
histogram(median(d_2s),'BinWidth',0.02);
hold on
x1 = xline(median(d_orig),'k','LineWidth',2);
hold on
xline(prctile(median(d_2s),95),'r')
xlabel('Median d')
legend([x1, x2],'observed d','tcrit 5%')
text(min(xlim), 31, ['p = ' num2str(round(p_median,2))]);
box off;
title('Median test, flipped');

% For completeness the last option
nexttile;
histogram(median(sort(d_2s)),'BinWidth',0.02);
hold on
x1 = xline(median(d_orig),'k','LineWidth',2);
hold on
xline(prctile(median(sort(d_2s)),95),'r')
xlabel('Median d')
legend([x1, x2],'observed d','tcrit 5%')
text(min(xlim), 31, ['p = ' num2str(round(p_median_sorted,2))]);
box off;
title('Median test, flipped, sorted');


% PLOT MEDIAN CURVES
% Original curve with median and prct
nexttile
plot(1:nSpec,sort(prctile(d_rand,50,2)),'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(prctile(d_rand,2.5,2)), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(prctile(d_rand,97.5,2)), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,d_orig, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); hold on;
scatter(24,d_orig(24),'black','filled'); hold on;
ylim([-0.5,0.5])
yline(0);

nexttile
plot(1:nSpec,sort(prctile(d_2s,50,2)),'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(prctile(d_2s,2.5,2)), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
% plot(1:nSpec,sort(prctile(d_2s,95,2)), 'Color', 'red', 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(prctile(d_2s,97.5,2)), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,d_orig, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); hold on;
scatter(24,d_orig(24),'black','filled'); hold on;
ylim([-0.5,0.5])
yline(0);

nexttile
plot(1:nSpec,prctile(sort(d_2s),50,2),'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,prctile(sort(d_2s),2.5,2), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
% plot(1:nSpec,prctile(sort(d_2s),95,2), 'Color', 'red', 'LineWidth', 1.5), hold on;
plot(1:nSpec,prctile(sort(d_2s),97.5,2), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,d_orig, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); hold on;
scatter(24,d_orig(24),'black','filled'); hold on;
ylim([-0.5,0.5])
yline(0);
title('Paper plot');



% They are not the same!!
% Median value of the 95 percentile across all randomizations
med_per95 = median(prctile(d_rand,95,2))
% Percentile 95 of the median across all randomizations
per95_med = prctile(median(d_rand),95)

% They are not the same!!
% Median value of the 95 percentile across all randomizations
med_per95 = median(prctile(d_2s,95,2))
% Percentile 95 of the median across all randomizations
per95_med = prctile(median(d_2s),95)

med_per95 = median(prctile(sort(d_2s),95,2))
% Percentile 95 of the median across all randomizations
per95_med = prctile(median(sort(d_2s)),95)
median(d_2s) == median(sort(d_2s));
% But
prctile(d_2s,95,2) ~= prctile(sort(d_2s),95,2);
median(d_2s,2) ~= median(sort(d_2s),2)

%% Figures for internal presentation
figure;
tiledlayout(1,2);
nexttile;
plot(1:nSpec,d_orig, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); hold on;
plot(1:nSpec,sort(prctile(d_rand,50,2)),'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(prctile(d_rand,2.5,2)), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(prctile(d_rand,97.5,2)), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
scatter(24,d_orig(24),'black','filled'); hold on;
ylim([-0.5,0.5])
yline(0);
title('All randomizations');
ylabel('Effect size (Cohens d)')
xlabel('Specifications sorted by effect size')
nexttile;
plot(1:nSpec,d_orig, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); hold on;
plot(1:nSpec,prctile(sort(d_2s),50,2),'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,prctile(sort(d_2s),2.5,2), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,prctile(sort(d_2s),97.5,2), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
scatter(24,d_orig(24),'black','filled'); hold on;
ylim([-0.5,0.5])
yline(0);
legend({'Observed data','Median under the null','2.5th and 97.5th under the null'})
title('All randomizations, flipped sign, sorted')