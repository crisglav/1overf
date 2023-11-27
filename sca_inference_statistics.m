% Inferential statistics for the specifiction curve analysis
%
% Based on Elisabeth's May script
% /rechenmagd4/Experiments/2020_10_alpha_peak_frequency/code/matlab/functions/APF_RQ2_S5_inferentialCurve_statistics.m
% 
% Cristina Gil Avila, 08.11.2023

% define and create output folder
datapath = '/rechenmagd3/Experiments/2023_1overf/results/sca';
figures_folder = '../results/figures/sca_inference/';
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
% Predominant direction of effect for the original data
pos = sum(d_orig > 0);
neg = sum(d_orig < 0);
if pos > neg
    sign_orig = 1;
elseif neg > pos
    sign_orig = -1;
end

% LOAD RANDOMIZED SPECIFICTIONS 
d_rand = nan(nSpec,nRand);
bf_rand = nan(nSpec,nRand);
sign = nan(1,nRand);
for iRand=1:nRand
    % load statistics of randomizations
    results = readtable(fullfile(datapath,sprintf('specs_bf_rand%.3d.txt',iRand)));
    [d_rand(:,iRand),ix] = sort(results.effect_size);
    bf_temp = results.bayes_factor;
    bf_rand(:,iRand) = bf_temp(ix);
    
    % Calculate the dominant sign of effects for each randomization
    pos = sum(d_rand(:,iRand) > 0);
    neg = sum(d_rand(:,iRand) < 0);
    if pos > neg
        sign(iRand) = 1;
    elseif neg > pos
        sign(iRand) = -1;
    end
end
% Flip the curves that do not have the dominant sign. The dominant sign is
% based on the original curve.
d_2s = d_rand.*sign*sign_orig;

%% INFERENTIAL TESTS

% MEDIAN TEST
% ===========================================================================
% Cristina's version (not flipping). Not sure if this approach is correct.
median_d_rand = sort(median(d_rand)); % sorting here is not needed for the test, just for plotting tcrit
tcrit_low = sum(median_d_rand <= prctile(median_d_rand,2.5));
tcrit_high = sum(median_d_rand <= prctile(median_d_rand,97.5));
t = min(sum(median_d_rand >= median(d_orig)),sum(median_d_rand <= median(d_orig)));
p_median_cris = 2*(t/nRand);
% Other optionis to do a z-test. Standarize the distribution of median values (z-score)
[z,mu,sigma] = zscore(median(d_rand));
z_orig = (median(d_orig)-mu)/sigma;
if z_orig > 0
    p_z = 2*sum(z>z_orig)/nRand;
else
    p_z = 2*sum(z<z_orig)/nRand;
end

% Original version (2 sided)
median_d_rand_2s = sort(median(d_2s)); % For the test data is not sorted (plot 2), but the papers' plot is sorted (plot 6)
% You have to decide whether to test if the median value is greater as
% observed or smaller as observed based on the orig_sign.
if sign_orig == 1
    % Number of simulations with median value at least as great as observed
    t = sum(median_d_rand_2s>=median(d_orig));
elseif sign_orig == -1 
    % Number of simulations with median value at least as smaller as observed
    t = sum(median_d_rand_2s<=median(d_orig));
end
p_median = t/nRand;

% One sided
if sign_orig == 1
    t = sum(median_d_rand >= median(d_orig));
elseif sign_orig == -1 
    t = sum(median_d_rand <= median(d_orig));
end
p_median_1s = t/nRand;

% Last option, for completeness (Same result as original version because when we compute the median the data has to be sorted)
median_d_rand_2s_sorted = sort(median(sort(d_2s)));
t = sum(median_d_rand_2s_sorted>=median(d_orig));
p_median_sorted = t/nRand;


% SHARE OF SIGNIFICANT RESULTS TEST
% ===========================================================================
% Extract the effects (d) of the 'significant' specifications
significant_d_orig = d_orig(bf_orig>3);
% Count only the significant specifications that have the same sign as the
% dominant curve effect
if ~isempty(significant_d_orig)
    BF_orig_count = sum((significant_d_orig * sign_orig)>0);
else
    BF_orig_count = 0;
end
% Extract the effects of 'significant' specifications for all
% randomizations
significant_d_rand = nan(nSpec,nRand); 
significant_d_rand(bf_rand>3) = d_2s(bf_rand>3);
% Count only the significant specifications that have the same sign as the
% dominant curve effect
BF_rand_count = sum((significant_d_rand .* sign) >0,1);

% One sided test. Check only if the share of specifications is higher than would be expected if all specifications had an effect of zero.
BF_rand_count = sort(BF_rand_count);
tcrit_high = sum(BF_rand_count <= prctile(BF_rand_count,95));
t = sum(BF_orig_count>BF_rand_count); 
p_share = 1-t/nRand;
test_pos_share = t>tcrit_high;

% Plot the share of significant results test
figure;
histogram(BF_rand_count)
hold on
x1 = xline(BF_orig_count,'k','LineWidth',2);
hold on
x2 = xline(prctile(BF_rand_count,5),'r');
xlabel('Number of significant specifications with an effect in the dominant direction')
ylabel('Randomizations')
legend([x1, x2],'observed number of significant specifications','critical statistic 5%')
title('Share of significant specifications test');



% AGGREGATE ALL BF TEST
% ==========================================================================
% We have to log(BF) before averaging to account for very small BF
% indicating strong effect in favor of the null
avgBF_orig = mean(log(bf_orig));
avgBF_rand = sort(mean(log(bf_rand),1));
tcrit_high = sum(avgBF_rand <= prctile(avgBF_rand,95));
% Indication of a positive effect (observed BF higher than the distribution of
% random BF)
t = sum(avgBF_orig>avgBF_rand);
% Indication of a null effect (obseved BF lower than the distribution of
% random BF)
% t = sum(avgBF_orig<avgBF_rand);
% How to make this test two-sided??? 
p_aggregateBF = 1-t/nRand;

figure;
histogram(avgBF_rand)
hold on
x1 = xline(avgBF_orig,'k','LineWidth',2);
hold on
x2 = xline(prctile(avgBF_rand,95),'r');
xlabel('Average logBF across specifications')
ylabel('Randomizations')
legend([x1, x2],'observed avg logBF across specifications','critical statistic 5%')
title('Aggregated BF test')



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