
% select electrode
iSession = 1;
electrode = 'global';

% load original r values and plot them
results = readtable(['/rechenmagd4/Experiments/2020_10_alpha_peak_frequency/results/stats/results_RQ2_session' num2str(iSession) '_all_specs_orig_stats.csv']);
if strcmp(electrode, 'global')
    rowsElec = find(results.electrodes == "global");
elseif strcmp(electrode, 'S1')
    rowsElec = find(results.electrodes ~= "global");
end
[rho_orig, ix]= sort(results.R(rowsElec,:));
medianR_orig = median(rho_orig);
Bayes10_orig = length(find(results.BF10(rowsElec,:) > 3));

% create some variables for collection of results
rhos_rands = [];
medianR_rands = [];
Bayes10_rands = [];
f = figure;

for iRand = 1:500

    % load statistics of randomizations
    results = readtable(['/rechenmagd4/Experiments/2020_10_alpha_peak_frequency/results/stats/results_RQ2_session' num2str(iSession) '_all_specs_rand' num2str(iRand) '_stats.csv']);

    % select data from current electrode of interest
    if strcmp(electrode, "global")
        rowsElec = find(results.electrodes == "global");
    elseif strcmp(electrode, "S1")
        rowsElec = find(results.electrodes ~= "global");
    end
    rhos_iRand = sort(results.R(rowsElec,:));
    Bayes10 = results.BF10(rowsElec,:);

    %             % ONLY DO THIS FOR TWO-SIDED TESTS!
    %             % check the dominant sign of effects - if majority is positive, reverse
    %             % signs of all estimates and then sort (see Simonsohn et al. 2020)
    %             nr_negativeR = length(find(rhos_iRand < 0));
    %             if nr_negativeR <= length(rhos_iRand)/2
    %                 rhos_iRand = sort(rhos_iRand*-1);
    %             end
    
    % collect all single correlation values, medians and count of Bayes10
    % factors larger than 3 for every randomization/shuffled data set
    rhos_rands = [rhos_rands, rhos_iRand];
    medianR_rands = [medianR_rands, median(rhos_iRand)];
    count_Bayes_H1 = length(find(Bayes10 > 3));
    Bayes10_rands = [Bayes10_rands, count_Bayes_H1];
    %                     plot(1:length(rhos_iRand),rhos_iRand, 'Color', [.7 .7 .7], 'LineWidth', 0.1); hold on;

end

% calculate 5, 50 and 95 percentiles of random spec curves
% here we dont do 2.5 and 97.5 CI since we are doing a one-sided
% test and want to check if the observed median is lower than the
% simulated ones with a 5% probability. strictly speaking, the uper
% interval does not really makes sense but would be weird not to
% plot
medianCurves = prctile(rhos_rands,50,2)';
lowPC = prctile(rhos_rands,5,2)';
highPC = prctile(rhos_rands,95,2)';

% plot
grey = [0.8 0.8 0.8];
plot(1:length(medianCurves),medianCurves,'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:length(medianCurves),lowPC, '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:length(medianCurves),highPC, '--', 'Color', grey, 'LineWidth', 1.5), hold on;

plot(1:length(rho_orig),rho_orig, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); hold on;
ylim([-0.3 0.3]);
test = yline(0);
ylabel("Pearson's R");
xlabel("specifications (nr, sorted by effect size)");
set(gca, 'TickDir', 'out');
box off
if strcmp(electrode, "global")
    xlim([1 72]);
elseif strcmp(electrode, "S1")
    xlim([1 144]);
end
yticks(-0.3:0.1:0.3);


% compute p values and add to title
% percentage of randomization with median effect size smaller than original
% negative correlation
p_median = 1-length(find(medianR_orig < medianR_rands))/length(medianR_rands);
% percentage of randomizations with more Bayes10 values > 3 than original
p_BF10 = 1-length(find(Bayes10_orig > Bayes10_rands))/length(medianR_rands);
%         title(electrode{1,1});
text(max(xlim)-2, -0.25, ['p median = ' num2str(round(p_median,2)) ', p BF10 = ' num2str(round(p_BF10, 2))], 'Horiz', 'right', 'Vert', 'top', 'LineWidth', 1);

orient(f,'landscape');
%         print('-fillpage', [output_folder 'RQ2_session' num2str(iSession) '_' electrode{1,1} '_infSpecCurve.pdf'], '-dpdf');
print([output_folder 'RQ2_session' num2str(iSession) '_' electrode{1,1} '_infSpecCurve.png'], '-dpng');
close;

%% Cristina
nRand = 500;
nSpec = 72;

% LOAD DATA
datapath = '/rechenmagd4/Experiments/2020_10_alpha_peak_frequency/results/stats/';
% Original results
results_orig = readtable(fullfile(datapath,'results_RQ2_session1_all_specs_orig_stats.csv'));
rowsElec = find(results_orig.electrodes == "global");
R_orig = sort(results_orig.R(rowsElec,:));
BF_orig = results_orig.BF10(rowsElec,:);
pos = sum(R_orig > 0);
neg = sum(R_orig < 0);
if pos > neg
    sign_orig = 1;
elseif neg > pos
    sign_orig = -1;
end
% Randomizations
R_rand = nan(nSpec,nRand);
BF_rand = nan(nSpec,nRand);
sign = nan(1,nRand);
for i=1:nRand
    % load statistics of randomizations
    results = readtable(fullfile(datapath,['results_RQ2_session1_all_specs_rand' num2str(i) '_stats.csv']));
    rowsElec = find(results.electrodes == "global");
    R_rand(:,i) = sort(results.R(rowsElec,:));
    BF_rand(:,i) = results.BF10(rowsElec,:);
    
    % Calculate the dominant sign of effects for each randomization
    pos = sum(R_rand(:,i) > 0)/nSpec;
%     neg = sum(R_rand(:,i) < 0);
    if pos > 0.5
        sign(i) = 1;
    elseif pos <= 0.5
        sign(i) = -1;
    end
end
% The original effect was negative, and we expect to see a negative effect.
% Reverse the sign of those randomizations where the majority of effects
% are positive, i.e. reverse the sign of all the randomizations that do not
% have the dominant sign of the original specification.
R_2s = R_rand.*sign*sign_orig;

%% Plot curves
figure;
grey = [0.8 0.8 0.8];

tiledlayout(2,3)

nexttile
plot(1:nSpec,R_rand), hold on;
plot(1:nSpec,R_orig, 'Color', [0 0.4470 0.7410], 'LineWidth', 2); hold on;
scatter(36,R_orig(36),'black','filled'); hold on;
yline(0,'LineWidth', 1.5);
ylim([-0.5,0.5])
title('All randomizations')

nexttile
plot(1:nSpec,R_2s), hold on;
plot(1:nSpec,R_orig, 'Color', [0 0.4470 0.7410], 'LineWidth', 2); hold on;
scatter(36,R_orig(36),'black','filled'); hold on;
yline(0,'LineWidth', 1.5);
ylim([-0.5,0.5])
title('All randomizations, flipped sign')

nexttile
plot(1:nSpec,sort(R_2s)), hold on;
plot(1:nSpec,R_orig, 'Color', [0 0.4470 0.7410], 'LineWidth', 2); hold on;
scatter(36,R_orig(36),'black','filled'); hold on;
yline(0,'LineWidth', 1.5);
ylim([-0.5,0.5])
title('All randomizations, flipped sign, sorted')

% Original curve with median and prct
nexttile
plot(1:nSpec,sort(prctile(R_rand,50,2)),'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(prctile(R_rand,2.5,2)), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(prctile(R_rand,97.5,2)), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,R_orig, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); hold on;
scatter(36,R_orig(36),'black','filled'); hold on;
ylim([-0.5,0.5])
yline(0);

nexttile
plot(1:nSpec,sort(prctile(R_2s,50,2)),'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(prctile(R_2s,2.5,2)), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(prctile(R_2s,97.5,2)), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,R_orig, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); hold on;
scatter(36,R_orig(36),'black','filled'); hold on;
ylim([-0.5,0.5])
yline(0);

nexttile
plot(1:nSpec,sort(prctile(sort(R_2s),50,2)),'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(prctile(sort(R_2s),2.5,2)), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(prctile(sort(R_2s),97.5,2)), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,R_orig, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); hold on;
scatter(36,R_orig(36),'black','filled'); hold on;
ylim([-0.5,0.5])
yline(0);


%%
% MEDIAN TEST
% Cristina's version
median_R_rand = sort(median(R_rand));
tcrit_low = sum(median_R_rand <= prctile(median_R_rand,2.5));
tcrit_high = sum(median_R_rand <= prctile(median_R_rand,97.5));
t = min(sum(median(R_orig)<median_R_rand),sum(median(R_orig)>median_R_rand));
p_median_cris = 2*(t/nRand);
test_pos = or(t<tcrit_low,t>tcrit_high);

% Original version
median_R_rand_2s = sort(median(R_2s)); % For the test data is not sorted (plot 2), but the papers' plot is sorted (plot 6)
% Proportion of simulations with median value at least as small as observed
t = sum(median_R_rand_2s<=median(R_orig));
p_median = t/nRand;

% Plot the median ttest Cristina's version
figure;
tiledlayout(1,3)
nexttile
histogram(median_R_rand);
hold on
x1 = xline(median(R_orig),'k','LineWidth',2);
hold on
x2 = xline(prctile(median_R_rand,2.5),'r');
hold on
xline(prctile(median_R_rand,97.5),'r')
xlabel('Median R')
ylabel('Randomizations')
box off;
text(min(xlim), max(ylim), ['p = ' num2str(round(p_median_cris,2))], 'Horiz', 'right', 'Vert', 'top', 'LineWidth', 1);
title('Median test');

% Plot the median ttest original version
nexttile
histogram(median_R_rand_2s,'BinWidth',0.02);
hold on
x1 = xline(median(R_orig),'k','LineWidth',2);
hold on
xline(prctile(median_R_rand_2s,5),'r')
box off;
text(min(xlim), -max(ylim), ['p = ' num2str(round(p_median,2))], 'Horiz', 'right', 'Vert', 'top', 'LineWidth', 1);
title('Median test, flipped sign');

% For completeness plot the third version with sorted values to check that
% it is something different
median_R_rand_2s_sorted = sort(median(sort(R_2s))); % For the test data is not sorted (plot 2), but the papers' plot is sorted (plot 6)
% Proportion of simulations with median value at least as small as observed
t = sum(median_R_rand_2s_sorted<=median(R_orig));
p_median_sorted = t/nRand;
nexttile
histogram(median_R_rand_2s_sorted,'BinWidth',0.02);
hold on
x1 = xline(median(R_orig),'k','LineWidth',2);
hold on
xline(prctile(median_R_rand_2s_sorted,5),'r')
xlabel('Median d')
ylabel('Randomizations')
legend([x1, x2],'observed d','critical statistic 5%')
box off;
text(min(xlim), max(ylim), ['p = ' num2str(round(p_median_sorted,2))], 'Horiz', 'right', 'Vert', 'top', 'LineWidth', 1);
title('Median test, flipped sign,sorted');





% SHARE OF SIGNIFICANT RESULTS TEST
% Extract the effects (R) of the 'significant' specifications
significant_R_orig = R_orig(BF_orig>3);
% Count only the significant specifications have the same sign as the
% dominant curve effect
if ~isempty(significant_R_orig)
    BF_orig_count = sum((significant_R_orig * sign_orig)>0);
else
    BF_orig_count = 0;
end
% Extract the effects of 'significant' specifications
significant_R_rand = nan(nSpec,nRand); 
significant_R_rand(BF_rand>3) = R_rand(BF_rand>3);
% Count only the significant specifications that have the same sign as the
% dominant curve effect
BF_rand_count = sum((significant_R_rand .* sign) >0,1);

% One sided test. Check only if the share of specifications is higher than would be expected if all specifications had an effect of zero.
BF_rand_count = sort(BF_rand_count);
tcrit_high = sum(BF_rand_count <= prctile(BF_rand_count,95));
t = sum(BF_orig_count>BF_rand_count); 
p_share = 1-t/nRand;
test_pos = t>tcrit_high;

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
% We have to log(BF) before averaging to account for very small BF
% indicating strong effect in favor of the null
avgBF_orig = mean(log(BF_orig));
avgBF_rand = sort(mean(log(BF_rand),1));
tcrit_high = sum(avgBF_rand <= prctile(avgBF_rand,95));
t = sum(avgBF_orig>avgBF_rand); 
p_aggregateBF = 1-t/nRand;
test_pos = t>tcrit_high;

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