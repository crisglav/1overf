% Inferential statistics for the specifiction curve analysis
%
% Based on Elisabeth's May script
% /rechenmagd4/Experiments/2020_10_alpha_peak_frequency/code/matlab/functions/APF_RQ2_S5_inferentialCurve_statistics.m
% 
% Cristina Gil Avila, 08.11.2023, TUM

% define and create output folder
datapath = '/rechenmagd3/Experiments/2023_1overf/results/sca/Rinterface';
figures_folder = '../results_v3/figures/sca_inference/';
if ~exist(figures_folder,'dir')
    mkdir(figures_folder);
end

% Hardcoded number of randomizations and specifications
nRand = 500;
nSpec = 48;
% Hardcoded direction of hypothesis
tail = 'both'; % 'right', 'both'

%% LOAD ORIGINAL CURVE
results_orig = readtable(fullfile(datapath,'stats_orig.csv'));
[d_orig, ix] = sort(results_orig.d);
bf_temp = results_orig.BF;
bf_orig = bf_temp(ix); % Sort according the effect size
pvalues_temp = results_orig.pvalue;
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
    results = readtable(fullfile(datapath,sprintf('stats_rand%.3d.csv',iRand)));
    [d_rand(:,iRand),ix] = sort(results.d);
    bf_temp = results.BF;
    bf_rand(:,iRand) = bf_temp(ix);
    pvalues_temp = results.pvalue;
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


% %% LOAD ELISABETHS DATA
% nRand = 500;
% nSpec = 72;
% 
% % Load original curve
% datapath = '/rechenmagd4/Experiments/2020_10_alpha_peak_frequency/results/stats/';
% results_orig = readtable(fullfile(datapath,'results_RQ2_session1_all_specs_orig_stats.csv'));
% rowsElec = find(results_orig.electrodes == "global");
% [d_orig, ix]= sort(results_orig.R(rowsElec,:));
% bf_temp = results_orig.BF10(rowsElec,:);
% bf_orig = bf_temp(ix);
% pvalues_temp = results_orig.pvalue(rowsElec,:);
% pvalues_orig = pvalues_orig(ix);
% % Predominant direction of effect for the original data
% pos = sum(d_orig > 0);
% neg = sum(d_orig < 0);
% if pos >= neg
%     sign_orig = 1;
% elseif neg > pos
%     sign_orig = -1;
% end
% 
% % LOAD RANDOMIZED SPECIFICTIONS
% % create some variables for collection of results
% d_rand  =  nan(nSpec,nRand);
% bf_rand    = nan(nSpec,nRand);
% 
% for iRand = 1:nRand
% 
%     % load statistics of randomizations
%     results = readtable([datapath 'results_RQ2_session1_all_specs_rand' num2str(iRand) '_stats.csv']);
%     rowsElec = find(results.electrodes == "global");
%     % sort according to effect size
%     [d_rand(:,iRand), ix] = sort(results.R(rowsElec,:));
%     % again, sort Bayes Factors in same way to keep assignment
%     bf_temp = results.BF10(rowsElec,:);
%     bf_rand(:,iRand)   = bf_temp(ix);
% 
%     % calculate the dominant sign of effects for each randomization
%     pos = sum(d_rand(:,iRand) > 0);
%     if (pos/nSpec) > 0.5
%         sign(iRand) = 1;
%     elseif (pos/nSpec) < 0.5
%         sign(iRand) = -1;
%     else
%         sign(iRand) = sign_orig; % if there is no dominant sign of effects, assume direction of original effect
%     end
% end
% d_2s = d_rand.*sign*sign_orig;
%% INFERENTIAL TESTS

% 1. TEST OF MEDIAN EFFECT SIZE
% ===================================================================
median_d_orig = median(d_orig);
median_d_rand = median(d_rand);
switch tail
    case 'left'
        % One-sided test
        % ==============
        % calculate p-value as percentage of randomizations with median
        % effects size smaller than original effect size
        t = sum(median_d_rand <= median_d_orig);
    case 'right'
        % One-sided test
        % ==============
        % calculate p-value as percentage of randomizations with median
        % effects size larger than original effect size
        t = sum(median_d_rand >= median_d_orig);
    case 'both'
        % Two-sided test
        % ==============
        % calculate p-value as percentage of randomizations with median
        % effects size more extreme than original effect size on both ends
        % of the distribution of median effect sizes
        t = sum(abs(median_d_rand) >= abs(median_d_orig));
end
p_median = t/nRand;

% % Test flipping the curves (2 sided)
% % ==============
% % This test is conceptually the same as the two-sided test.
% median_d_rand_2s = median(d_2s);
% if sign_orig == 1
%     t = sum(median_d_rand_2s >= median_d_orig);
% elseif sign_orig == -1 
%     t = sum(median_d_rand_2s <= median_d_orig);
% end
% p_median_flip = t/nRand;


% 2. TEST OF SHARE OF 'SIGNIFICANT' RESULTS (HERE SIGNIFICANT: BF > 3)
% ====================================================================
% Extract the effects (d) of the 'significant' specifications for the
% original curve and the randomizations
significant_d_orig = d_orig(bf_orig > 3);
significant_d_rand = nan(nSpec,nRand); 
significant_d_rand(bf_rand > 3) = d_rand(bf_rand > 3);

switch tail
    case 'left'
        % One-sided test
        % ==============
        % Count the significant specifications that have a negative sign
        if ~isempty(significant_d_orig)
            BF_orig_count = sum(significant_d_orig < 0);
        else
            BF_orig_count = 0;
        end
        BF_rand_count = sum((significant_d_rand < 0),1);
        % Percentage of randomizations that have a number of negative significant
        % specifications higher than the original numer of negative signficiant
        % specifications
        t = sum(BF_rand_count >= BF_orig_count);

    case 'right'
        % One-sided test
        % ==============
        % Count the significant specifications that have a positive sign
        if ~isempty(significant_d_orig)
            BF_orig_count = sum(significant_d_orig > 0);
        else
            BF_orig_count = 0;
        end
        BF_rand_count = sum((significant_d_rand > 0),1);
        % Percentage of randomizations that have a number of positive significant
        % specifications higher than the original numer of positive signficiant
        % specifications
        t = sum(BF_rand_count >= BF_orig_count);

    case 'both'
        % Two-sided test
        % ==============
        % Count only the significant specifications that have the same sign as the
        % dominant curve effect
        if ~isempty(significant_d_orig)
            BF_orig_count = sum((significant_d_orig * sign_orig) > 0);
        else
            BF_orig_count = 0;
        end
        % For each randomization, count the significant specifications
        % that have the same sign as the dominant curve effect of that
        % randomization (!)
        BF_rand_count = sum((significant_d_rand .* sign) > 0,1);

% 
%         % Count the significant specifications that have either a positive or negative sign
%         if ~isempty(significant_d_orig)
%             BF_orig_count = sum(abs(significant_d_orig) > 0);
%         else
%             BF_orig_count = 0;
%         end
%         BF_rand_count = sum(abs(significant_d_rand) > 0 ,1);


        % Percentage of randomizations that have a number of significant specifications higher
        % than the original number of signficiant specifications (either
        % positive or negative)
        t = sum(BF_rand_count >= BF_orig_count);
end
p_share = t/nRand;

% % OLD VERSION
% % Count only the significant specifications that have the same sign as the
% % dominant curve effect
% if ~isempty(significant_d_orig)
%     BF_orig_count = sum((significant_d_orig * sign_orig) > 0);
% else
%     BF_orig_count = 0;
% end
% 
% % One-sided test
% % ==============
% % Count the significant specifications that have the same sign
% % as the dominant curve effect from the original curve (!)
% BF_rand_count = sum((significant_d_rand * sign_orig) > 0,1);
% t = sum(BF_rand_count>=BF_orig_count); 
% p_share_1s = t/nRand;
% 
% % Two-sided test
% % ==============
% % For each randomization, count the significant specifications
% % that have the same sign as the dominant curve effect of that
% % randomization (!)
% BF_rand_count = sum((significant_d_rand .* sign) > 0,1);
% t = sum(BF_rand_count>=BF_orig_count); 
% p_share_2s = t/nRand;
% 
% % Test flipping the curves (2 sided)
% % ==============
% % This test is conceptually the same as the two-sided test.
% % Calculate the dominant sign based on the flipped curves (edom)
% edom = nan(nSpec,nRand);
% if sign_orig == 1
%     edom(d_2s >= 0) = 1;
%     edom(d_2s < 0) = 0;
% else
%     edom(d_2s >= 0) = 0;
%     edom(d_2s < 0) = 1;
% end
% sig_freq = bf_rand > 3;
% sig_dom = and(sig_freq,edom);
% sig_dom_freq = sum(sig_dom);
% t = sum(sig_dom_freq >= BF_orig_count);
% p_share_flip = t/nRand;


% 3. TEST OF AGGREGATED P-VALUES
% =================================================================
% For each randomization, average pvalues following the Stouffer's method.
% First get z-scores for each pvalue and average across specifications
% Note that all versions of the test are conceptually the same.

switch tail
    case 'left'
        % One-sided test
        % ==============
        % calculate p-value as percentage of randomizations with average z value
        % smaller than original average z value
        z_temp = norminv(pvalues_orig);
        z_orig = sum(z_temp)/sqrt(nSpec);

        z_temp = norminv(pvalues_rand);
        z_rand = sum(z_temp,1)./sqrt(nSpec);

        t = sum(z_rand <= z_orig);
    case 'right'
        % One-sided test
        % ==============
        % calculate p-value as percentage of randomizations with average z value
        % larger than original average z value
        z_temp = norminv(1-pvalues_orig);
        z_orig = sum(z_temp)/sqrt(nSpec);

        z_temp = norminv(1-pvalues_rand);
        z_rand = sum(z_temp,1)./sqrt(nSpec);

        t = sum(z_rand >= z_orig);
    case 'both'
        % Two-sided test
        % ==============
        % calculate p-value as percentage of randomizations with average z value
        % more extreme than original average z value on both ends of the distribution of z
        % values rand. Pvalues are diveded by two because the original t-tests were 
        % two-sided.  We do 1-pvalue to get positive zscores.
        z_temp = norminv(1-(pvalues_orig./2));
        z_orig = sum(z_temp)/sqrt(nSpec);

        z_temp = norminv(1-(pvalues_rand./2));
        z_rand = sum(z_temp,1)./sqrt(nSpec);

        t = sum(z_rand >= z_orig);
end
p_aggregate = t/nRand;


% % Test flipping the curves (2 sided)
% % ==============
% pvalues_2s = pvalues_rand.*sign;
% z_temp = norminv(pvalues_2s);
% z_rand_2s = sum(z_temp,1)./sqrt(nSpec);
% 
% if sign_orig == 1
%     t = sum(z_rand_2s >= z_orig);
% elseif sign_orig == -1 
%     t = sum(z_rand_2s <= z_orig);
% end
% p_aggregate_flip = t/nRand;

%% PLOT MEDIAN CURVES
% Calculate median curves and upper and lower percentile curves
% (i.e., confidence intervals)
% For a two-sided test, do 2.5 and 97.5 for CI
% For a one-sided test, do 5 and 95 because we want to check if the
% observed median is lower than the randomized ones with a 5%
% probability. Strictly speaking, plotting the upper/lower interval for a (left/right test) does
% not really make sense but not plotting would be weird.

% NOTE: curves of all randomizations are sorted (important for
% plotting)
switch tail
    case {'left','right'}
        lowPC = prctile(d_rand,5,2);
        highPC = prctile(d_rand,95,2);
    case 'both'
        lowPC = prctile(d_rand,2.5,2);
        highPC = prctile(d_rand,97.5,2);
end

grey = [0.8 0.8 0.8];
figure;
plot(1:nSpec,d_orig, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); hold on;
plot(1:nSpec,sort(prctile(d_rand,50,2)),'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(lowPC), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(highPC), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
xlim([1 nSpec]);
ylim([-0.5,0.5])
yline(0);
title('Null distribution of specification curves');
ylabel('Effect size (Cohens d)')
xlabel('Specifications (nr, sorted by effect size)')
set(gca, 'TickDir', 'out');
box off
legend({'Observed data','Median under the null','2.5th and 97.5th under the null'})
% text(max(xlim)-2, -0.25, ['p median/BF count/aggr p = ' num2str(round(p_median,3)) '/' num2str(round(p_share, 3)) '/' num2str(round(p_aggregate,3))], 'Horiz', 'right', 'Vert', 'top', 'LineWidth', 1);

% Histogram of the median values
figure;
histogram(median_d_rand,'EdgeColor','none','FaceColor',grey,'Orientation','horizontal')
ylim([-0.5,0.5])
ylabel('Effect size (Cohens d)')
xlabel('Count randomizations')
yline(median_d_orig,'color',[0 0.4470 0.7410],'LineWidth', 1.5)
title('Histogram of median values')
set(gca, 'TickDir', 'out');
box off

% Significant specifications
figure;
plot(1:nSpec,sort(d_rand),'Color', grey, 'LineWidth', 0.5), hold on;
plot(1:nSpec,d_orig, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); hold on;
scatter(find(bf_orig > 3),significant_d_orig,30,[0 0.4470 0.7410],'filled');
for i=1:nRand
    if ~isempty(find(bf_rand(:,i) > 3, 1))
        scatter(1:nSpec,significant_d_rand(:,i),5,'k','filled'); hold on;
    end
end
xlim([1 nSpec]);
ylim([-0.5,0.5])
yline(0);
ylabel('Effect size (Cohens d)')
xlabel('Specifications (nr, sorted by effect size)')
title('Significant specifications')
box off
set(gca, 'TickDir', 'out');
