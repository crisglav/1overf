% Inferential statistics for the specifiction curve analysis
%
% Based on Elisabeth's May script
% /rechenmagd4/Experiments/2020_10_alpha_peak_frequency/code/matlab/functions/APF_RQ2_S5_inferentialCurve_statistics.m
% 
% Cristina Gil Avila, 08.11.2023, TUM

clear all, close all,

% define and create output folder
sca_path = '/rechenmagd3/Experiments/2023_1overf/results/sca/';
randomizations_path = '/rechenmagd3/Experiments/2023_1overf/results/sca/randomizations_e1_pa';
figures_path = '../../results/figures/sca_inference/';
addpath('/rechenmagd3/Experiments/2023_1overf/code/')
if ~exist(figures_path,'dir')
    mkdir(figures_path);
end

% Direction of hypothesis has to be specified by the researcher
tail = 'both'; % 'right', 'left'

%% LOAD THE DATA
nSpec = 48;
nRand = 500;

% LOAD ORIGINAL CURVE
results_orig = readtable(fullfile(randomizations_path,'stats_orig.csv'));
[d_orig, ix] = sort(results_orig.R);
bf_temp = results_orig.BF;
bf_orig = bf_temp(ix); % Sort according the effect size
pvalues_temp = results_orig.pvalue;
pvalues_orig = pvalues_temp(ix);

% LOAD RANDOMIZED SPECIFICTIONS 
d_rand = nan(nSpec,nRand);
bf_rand = nan(nSpec,nRand);
pvalues_rand = nan(nSpec,nRand);
for iRand=1:nRand
    % load statistics of randomizations
    results = readtable(fullfile(randomizations_path,sprintf('stats_rand%.3d.csv',iRand)));
    [d_rand(:,iRand),ix] = sort(results.R);
    bf_temp = results.BF;
    bf_rand(:,iRand) = bf_temp(ix);
    pvalues_temp = results.pvalue;
    pvalues_rand(:,iRand) = pvalues_temp(ix);   
end

%% INFERENTIAL TESTS
[p_median,p_share,p_aggregate] = sca_inference_tests(d_orig,bf_orig,pvalues_orig,d_rand,bf_rand,pvalues_rand,tail);


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
f = figure;
plot(1:nSpec,d_orig, 'Color', '#5400a2', 'LineWidth', 1.5); hold on;
plot(1:nSpec,sort(prctile(d_rand,50,2)),'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(lowPC), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
plot(1:nSpec,sort(highPC), '--', 'Color', grey, 'LineWidth', 1.5), hold on;
xlim([1 nSpec]);
ylim([-0.5,0.5])
yline(0);
title('Null distribution of specification curves');
ylabel('Effect size (R)')
xlabel('Specifications (nr, sorted by effect size)')
set(gca, 'TickDir', 'out');
box off
legend({'Observed data','Median under the null','2.5th and 97.5th under the null'})
text(max(xlim)-2, -0.4, sprintf('p_{median} = %.2f      p_{BFcount} = %.2f      p_{aggr} = %.2f', p_median, p_share,p_aggregate), 'Horiz', 'right', 'Vert', 'top', 'LineWidth', 1);
saveas(f,fullfile(figures_path,'e1_sca_inference_pa.svg'));

% Save data in csv
t = table(d_orig,sort(prctile(d_rand,50,2)),lowPC,highPC,'VariableNames',{'d_orig','d_median','lowPC','highPC'});
writetable(t, fullfile(sca_path,'e1_sca_inference_pa.csv'));

% % Histogram of the median values
% figure;
% histogram(median_d_rand,'EdgeColor','none','FaceColor',grey,'Orientation','horizontal')
% ylim([-0.5,0.5])
% ylabel('Effect size (Cohens d)')
% xlabel('Count randomizations')
% yline(median_d_orig,'color',[0 0.4470 0.7410],'LineWidth', 1.5)
% title('Histogram of median values')
% set(gca, 'TickDir', 'out');
% box off

% % Significant specifications
% figure;
% plot(1:nSpec,sort(d_rand),'Color', grey, 'LineWidth', 0.5), hold on;
% plot(1:nSpec,d_orig, 'Color', '#5400a2', 'LineWidth', 1.5); hold on;
% scatter(find(bf_orig > 3),significant_d_orig,30,[0 0.4470 0.7410],'filled');
% for i=1:nRand
%     if ~isempty(find(bf_rand(:,i) > 3, 1))
%         scatter(1:nSpec,significant_d_rand(:,i),5,'k','filled'); hold on;
%     end
% end
% xlim([1 nSpec]);
% ylim([-0.5,0.5])
% yline(0);
% ylabel('Effect size (Cohens d)')
% xlabel('Specifications (nr, sorted by effect size)')
% title('Significant specifications')
% box off
% set(gca, 'TickDir', 'out');
