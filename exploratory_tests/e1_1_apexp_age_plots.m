% Plot relationship between aperiodic exponents and age
%
% Cristina Gil, TUM, 04.03.2024

% Load data
datapath = '../../results/statistics/exp_PFC_real.csv';
data = readtable(datapath);

% Fit a linear model to get the scatter plot nice
model = fitlm(data.age, data.exp_PFC);

f = figure('Units','centimeters','Position',[25 25 11 9]);
h = plot(model);
h(1).Marker = 'o';
h(1).MarkerFaceColor = '#94B74E';
h(1).MarkerEdgeColor = 'none';
h(2).Color = [0.5 0.5 0.5];
h(2).LineWidth = 1.5;
h(3).Color = [0.5 0.5 0.5];
h(3).LineWidth = 1.5;
h(4).Color = [0.5 0.5 0.5];
h(4).LineWidth = 1.5;
xlabel('Age');
ylabel('Aperiodic exponents');
title('');
xlim([18 86]);
box off;
legend off;
grid on;
% Save figure
saveas(f,'../../results/figures/e1_age.svg');