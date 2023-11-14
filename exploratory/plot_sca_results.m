addpath('/rechenmagd4/toolboxes_and_functions/plotting_functions')

load('/rechenmagd3/Experiments/2023_1overf/results/sca/specs_ap_exp.mat');
c = lines(1);
figure;
tlc = tiledlayout(3,1);
nexttile
stdshade(exp',0.2,c)
xlabel('Specifications')
ylabel('Aperiodic exponent')

nexttile;
stdshade(r_squared',0.2,c)
xlabel('Specifications')
ylabel('Model R squared')

nexttile;
stdshade(mae',0.2,c)
xlabel('Specifications')
ylabel('MAE')

% figure;
% for i=1:264
%     scatter(1:264,exp(i,:));
%     hold on
% end