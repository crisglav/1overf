%% Example plot with two PSD with different slopes
% This is to show what would be a PSD with higher E/I ratio and lower E/I
% ration.

expPath = '../../results/features/electrode_space';
bidsID = 'sub-035_task-closed';
load(fullfile(expPath,[bidsID '_fm.mat']));

% Identify in the plot the two lines with flatter and steeper slope
figure;
LineList = plot(fm.freq_orig,log10(fm.pow_orig));
set(LineList, 'ButtonDownFcn',{@myLineCallback, LineList})

%%
figure;
tiledlayout(2,2);
nexttile
plot(fm.freq_orig,log10(fm.pow_orig(51,:)));
hold on
plot(fm.freqs,fm.ap_fit_group(51,:));
ylim([-3 1.5])
title(label(51));
ylabel('log Power');
xlabel('Frequency');
str = sprintf('ap exp = %0.2f', fm.group_results{51}.aperiodic_params(end));
text(40,1,str);
nexttile
plot(fm.freq_orig,log10(fm.pow_orig(61,:)));
hold on
plot(fm.freqs,fm.ap_fit_group(61,:));
ylim([-3 1.5])
title(label(61));
str = sprintf('ap exp = %0.2f', fm.group_results{61}.aperiodic_params(end));
text(40,1,str);
xlabel('Frequency');



nexttile
plot(log10(fm.freq_orig),log10(fm.pow_orig(51,:)));
hold on
plot(log10(fm.freqs),fm.ap_fit_group(51,:));
ylim([-3 1.5])
ylabel('log Power');
xlabel('log Frequency');
nexttile
plot(log10(fm.freq_orig),log10(fm.pow_orig(61,:)));
hold on
plot(log10(fm.freqs),fm.ap_fit_group(61,:));
ylim([-3 1.5])
xlabel('log Frequency');

function myLineCallback(LineH, EventData, LineList)
% disp(LineH);
% disp(get(LineH,'YData'));
disp(find(LineList == LineH));
set(LineList,'LineWidth',0.5);
set(LineH, 'LineWidth', 2.5);
uistack(LineH, 'top');
end