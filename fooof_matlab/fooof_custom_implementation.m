% Script testing the reimplementation of fooof in Matlab
%
% Original python code by Tom Donoghue et. al in https://github.com/fooof-tools/fooof
%
% Cristina Gil Avila, TUM, 01.09.2023
close all
clear all

% FOOOF settings
self.freq_range =        [3 40]; % [3 40]
self.peak_width_limits = [1 8]; % [1 8] % [0.5 12];
self.max_n_peaks       = 6; % 3;     % inf;
self.min_peak_height   = 0.15; % 0.15;  % 0;
self.peak_threshold    = 2;
self.aperiodic_mode    = 'fixed'; %'knee'

% Fixed settings
self.ap_percentile_thresh = 0.025;


%% Loading demo data
power_path = '../../results/power/PFC/';
power_files = dir(fullfile(power_path,'*.mat'));
file = power_files(1);
load(fullfile(file.folder,file.name));

% Average across channels
avgpow = mean(pow,1);

%%
% pow_orig = avgpow;
% freq_orig = freq;
% 
% f = and(freq >= self.freq_range(1),freq<=self.freq_range(2));
% freq = freq_orig(f);
% pow = pow_orig(f);
% pow = log10(pow);
% 
% % Add data to the configuration settings
% self.freq = freq;
% self.pow = pow;

%% Fit fooof
self = fit_fooof(self);

%% Playzone
% I think it would be interesting to know how well the model fitted in the freq_range extrapolates to the full spectrum

% Add the original frequency to the model
self.freq_orig = freq_orig;
self.pow_orig = pow_orig;

% Calculate peak fit on the full frequency range
gaussian_params = reshape(self.gaussian_params',1,[]);
self.peak_fit_full = gaussian_function(gaussian_params,self.freq_orig);

% Run the aperiodic fit on the whole power spectrum
self.ap_fit_full = aperiodic_function(self.aperiodic_params(1), self.aperiodic_params(2), self.freq_orig);
self.spectrum_flat_full = self.pow_orig - self.ap_fit_full;

% Create full power spectrum model
self.fooofed_spectrum_full = self.peak_fit_full + self.ap_fit_full;

% ===== Calculate R^2 and error of the model fit ======
r2 = corr(log10(self.pow_orig)', self.fooofed_spectrum_full')^2;
mae = mean(abs(log10(self.pow_orig) - self.fooofed_spectrum_full));
mse = mean((log10(self.pow_orig) - self.fooofed_spectrum_full).^2);
rmse = sqrt(mean((log10(self.pow_orig) - self.fooofed_spectrum_full).^2));
%% Generate report

% Save fooof model


%% More plotting
% Overlap the fooofed spectrum to the original spectrum in full frequency range
padded_ap_fit = nan(size(freq_orig));
padded_ap_fit(f) = self.ap_fit;
padded_fooofed_spectrum = nan(size(freq_orig));
padded_fooofed_spectrum(f) = self.fooofed_spectrum;

fig = figure('Units','centimeters','Position',[15 15 35 16]);
tlc = tiledlayout(1,2);
nexttile;
plot(freq_orig,log10(pow_orig),'k');
hold on
plot(freq_orig,padded_ap_fit,'--b')
hold on
plot(freq_orig,padded_fooofed_spectrum,'r');
legend ({'Original spectrum','Full model fit','Aperiodic fit'});
xlabel('Frequency');
ylabel('logPower');
title('Model fit 3 - 40 Hz');
str1 = sprintf('r2 = %0.3f ',self.r2);
str2 = sprintf('mae = %0.3f',self.error_mae);
annotation(fig,'textbox',[0.15 0.05 .3 .3],'String',{str1,str2},'FitBoxToText','on');


nexttile
plot(freq_orig,log10(pow_orig),'k');
hold on
plot(freq_orig,self.ap_fit_full,'--b')
hold on
plot(freq_orig,self.fooofed_spectrum_full,'r');
legend ({'Original spectrum','Full model fit','Aperiodic fit'});
xlabel('Frequency');
ylabel('logPower');
title('Model fit 1 - 100 Hz');
str1 = sprintf('r2 = %0.3f ',r2);
str2 = sprintf('mae = %0.3f',mae);
annotation(fig,'textbox',[0.55 0.05 .3 .3],'String',{str1,str2},'FitBoxToText','on');

%% Some plotting
% Reproducing figures of tutorial 3 of fooof (jupyter notebook)
fig2 = figure('Units','centimeters','Position',[15 15 35 35]);
tlc = tiledlayout(3,3);

% Figure1, non-log scale, full frequency range
nexttile;
plot(freq_orig,log10(pow_orig),'k');
xlabel('Frequency');
ylabel('Power');
title('Original PSD');

% Figure 2: Initial robust aperiodic fit
nexttile;
plot(self.freq,self.pow,'k');
hold on
plot(self.freq,self.initial_ap_fit,'--b');
legend({'Original power spectrum','Initial robust aperiodic fit'});
xlabel('Freqency');
ylabel('logPower');
title('Step 1: Initial robust aperiodic fit');

% Figure 3: Flat spectrum
nexttile;
plot(self.freq,self.inital_spectrum_flat,'k');
legend ('Flattened Spectrum');
xlabel('Frequency');
ylabel('logPower');
title('Step 2: Removed aperiodic component from original PSD');

% Figure 5: Gaussian periodic fit
nexttile;
plot(self.freq,self.peak_fit,'g');
legend ('Final periodic fit');
xlabel('Frequency');
ylabel('logPower');
title('Step 3: Result of multi gaussian fit');

% Figure 6: Peak removed spectrum
nexttile;
plot(self.freq,self.spectrum_peak_rm,'k');
legend ('Peak removed spectrum');
xlabel('Frequency');
ylabel('logPower');
title('Step 4: Remove gaussians from original PSD');

% Figure 7: Final aperiodic fit
nexttile;
plot(self.freq,self.spectrum_peak_rm,'k');
hold on
plot(self.freq,self.ap_fit,'--b')
legend ({'Peak removed spectrum','Final Aperiodic Fit'});
xlabel('Frequency');
ylabel('logPower');
title('Step 5: Re-fit aperiodic component');

% Figure 8: Full model fit
nexttile;
plot(self.freq,self.fooofed_spectrum,'r');
legend ('Full model');
xlabel('Frequency');
ylabel('logPower');
title('Combined model: Gaussian + aperiodic fit');

% Figure 9
nexttile;
plot(self.freq,self.pow,'k');
hold on
plot(self.freq,self.ap_fit,'--b')
hold on
plot(self.freq,self.fooofed_spectrum,'r');
legend ({'Original spectrum','Full model fit','Aperiodic fit'});
xlabel('Frequency');
ylabel('logPower');
title('Assess goodness of fit');

%% OOP
fm = fooof();
fm = add_data(fm,freq,avgpow,[2 40]);
fm = fit(fm);

% Plot
figure
plot(fm.freqs,fm.power_spectrum,'k');
hold on
plot(fm.freqs,fm.ap_fit,'--b')
hold on
plot(fm.freqs,fm.fooofed_spectrum,'r');
legend ({'Original spectrum','Aperiodic fit','Full model fit'});
xlabel('Frequency');
ylabel('logPower');
title('Assess goodness of fit');