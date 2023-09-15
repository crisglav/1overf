% Script testing the reimplementation of fooof in Matlab
%
% Original python code by Tom Donoghue et. al in https://github.com/fooof-tools/fooof
%
% Cristina Gil Avila, TUM, 01.09.2023
close all
clear all

% FOOOF settings
self.freq_range        = [2 40];
self.peak_width_limits = [0.5 12]; % [1 8] % [0.5 12];
self.max_n_peaks       = inf; % 3;     % inf;
self.min_peak_height   = 0; % 0.15;  % 0;
self.peak_threshold    = 2;
self.aperiodic_mode    = 'knee';

% Fixed settings
self.ap_percentile_thresh = 0.025;


%% Load subset data and fooof it
power_path = '../../results/power/PFC/';
power_files = dir(fullfile(power_path,'*.mat'));
subset = power_files(7:12);
fooof = cell(1,length(subset));

for iFile =1:length(subset)
    file = subset(iFile);
    load(fullfile(file.folder,file.name));

    % Average across channels
    avgpow = mean(pow,1);
    pow_orig = avgpow;
    freq_orig = freq;

    f = and(freq >= self.freq_range(1),freq<=self.freq_range(2));
    freq = freq_orig(f);
    pow = pow_orig(f);
    pow = log10(pow);

    % Add data to the configuration settings
    self.freq = freq;
    self.pow = pow;

    % Fit fooof
    self = fit_fooof(self);
    
%     % How does the model fitted in the freq_range extrapolates to the full spectrum?
%     % Add the original frequency to the model
%     self.freq_orig = freq_orig;
%     self.pow_orig = pow_orig;
% 
%     % Calculate peak fit on the full frequency range
%     gaussian_params = reshape(self.gaussian_params',1,[]);
%     self.peak_fit_full = gaussian_function(gaussian_params,self.freq_orig);
% 
%     % Run the aperiodic fit on the whole power spectrum
%     self.ap_fit_full = exp_function(self.aperiodic_params, self.freq_orig);
%     self.spectrum_flat_full = self.pow_orig - self.ap_fit_full;
% 
%     % Create full power spectrum model
%     self.fooofed_spectrum_full = self.peak_fit_full + self.ap_fit_full;
% 
%     % ===== Calculate R^2 and error of the model fit ======
%     self.r2_full = corr(log10(self.pow_orig)', self.fooofed_spectrum_full')^2;
%     self.mae_full = mean(abs(log10(self.pow_orig) - self.fooofed_spectrum_full));
%     self.mse_full = mean((log10(self.pow_orig) - self.fooofed_spectrum_full).^2);
%     self.rmse_full = sqrt(mean((log10(self.pow_orig) - self.fooofed_spectrum_full).^2));
    
    fooof{iFile} = self;
end

%% Plotting
% Frequency range as specified
fig = figure('Units','centimeters','Position',[14.4992    4.9530   33.8243   15.1765]);
tlc = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
for iFile=1:length(subset)
    ax = nexttile;
    plot(fooof{iFile}.freq,fooof{iFile}.pow,'k');
    hold on
    plot(fooof{iFile}.freq,fooof{iFile}.fooofed_spectrum,'r');
    hold on
    plot(fooof{iFile}.freq,fooof{iFile}.ap_fit,'--b')
    % legend ({'Original spectrum','Aperiodic fit','Full model fit'});
    xlabel('Frequency');
    ylabel('logPower');
    str = strsplit(subset(iFile).name,'_');
    title(str{1});
    str1 = sprintf('r2 = %0.3f ',fooof{iFile}.r2);
    str2 = sprintf('mae = %0.3f',fooof{iFile}.error_mae);
    text(ax,'Units', 'Normalized', 'Position', [0.05, 0.1],'String',{str1,str2});

end
title(tlc,sprintf('Model %d - %d Hz',self.freq_range(1),self.freq_range(2)));

% % How does the model fitted in the freq_range extrapolates to the full spectrum?
% fig = figure('Units','centimeters','Position',[14.4992    4.9530   33.8243   15.1765]);
% tlc = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
% for iFile=1:length(subset)
%     ax = nexttile;
%     plot(fooof{iFile}.freq_orig,log10(fooof{iFile}.pow_orig),'k');
%     hold on
%     plot(fooof{iFile}.freq_orig,fooof{iFile}.fooofed_spectrum_full,'r');
%     hold on
%     plot(fooof{iFile}.freq_orig,fooof{iFile}.ap_fit_full,'--b')
%     % legend ({'Original spectrum','Aperiodic fit','Full model fit'});
%     xlabel('Frequency');
%     ylabel('logPower');
%     str = strsplit(subset(iFile).name,'_');
%     title(str{1});
%     str1 = sprintf('r2 = %0.3f ',fooof{iFile}.r2_full);
%     str2 = sprintf('mae = %0.3f',fooof{iFile}.mae_full);
%     text(ax,'Units', 'Normalized', 'Position', [0.05, 0.1],'String',{str1,str2});
% end
% title(tlc,sprintf('Model %d - %d Hz, fitting 1 - 100 Hz',self.freq_range(1),self.freq_range(2)));

% %% Log-log plot
% fig = figure('Units','centimeters','Position',[14.4992    4.9530   33.8243   15.1765]);
% tlc = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
% for iFile=1:length(subset)
%     ax = nexttile;
%     plot(log10(fooof{iFile}.freq_orig),log10(fooof{iFile}.pow_orig),'k');
%     hold on
%     plot(log10(fooof{iFile}.freq_orig),fooof{iFile}.fooofed_spectrum_full,'r');
%     hold on
%     plot(log10(fooof{iFile}.freq_orig),fooof{iFile}.ap_fit_full,'--b')
%     % legend ({'Original spectrum','Aperiodic fit','Full model fit'});
%     xlabel('logFrequency');
%     ylabel('logPower');
%     str = strsplit(subset(iFile).name,'_');
%     title(str{1});
% end
% title(tlc,sprintf('Model %d - %d Hz, fitting 1 - 100 Hz',self.freq_range(1),self.freq_range(2)));

% %% More plotting
% % Overlap the fooofed spectrum to the original spectrum in full frequency range
% padded_ap_fit = nan(size(freq_orig));
% padded_ap_fit(f) = self.ap_fit;
% padded_fooofed_spectrum = nan(size(freq_orig));
% padded_fooofed_spectrum(f) = self.fooofed_spectrum;
% 
% fig = figure('Units','centimeters','Position',[15 15 35 16]);
% tlc = tiledlayout(1,2);
% nexttile;
% plot(freq_orig,log10(pow_orig),'k');
% hold on
% plot(freq_orig,padded_ap_fit,'--b')
% hold on
% plot(freq_orig,padded_fooofed_spectrum,'r');
% legend ({'Original spectrum','Full model fit','Aperiodic fit'});
% xlabel('Frequency');
% ylabel('logPower');
% title('Model fit 3 - 40 Hz');
% str1 = sprintf('r2 = %0.3f ',self.r2);
% str2 = sprintf('mae = %0.3f',self.error_mae);
% annotation(fig,'textbox',[0.15 0.05 .3 .3],'String',{str1,str2},'FitBoxToText','on');
% 
% 
% nexttile
% plot(freq_orig,log10(pow_orig),'k');
% hold on
% plot(freq_orig,self.ap_fit_full,'--b')
% hold on
% plot(freq_orig,self.fooofed_spectrum_full,'r');
% legend ({'Original spectrum','Full model fit','Aperiodic fit'});
% xlabel('Frequency');
% ylabel('logPower');
% title('Model fit 3 - 40 Hz');
% str1 = sprintf('r2 = %0.3f ',r2);
% str2 = sprintf('mae = %0.3f',mae);
% annotation(fig,'textbox',[0.55 0.05 .3 .3],'String',{str1,str2},'FitBoxToText','on');