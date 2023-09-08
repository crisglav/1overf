function self = fit_fooof(self)
% Fit FOOOF model to power spectrum

% ===== Robust aperiodic fit =====
self.aperiodic_params = robust_ap_fit(self,self.freq,self.pow);
self.ap_fit = exp_function(self.aperiodic_params, self.freq);
self.initial_ap_fit = self.ap_fit; % For visualization later on because it is overwritten

% ===== Flatten the power spectrum using fit aperiodic fit =====
self.spectrum_flat = self.pow - self.ap_fit;
self.inital_spectrum_flat = self.spectrum_flat; % For visualization

% ===== Find peaks and fit them with gaussians ======
self.gaussian_params = fit_peaks(self,self.spectrum_flat);

% Calculatte the peak fit
gaussian_params = reshape(self.gaussian_params',1,[]);
self.peak_fit = gaussian_function(gaussian_params,self.freq);

% ===== Create peak-removed (but not flattened) power spectrum =====
self.spectrum_peak_rm = self.pow - self.peak_fit;

% ===== Run the final aperiodic fit on peak-removed power spectrum ======
% This overwrites the previous aperiodic fit, and recomputes flattened spectrum
self.aperiodic_params = simple_ap_fit(self,self.freq,self.spectrum_peak_rm);
self.ap_fit = exp_function(self.aperiodic_params, self.freq);
self.spectrum_flat = self.pow - self.ap_fit;

% ===== Create full power spectrum model ====
self.fooofed_spectrum = self.peak_fit + self.ap_fit;

% ===== Convert gaussian definitions to peak parameters =====
self.peak_params = create_peak_params(self,self.gaussian_params);

% ===== Calculate R^2 and error of the model fit ======
self.r2 = calc_r_squared(self);
[self.error_mae, self.error_mse, self.error_rmse] = calc_error(self);
end