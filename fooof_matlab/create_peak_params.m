function peak_params = create_peak_params(self,gaussian_params)

peak_params = zeros(size(gaussian_params));
n_peaks = size(gaussian_params,1);
for i=1:n_peaks
    % Get the index of the power spectrum at the frequency closest to the cf of the peak
    [~,ix]= min(abs(self.freq - gaussian_params(i,1)));
    peak_params(i,1) = gaussian_params(i,1);
    peak_params(i,2) = self.fooofed_spectrum(ix) - self.ap_fit(ix);
    peak_params(i,3) = gaussian_params(i,3)*2;
    
end
end