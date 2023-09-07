function aperiodic_params = robust_ap_fit(self,freq,pow)

% Do initial aperidoic fit
popt = simple_ap_fit(self,freq,pow);
initial_fit = aperiodic_function(popt(1), popt(2), freq); 

% Flatten powe spectrum based on initial aperiodic fit
flatspec = pow - initial_fit;

% Flatten outliers, defined as any points that drop below 0
flatspec(flatspec < 0) = 0;

% Use percentile threshold, intermos of points to extract and re-fit
perc_thresh = prctile(flatspec, self.ap_percentile_thresh);
perc_mask = flatspec <= perc_thresh;
freqs_ignore = freq(perc_mask);
spectrum_ignore = pow(perc_mask);

% Get bounds for aperiodic fitting, dropping knee bound if not set to fit
% knee - not implemented

% Second aperiodic fit, using results of first fit as guess parameters
f = fittype('aperiodic_function(of,ex,x)');
popt = fit(freqs_ignore',spectrum_ignore',f,'StartPoint',popt);
aperiodic_params = [popt.of, popt.ex];

end