function aperiodic_params = robust_ap_fit(self,freq,pow)

% Do initial aperidoic fit
guess = simple_ap_fit(self,freq,pow);
initial_fit = exp_function(guess, freq);
% initial_fit = aperiodic_function(self,guess, freq); 

% Flatten powe spectrum based on initial aperiodic fit
flatspec = pow - initial_fit;

% Flatten outliers, defined as any points that drop below 0
flatspec(flatspec < 0) = 0;

% Use percentile threshold, intermos of points to extract and re-fit
perc_thresh = prctile(flatspec, self.ap_percentile_thresh);
perc_mask = flatspec <= perc_thresh;
freqs_ignore = freq(perc_mask);
spectrum_ignore = pow(perc_mask);

% Second aperiodic fit, using results of first fit as guess parameters
options = optimoptions('lsqcurvefit');
options.MaxFunctionEvaluations = 5000;
options.FunctionTolerance = 1e-8;
switch self.aperiodic_mode
    case 'knee'
        lo_bound = [-inf, -inf -inf];
        hi_bound = [inf, inf, inf];
    otherwise % 'fixed'
        % Aperiodic bounds (hard coded)
        lo_bound = [-inf, -inf];
        hi_bound = [inf, inf];
end
aperiodic_params = lsqcurvefit(@exp_function,guess,freqs_ignore',spectrum_ignore',lo_bound,hi_bound,options);

end