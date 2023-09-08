function aperiodic_params = simple_ap_fit(self,freq,pow)

% Guess initial values
off_guess = pow(1);
knee_guess = 0;
exp_guess = abs((pow(end)-pow(1))/(log10(freq(end))-log10(freq(1))));

% Fitting options
options = optimoptions('lsqcurvefit');
options.MaxFunctionEvaluations = 5000;
options.FunctionTolerance = 1e-8;
switch self.aperiodic_mode
    case 'knee'
        lo_bound = [-inf, -inf -inf];
        hi_bound = [inf, inf, inf];
        guess = [off_guess, knee_guess, exp_guess];
    otherwise
        lo_bound = [-inf, -inf];
        hi_bound = [inf, inf];
        guess = [off_guess exp_guess];        
end
% Fitting
aperiodic_params = lsqcurvefit(@exp_function,guess,freq',pow',lo_bound,hi_bound,options);
% The fitting can sometimes return aperiodic params in form of a complex number with imaginary part equal to zero
aperiodic_params = real(aperiodic_params);
end