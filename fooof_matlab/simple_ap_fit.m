function aperiodic_params = simple_ap_fit(self,freq,pow)
% The knee version is not implemented

% Guess initial values
off_guess = pow(1);
exp_guess = abs((pow(end)-pow(1))/(log10(freq(end))-log10(freq(1))));

% Fitting
f = fittype('aperiodic_function(of,ex,x)');
popt = fit(freq',pow',f,'StartPoint',[off_guess exp_guess]);
aperiodic_params = [popt.of, popt.ex];
end