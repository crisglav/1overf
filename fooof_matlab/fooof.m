classdef fooof
    
    properties(Access = public)
        % Parameters (with default settings)
        freq_range        = [2 40];
        peak_width_limits = [0.5 12];
        max_n_peaks       = inf;
        min_peak_height   = 0;
        peak_threshold    = 2;
        aperiodic_mode    = 'fixed';
        
        % Data attributes
        pow_orig
        freq_orig
        power_spectrum
        freqs
        freq_res
        inital_spectrum_flat
        spectrum_flat
        spectrum_peak_rm
        
        % Model component attributes
        initial_ap_fit
        ap_fit
        peak_fit        
        fooofed_spectrum
        
        % Model parameters
        aperiodic_params
        gaussian_params
        peak_params
        error_mae
        error_mse
        error_rmse
        r2
        
    end
    
    properties(Access = private)
        % Aperiodic computation parameters
        ap_percentile_thresh = 0.025;
        
        % Fitting peak parameters
        bw_std_edge = 1;
        gauss_overlap_thresh = 0.75;
        cf_bound = 1.5;
        gauss_std_limits;       
    end
        
    methods
        function obj = fooof(varargin)
            % constructor of fooof object            
            if nargin > 0
                p = inputParser;
                addParameter(p,'freq_range',[2, 40],@(x)validateattributes(x,{'numeric'},{'numel',2}));
                addParameter(p,'peak_width_limits',[0.5, 12],@(x)validateattributes(x,{'numeric'},{'numel',2}));
                addParameter(p,'max_n_peaks',inf,@(x)validateattributes(x,{'numeric'},{'scalar'}));
                addParameter(p,'min_peak_height',0,@(x)validateattributes(x,{'numeric'},{'scalar'}));
                addParameter(p,'peak_threshold',2,@(x)validateattributes(x,{'numeric'},{'scalar'}));
                addParameter(p,'aperiodic_mode','fixed', @(x) any(validatestring(x,{'fixed','knee'})));

                parse(p,varargin{:})
                
                obj.freq_range = p.Results.freq_range;
                obj.peak_width_limits = p.Results.peak_width_limits;
                obj.max_n_peaks = p.Results.max_n_peaks;
                obj.min_peak_height = p.Results.min_peak_height;
                obj.peak_threshold = p.Results.peak_threshold;
                obj.ap_percentile_thresh = 0.025;
            end
        end
        
        function obj = add_data(obj,freq,pow,varargin)
            % Add data to fooof model. Powspectrum must be in linear scale.                      
            % Trim power spectra if freq range is specified. Otherwise keep all the data
            obj.freq_orig = freq;
            obj.pow_orig = pow;
            if ~isempty(varargin)
                range = varargin{1};
                mask = and(freq >= range(1),freq<=range(2));
                obj.freqs = freq(mask);
                obj.power_spectrum = log10(pow(mask));
            else
                obj.freqs = freq;
                obj.power_spectrum = log10(pow);
            end
            obj.freq_res = freq(2) - freq(1);

        end
        
        function obj = fit(obj)
            % Fit FOOOF model to power spectrum
            
            % ===== Robust aperiodic fit =====
            obj.aperiodic_params = robust_ap_fit(obj,obj.freqs,obj.power_spectrum);
            obj.ap_fit = exp_function(obj.aperiodic_params, obj.freqs);
            obj.initial_ap_fit = obj.ap_fit; % For visualization later on because it is overwritten
            
            % ===== Flatten the power spectrum using fit aperiodic fit =====
            obj.spectrum_flat = obj.power_spectrum - obj.ap_fit;
            obj.inital_spectrum_flat = obj.spectrum_flat; % For visualization
            
            % ===== Find peaks and fit them with gaussians ======
            obj.gaussian_params = fit_peaks(obj,obj.spectrum_flat);
            
            % Calculatte the peak fit
            gaussian_params = reshape(obj.gaussian_params',1,[]);
            obj.peak_fit = gaussian_function(gaussian_params,obj.freqs);
            
            % ===== Create peak-removed (but not flattened) power spectrum =====
            obj.spectrum_peak_rm = obj.power_spectrum - obj.peak_fit;
            
            % ===== Run the final aperiodic fit on peak-removed power spectrum ======
            % This overwrites the previous aperiodic fit, and recomputes flattened spectrum
            obj.aperiodic_params = simple_ap_fit(obj,obj.freqs,obj.spectrum_peak_rm);
            obj.ap_fit = exp_function(obj.aperiodic_params, obj.freqs);
            obj.spectrum_flat = obj.power_spectrum - obj.ap_fit;
            
            % ===== Create full power spectrum model ====
            obj.fooofed_spectrum = obj.peak_fit + obj.ap_fit;
            
            % ===== Convert gaussian definitions to peak parameters =====
            obj.peak_params = create_peak_params(obj,obj.gaussian_params);
            
            % ===== Calculate R^2 and error of the model fit ======
            obj.r2 = calc_r_squared(obj);
            [obj.error_mae, obj.error_mse, obj.error_rmse] = calc_error(obj);
        end
        
        
        function aperiodic_params = simple_ap_fit(obj,freq,pow)
            % Fit the aperiodic component of the power spectra
        
            % Guess initial values
            off_guess = pow(1);
            knee_guess = 0;
            exp_guess = abs((pow(end)-pow(1))/(log10(freq(end))-log10(freq(1))));
            
            % Fitting options
            options = optimoptions('lsqcurvefit');
            options.MaxFunctionEvaluations = 5000;
            options.FunctionTolerance = 1e-8;
            switch obj.aperiodic_mode
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

        function aperiodic_params = robust_ap_fit(obj,freq,pow)
            % Do initial robust aperidoic fit
            guess = simple_ap_fit(obj,freq,pow);
            initial_fit = exp_function(guess, freq);
            
            % Flatten powe spectrum based on initial aperiodic fit
            flatspec = pow - initial_fit;
            
            % Flatten outliers, defined as any points that drop below 0
            flatspec(flatspec < 0) = 0;
            
            % Use percentile threshold, intermos of points to extract and re-fit
            perc_thresh = prctile(flatspec, obj.ap_percentile_thresh);
            perc_mask = flatspec <= perc_thresh;
            freqs_ignore = freq(perc_mask);
            spectrum_ignore = pow(perc_mask);
            
            % Second aperiodic fit, using results of first fit as guess parameters
            options = optimoptions('lsqcurvefit');
            options.MaxFunctionEvaluations = 5000;
            options.FunctionTolerance = 1e-8;
            switch obj.aperiodic_mode
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
        
        function gaussian_params = fit_peaks(obj,flat_iter)
            obj.gauss_std_limits = obj.peak_width_limits/2;
            if(isinf(obj.max_n_peaks))
                guess = zeros(100,3);
            else
                guess = zeros(obj.max_n_peaks,3);
            end
            
            lguess = 1;
            while lguess <= obj.max_n_peaks
                
                % Find candidate peak - the maximum point of the flattened spectrum
                [max_height, max_ind] = max(flat_iter);
                
                % Stop searching for peaks once height drops below height threshold
                if max_height <= obj.peak_threshold * std(flat_iter,1)
                    lguess = lguess -1;
                    break
                end
                
                % set the guess parameters for gaussian fitting, pecifying the mean and
                % height
                guess_freq = obj.freqs(max_ind);
                guess_height = max_height;
                
                % Stop fitting process if candidate peak drops below minimum height
                if not(guess_height > obj.min_peak_height)
                    lguess = lguess -1;
                    break
                end
                
                % Data driven first guess at standard deviation
                % Find half heiht index on each side of the center frequency
                half_height = 0.5 * max_height;
                
                le_ind = find(flat_iter(1:max_ind) <= half_height,1,'last');
                ri_ind = max_ind + find(flat_iter(max_ind:end) <= half_height,1) -1;
                if isempty(le_ind), le_ind = nan; end
                if isempty(ri_ind), ri_ind = nan; end
                
                try
                    % Guess the bandwidth of the peak based on the short side
                    short_side = min(max_ind - le_ind, ri_ind - max_ind);
                    
                    % Full width half max
                    fwhm = short_side * 2 * obj.freq_res;
                    guess_std = fwhm/(2*sqrt(2*log(2)));
                catch
                    guess_std = mean(obj.peak_width_limits);
                end
                
                % Check that the guess value isn't outside preset limits - restrict if
                % so
                if guess_std < obj.gauss_std_limits(1)
                    guess_std = obj.gauss_std_limits(1);
                end
                if guess_std > obj.gauss_std_limits(2)
                    guess_std = obj.gauss_std_limits(2);
                end
                
                % Collect guess parameters and substract this guess gaussian from thee
                % data
                guess(lguess,:) = [guess_freq, guess_height, guess_std];
                
                % Fit the gaussian function
                g = guess_height * exp(-(obj.freqs - guess_freq).^2 / (2*guess_std^2));
                
                % Substract the modelled gaussian from the power spectrum
                flat_iter = flat_iter - g;
                lguess = lguess+1;
            end
            if(lguess > obj.max_n_peaks)
                lguess = obj.max_n_peaks;
            end
            % Remove rows without peaks
            guess = guess(1:lguess,:);
            % Check the peaks based on edges, and on overlap, dropping any that violate
            % requirements
            guess = drop_peak(obj,guess);
            guess = drop_peak_overlap(obj,guess);
            
            if ~isempty(guess)
                gaussian_params = fit_peak_gauss(obj,guess);
                [~, ix] = sort(gaussian_params(:,1));
                gaussian_params = gaussian_params(ix,:);
            else
                gaussian_params = [];
            end
            
        end
        
        function guess = drop_peak(self,guess)
            % Center frequency, Bandwidth
            cf = guess(:,1);
            bw = guess(:,3) * self.bw_std_edge;
            
            % Check if peaks within drop threshold from the edge of the freq range
            keep_peaks = and(abs(cf - self.freq_range(1)) > bw, abs(cf - self.freq_range(2)) > bw);
            
            % Remove peaks that do not fullfil the frequency edge criterion
            guess = guess(keep_peaks,:);
        end
        
        function guess = drop_peak_overlap(self,guess)
            % Sort the peak guesses by increasing frequency
            [~, ix] = sort(guess(:,1));
            guess = guess(ix,:);
            
            % Calculate the std bounds for checking the amount of overlap
            % The bounds are the gaussian frequency +/- gaussian standard deviation
            bounds = [guess(:,1) - guess(:,3)*self.gauss_overlap_thresh, guess(:,1) + guess(:,3)*self.gauss_overlap_thresh] ;
            
            % Loop through peak bounds, comparing current bound to that of the next
            % peak. If the left peak's upper bound extends pass the reight peaks
            % lower bound, then drop the gaussian with the lower height.
            drop_ind = zeros(1,size(bounds,1));
            for i = 1:size(bounds,1)-1
                if bounds(i,2) > bounds(i+1,1)
                    
                    [~,ix] = min([guess(i,2),guess(i+1,2)]);
                    drop_ind(i+ix-1) = 1;
                end
            end
            guess = guess(~drop_ind,:);
            
        end
        
        function gauss_params = fit_peak_gauss(self,guess)
            n_peaks = size(guess,1);
            % Set the bounds for CF
            lo_bound = [guess(:,1) - 2*self.cf_bound *guess(:,3), zeros(n_peaks,1), repmat(self.gauss_std_limits(1),n_peaks,1)];
            hi_bound = [guess(:,1) + 2*self.cf_bound*guess(:,3),inf(n_peaks,1),repmat(self.gauss_std_limits(2),n_peaks,1)];
            
            % Check that CF bounds are within frequency range
            ix = lo_bound(:,1) < self.freq_range(1);
            lo_bound(ix,1) = self.freq_range(1);
            ix = hi_bound(:,1) > self.freq_range(2);
            hi_bound(ix,1) = self.freq_range(2);
            
            % Fit the multiple gaussians
            beta0 = reshape(guess',1,[]); % Starting point
            lb = reshape(lo_bound',1,[]); % Lower bound
            up = reshape(hi_bound',1,[]); % Upper bound
            options = optimoptions('lsqcurvefit');
            options.MaxFunctionEvaluations = 5000;
            options.FunctionTolerance = 1e-8;
            beta = lsqcurvefit(@gaussian_function,beta0,self.freqs',self.spectrum_flat',lb,up,options);
            gauss_params = reshape(beta,3,[])';
        end
        
        function peak_params = create_peak_params(obj,gaussian_params)
            peak_params = zeros(size(gaussian_params));
            n_peaks = size(gaussian_params,1);
            for i=1:n_peaks
                % Get the index of the power spectrum at the frequency closest to the cf of the peak
                [~,ix]= min(abs(obj.freqs - gaussian_params(i,1)));
                peak_params(i,1) = gaussian_params(i,1);
                peak_params(i,2) = obj.fooofed_spectrum(ix) - obj.ap_fit(ix);
                peak_params(i,3) = gaussian_params(i,3)*2;
            end
        end
        
        function [mae, mse, rmse] = calc_error(obj)
            % Mean absolute error
            mae = mean(abs(obj.power_spectrum - obj.fooofed_spectrum));
            % Mean squared error
            mse = mean((obj.power_spectrum - obj.fooofed_spectrum).^2);
            % Root mean squared error
            rmse = sqrt(mean((obj.power_spectrum - obj.fooofed_spectrum).^2));
        end
        
        function r2 = calc_r_squared(obj)
            r2 = corr(obj.power_spectrum', obj.fooofed_spectrum')^2;
        end

    end

end