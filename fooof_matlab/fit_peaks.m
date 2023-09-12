function gaussian_params = fit_peaks(self,flat_iter)

% Fitting peak parameters
self.bw_std_edge = 1;
self.gauss_overlap_thresh = 0.75;
self.cf_bound = 1.5;
self.gauss_std_limits = self.peak_width_limits/2;

freq_res = self.freq(2) - self.freq(1);

if(isinf(self.max_n_peaks))
    guess = zeros(100,3);
else
    guess = zeros(self.max_n_peaks,3);
end

lguess = 1;
while lguess <= self.max_n_peaks
    
    % Find candidate peak - the maximum point of the flattened spectrum
    [max_height, max_ind] = max(flat_iter);
    
    % Stop searching for peaks once height drops below height threshold
    if max_height <= self.peak_threshold * std(flat_iter,1)
        lguess = lguess -1;
        break
    end
    
    % set the guess parameters for gaussian fitting, pecifying the mean and
    % height
    guess_freq = self.freq(max_ind);
    guess_height = max_height;
    
    % Stop fitting process if candidate peak drops below minimum height
    if not(guess_height > self.min_peak_height)
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
        fwhm = short_side * 2 * freq_res;
        guess_std = fwhm/(2*sqrt(2*log(2)));
    catch
        guess_std = mean(self.peak_width_limits);
    end
    
    % Check that the guess value isn't outside preset limits - restrict if
    % so
    if guess_std < self.gauss_std_limits(1)
        guess_std = self.gauss_std_limits(1);
    end
    if guess_std > self.gauss_std_limits(2)
        guess_std = self.gauss_std_limits(2);
    end
    
    % Collect guess parameters and substract this guess gaussian from thee
    % data
    guess(lguess,:) = [guess_freq, guess_height, guess_std];
    
    % Fit the gaussian function
    g = guess_height * exp(-(self.freq - guess_freq).^2 / (2*guess_std^2));
    
    % Substract the modelled gaussian from the power spectrum
    flat_iter = flat_iter - g;
    lguess = lguess+1;
end
if(lguess > self.max_n_peaks)
    lguess = self.max_n_peaks;
end
% Remove rows without peaks
guess = guess(1:lguess,:);
% Check the peaks based on edges, and on overlap, dropping any that violate
% requirements
guess = drop_peak(self,guess);
guess = drop_peak_overlap(self,guess);

if ~isempty(guess)
    gaussian_params = fit_peak_gauss(self,guess);
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
            
    beta0 = reshape(guess',1,[]); % Starting point
    lb = reshape(lo_bound',1,[]); % Lower bound
    up = reshape(hi_bound',1,[]); % Upper bound
    options = optimoptions('lsqcurvefit');
    options.MaxFunctionEvaluations = 5000;
    options.FunctionTolerance = 1e-8;
    beta = lsqcurvefit(@gaussian_function,beta0,self.freq',self.spectrum_flat',lb,up,options);  
    gauss_params = reshape(beta,3,[])';
end