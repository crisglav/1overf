classdef fooof
    
    properties(Access = public)
        % All the default settings
        freq_range =        [2 40];
        peak_width_limits = [0.5 12];
        max_n_peaks       = inf; 
        min_peak_height   = 0;
        peak_threshold    = 2;
    end
    
    properties(Access = private)
        ap_percentile_thresh = 0.025;
    end
end