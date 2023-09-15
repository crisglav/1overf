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
    
    methods
        % constructor
        function fm = fooof(freq_range,peak_width_limits,max_n_peaks,min_peak_height,peak_threshold)
            if nargin > 0
                fm.freq_range = freq_range;
                fm.peak_width_limits = peak_width_limits;
                fm.max_n_peaks = max_n_peaks;
                fm.min_peak_height = min_peak_height;
                fm.peak_threshold = peak_threshold;
            end
        end
        
    end
end