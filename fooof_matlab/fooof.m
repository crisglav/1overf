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
        function fm = fooof(varargin)
            if nargin > 0
                p = inputParser;
                addParameter(p,'freq_range',[2, 40],@(x)validateattributes(x,{'numeric'},{'numel',2}));
                addParameter(p,'peak_width_limits',[0.5, 12],@(x)validateattributes(x,{'numeric'},{'numel',2}));
                addParameter(p,'max_n_peaks',inf,@(x)validateattributes(x,{'numeric'},{'scalar'}));
                addParameter(p,'min_peak_height',0,@(x)validateattributes(x,{'numeric'},{'scalar'}));
                addParameter(p,'peak_threshold',2,@(x)validateattributes(x,{'numeric'},{'scalar'}));
                addParameter(p,'aperiodic_mode','fixed', @(x) any(validatestring(x,{'fixed','knee'})));

                parse(p,varargin{:})
                
                fm.freq_range = p.Results.freq_range;
                fm.peak_width_limits = p.Results.peak_width_limits;
                fm.max_n_peaks = p.Results.max_n_peaks;
                fm.min_peak_height = p.Results.min_peak_height;
                fm.peak_threshold = p.Results.peak_threshold;
            end
        end
        
    end
    
    
    
end