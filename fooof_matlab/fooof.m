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
        function self = fooof(varargin)
            if nargin > 0
                p = inputParser;
                addParameter(p,'freq_range',[2, 40],@(x)validateattributes(x,{'numeric'},{'numel',2}));
                addParameter(p,'peak_width_limits',[0.5, 12],@(x)validateattributes(x,{'numeric'},{'numel',2}));
                addParameter(p,'max_n_peaks',inf,@(x)validateattributes(x,{'numeric'},{'scalar'}));
                addParameter(p,'min_peak_height',0,@(x)validateattributes(x,{'numeric'},{'scalar'}));
                addParameter(p,'peak_threshold',2,@(x)validateattributes(x,{'numeric'},{'scalar'}));
                addParameter(p,'aperiodic_mode','fixed', @(x) any(validatestring(x,{'fixed','knee'})));

                parse(p,varargin{:})
                
                self.freq_range = p.Results.freq_range;
                self.peak_width_limits = p.Results.peak_width_limits;
                self.max_n_peaks = p.Results.max_n_peaks;
                self.min_peak_height = p.Results.min_peak_height;
                self.peak_threshold = p.Results.peak_threshold;
            end
        end
        
        function self = add_data(self,freq,pow,varargin)
            % Add data to fooof model. Powspectrum must be in linear scale.
            p = inputParser;
            addRequired(p,'freq',@(x)validateattributes(x,{'numeric'},{'vector'}));
            addRequired(p,'pow',@(x)validateattributes(x,{'numeric'},{'vector'}));

            parse(p,varargin{:})
            
        end
        
    end
    
    
    
end