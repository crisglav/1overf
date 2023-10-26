classdef fooofGroup < fooof
    properties
        power_spectra
        group_results
    end
    
    methods
        function obj = fooofGroup(varargin)
            % constructor of fooofGroup object. It is derived from class fooof
            obj@fooof(varargin{:});
        end     
        
        function obj = add_data(obj,freq,pow,varargin)
            % Add data to fooof model. Powspectrum must be in linear scale.                      
            % Trim power spectra if freq range is specified. Otherwise keep all the data
            if size(pow,1) == 1
                warning('Only one power spectrum given. Use a fooof object.')
            end
            obj.freq_orig = freq;
            obj.pow_orig = pow;
            if ~isempty(varargin)
                range = varargin{1};
                mask = and(freq >= range(1),freq<=range(2));
                obj.freqs = freq(mask);
                obj.power_spectra = log10(pow(:,mask));
            else
                obj.freqs = freq;
                obj.power_spectra = log10(pow);
            end
            obj.freq_res = freq(2) - freq(1);
            obj.group_results = cell(1,size(obj.power_spectra,1));

        end
        
        function obj = fit(obj)
            for i=1:size(obj.power_spectra,1)
                obj.power_spectrum = obj.power_spectra(i,:);
                obj = fit@fooof(obj);
                obj.group_results{i} = get_results(obj);  
            end
        end
        
        function results = get_results(obj)
            results = get_results@fooof(obj);
        end
        
    end
end

