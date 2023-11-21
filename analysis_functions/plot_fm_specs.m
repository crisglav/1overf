function fig = plot_fm_specs(fm,spec,varargin)
spec_id = spec.spec_id;
epoch_length = spec.epoch_length{:};
taper = spec.taper{:};
fooof_range = spec.fooof_range{:};
fooof_knee = spec.fooof_knee{:};
if nargin > 0
    p = inputParser;
    % TO DO: Do extra checks to the freq range (freq2 > freq1)
    addParameter(p,'fig_save',false,@islogical);
    addParameter(p,'file_name','');
    addParameter(p,'file_path','');
    addParameter(p,'loglog',false,@islogical);
    parse(p,varargin{:})
    
    fig_save = p.Results.fig_save;
    file_name = p.Results.file_name;
    file_path = p.Results.file_path;
    loglog = p.Results.loglog;
end

% Overlap the fooofed spectrum to the original spectrum in full frequency range
padded_ap_fit = nan(size(fm.freq_orig));
fmask = and(fm.freq_orig >= fm.freq_range(1),fm.freq_orig<=fm.freq_range(2));
padded_ap_fit(fmask) = fm.ap_fit;
padded_fooofed_spectrum = nan(size(fm.freq_orig));
padded_fooofed_spectrum(fmask) = fm.fooofed_spectrum;

fig = figure('Units','centimeters','Position',[15 10 20 15],'visible','off');
ax = gca;
if loglog
    sp = plot(log10(fm.freq_orig),log10(fm.pow_orig),'k');
    hold on
    ap = plot(log10(fm.freq_orig),padded_ap_fit,'--b');
    hold on
    fo = plot(log10(fm.freq_orig),padded_fooofed_spectrum,'r');
    xlabel('logFrequency');
else
    sp = plot(fm.freq_orig,log10(fm.pow_orig),'k');
    hold on
    ap = plot(fm.freq_orig,padded_ap_fit,'--b');
    hold on
    fo = plot(fm.freq_orig,padded_fooofed_spectrum,'r');
    ylim([-1.5 1.5]);
    xlabel('Frequency');
end
legend ([sp(1),ap,fo],{'Original spectrum','Aperiodic fit','Full model fit'});
ylabel('logPower');

s = sprintf('SPEC %d. EPOCH LENGTH: %s TAPER: %s FREQ RANGE: %s KNEE: %s',spec_id, epoch_length, taper, fooof_range,fooof_knee);
title(s);
% Extract error
if isprop(fm,'group_results')
    r_squared = mean(cellfun(@(x) x.r_squared, fm.group_results)); 
    mae = mean(cellfun(@(x) x.error(1), fm.group_results));
else
    r_squared = fm.r_squared;
    mae = fm.error(1);
end
str1 = sprintf('r2 = %0.3f ',r_squared);
str2 = sprintf('mae = %0.3f',mae);
text('Units', 'Normalized', 'Position', [0.05, 0.1],'String',{str1,str2});

if(fig_save)
    if(isempty(file_name))
        error('Specify a valid file_name');
    end
    exportgraphics(ax,fullfile(file_path,[file_name num2str(spec_id) '.pdf']));
end
end

