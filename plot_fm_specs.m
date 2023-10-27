function fig = plot_fm_specs(fm,spec,varargin)
spec_id = spec{1};
taper = spec{2};
padding = spec{3};
average_psd = spec{4};
fooof_range = spec{5};
fooof_knee = spec{6};
if nargin > 0
    p = inputParser;
    % TO DO: Do extra checks to the freq range (freq2 > freq1)
    addParameter(p,'fig_save',false,@islogical);
    addParameter(p,'file_name','');
    addParameter(p,'file_path','');
    parse(p,varargin{:})
    
    fig_save = p.Results.fig_save;
    file_name = p.Results.file_name;
    file_path = p.Results.file_path;
end

% Overlap the fooofed spectrum to the original spectrum in full frequency range
padded_ap_fit = nan(size(fm.freq_orig));
fmask = and(fm.freq_orig >= fm.freq_range(1),fm.freq_orig<=fm.freq_range(2));
padded_ap_fit(fmask) = fm.ap_fit;
padded_fooofed_spectrum = nan(size(fm.freq_orig));
padded_fooofed_spectrum(fmask) = fm.fooofed_spectrum;

fig = figure('Units','centimeters','Position',[15 10 20 15]);
ax = gca;
sp = plot(fm.freq_orig,log10(fm.pow_orig),'k');
hold on
ap = plot(fm.freq_orig,padded_ap_fit,'--b');
hold on
fo = plot(fm.freq_orig,padded_fooofed_spectrum,'r');
ylim([-1.5 1.5]);
legend ([sp,ap,fo],{'Original spectrum','Aperiodic fit','Full model fit'});
xlabel('Frequency');
ylabel('logPower');

sprintf('Spec %d: taper: %c padding: %d freq_range: %d - %d knee: ',spec_id, taper,padding, fooof_range,fooof_knee)

title(sprintf('Model fit %d - %d Hz',fm.freq_range));
str1 = sprintf('r2 = %0.3f ',fm.r_squared);
str2 = sprintf('mae = %0.3f',fm.error(1));
text('Units', 'Normalized', 'Position', [0.05, 0.1],'String',{str1,str2});

if(fig_save)
    if(isempty(file_name))
        error('Specify a valid file_name');
    end
    exportgraphics(ax,fullfile(file_path,[file_name '.pdf']));
end
end

