% Auxiliary function to test the python code
% Loads computed power, extracts power in PFC and saves it 

addpath('C:\Users\Mitarbeiter\fieldtrip');
ft_defaults;
load('../params.mat');

power_path = '../../../results/features/power';
mkdir(fullfile(power_path,'PFC'));

files = dir(fullfile(power_path,'*.mat'));
n = length(files);

for iRec = 1:n
    load(fullfile(files(iRec).folder,files(iRec).name));
    % Extract power at the PFC and average across ROI virtual channels
    cfg = [];
    cfg.channel = find(params.PFC_mask);
    power_PFC = ft_selectdata(cfg,power);
    freq = power_PFC.freq;
    pow = power_PFC.powspctrm;
    save(fullfile(power_path,'PFC',files(iRec).name),'freq','pow');
end