% Check the stability of the estimation of the aperiodic exponents
%
% Cristina Gil, TUM, 20.02.2024


% Load randomizations
sca_path = '/rechenmagd3/Experiments/2023_1overf/results/sca';
nRand = 500;

% Load original specification curve
load(fullfile(sca_path,'specs_ap_exp.mat'));

% Keep the original (correct) masks of patients and healthy participants
hc_mask_orig = hc_mask;
pa_mask_orig = pa_mask;
exp_orig = exp;

% Median aperiodic exponents of the original curve
exp_hc = exp(:,hc_mask_orig);
exp_pa = exp(:,pa_mask_orig);

% Median of the the first specification (main analysis), no randomization
median_hc_orig = median(exp_hc(1,:));
median_pa_orig = median(exp_pa(1,:));
    
% Median of the the first specification (main analysis), all randomizations
std_hc = nan(1,nRand);
std_pa = nan(1,nRand);
exp_all = nan(48,264,nRand);

for iRand=1:nRand
    load(fullfile(sca_path,sprintf("specs_ap_exp_rand%0.3d.mat",iRand)));
    exp_all(:,:,iRand) = exp;
    
    exp_hc = exp(:,hc_mask_orig);
    exp_pa = exp(:,pa_mask_orig);
    
    % Median of the the first specification (main analysis)
    std_hc(iRand) = std(exp_hc(1,:));
    std_pa(iRand) = std(exp_pa(1,:));

end



