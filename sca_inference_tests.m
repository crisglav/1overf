function [p_median,p_share,p_aggregate] = sca_inference_tests(effect_orig,bf_orig,p_orig,effect_rand,bf_rand,p_rand,tail)
% Compute inference tests based on specification curve analysis
% 
% Input: 
% - effect_orig: effect sizes from the observed specification curve. Vector 1 x nSpec
% - bf_orig: bayes factor from the observed specification curve. Vector 1 x nSpec
% - p_orig: pvalues from the observed specification curve. Vector 1 x nSpec
% - effect_rand: effect sizes from the randomized specification curves. Matrix nSpec x nRand
% - bf_rand: bayes factor from the randomized specification curves. Matrix nSpec x nRand
% - p_rand: pvalues sizes from the randomized specification curves. Matrix nSpec x nRand
% - tail: indicate whether the a priori hypothesis was left sided, right sided or two-sided. String "left", "right", "both"
%
% Output:
% - p_median
% - p_share
% - p_aggregate
%
% Tests based on the paper 'Specification curve analysis',Simonsohn et al. 2020, Nat Hum Behav
% Code by Elisabeth May, Cristina Gil and Felix Bott, Technical University of Munich, 18.12.2023


% Check that the input is as indicated in the documentation

% Number of specifications and randomizations
[nSpec, nRand] = size(effect_rand); 


% Calculate the dominant sign of the original and randomized curves
% ===================================================================
% The dominant sign of the overall specification curve will be negative if 
% the majority of specifications have a negative effect and positive otherwise. 
% If there are the same number of positive and negative specifications the
% dominant sign will be positive.

% Original curve
pos = sum(effect_orig > 0);
if (pos/nSpec) < 0.5
    sign_orig = -1;
else
    sign_orig = 1;
end
% Randomizations
pos = sum(effect_rand > 0);
sign_rand = ones(1,nRand);
sign_rand((pos/nSpec) < 0.5) = -1;
sign_rand((pos/nSpec) == 0.5) = sign_orig; % if there is no dominant sign of effects, assume direction of original effect


% INFERENTIAL TESTS

% 1. TEST OF MEDIAN EFFECT SIZE
% ===================================================================
median_orig = median(effect_orig);
median_rand = median(effect_rand);
switch tail
    case 'left'
        % One-sided test
        % ==============
        % calculate p-value as percentage of randomizations with median
        % effects size smaller than original effect size
        t = sum(median_rand <= median_orig);
    case 'right'
        % One-sided test
        % ==============
        % calculate p-value as percentage of randomizations with median
        % effects size larger than original effect size
        t = sum(median_rand >= median_orig);
    case 'both'
        % Two-sided test
        % ==============
        % calculate p-value as percentage of randomizations with median
        % effects size more extreme than original effect size on both ends
        % of the distribution of median effect sizes
        t = sum(abs(median_rand) >= abs(median_orig));
end
p_median = t/nRand;

% 2. TEST OF SHARE OF 'SIGNIFICANT' RESULTS (HERE SIGNIFICANT: BF > 3)
% ====================================================================
% Extract the effect sizes of the 'significant' specifications for the
% original curve and the randomizations
significant_orig = effect_orig(bf_orig > 3);
significant_rand = nan(nSpec,nRand); 
significant_rand(bf_rand > 3) = effect_rand(bf_rand > 3);

switch tail
    case 'left'
        % One-sided test
        % ==============
        % Count the significant specifications that have a negative sign
        if ~isempty(significant_orig)
            BF_orig_count = sum(significant_orig < 0);
        else
            BF_orig_count = 0;
        end
        BF_rand_count = sum((significant_rand < 0),1);
        % Percentage of randomizations that have a number of negative significant
        % specifications higher than the original numer of negative signficiant
        % specifications
        t = sum(BF_rand_count >= BF_orig_count);

    case 'right'
        % One-sided test
        % ==============
        % Count the significant specifications that have a positive sign
        if ~isempty(significant_orig)
            BF_orig_count = sum(significant_orig > 0);
        else
            BF_orig_count = 0;
        end
        BF_rand_count = sum((significant_rand > 0),1);
        % Percentage of randomizations that have a number of positive significant
        % specifications higher than the original numer of positive signficiant
        % specifications
        t = sum(BF_rand_count >= BF_orig_count);

    case 'both'
        % Two-sided test
        % ==============
        % Count only the significant specifications that have the same sign as the
        % dominant curve effect
        if ~isempty(significant_orig)
            BF_orig_count = sum((significant_orig * sign_orig) > 0);
        else
            BF_orig_count = 0;
        end
        % For each randomization, count the significant specifications
        % that have the same sign as the dominant curve effect of that
        % randomization.
        BF_rand_count = sum((significant_rand .* sign_rand) > 0,1);
        % Percentage of randomizations that have a number of significant specifications higher
        % than the original number of signficiant specifications (either
        % positive or negative)
        t = sum(BF_rand_count >= BF_orig_count);
end
p_share = t/nRand;

% 3. TEST OF AGGREGATED P-VALUES
% =================================================================
% For each randomization, average pvalues following the Stouffer's method.
% First get z-scores for each pvalue and average across specifications
% Note that all versions of the test are conceptually the same, so the
% 'left', 'right' and 'both' options yield the same p-value.

switch tail
    case 'left'
        % One-sided test
        % ==============
        % calculate p-value as percentage of randomizations with average z value
        % smaller than original average z value
        z_temp = norminv(p_orig);
        z_orig = sum(z_temp)/sqrt(nSpec);

        z_temp = norminv(p_rand);
        z_rand = sum(z_temp,1)./sqrt(nSpec);

        t = sum(z_rand <= z_orig);

    case 'right'
        % One-sided test
        % ==============
        % calculate p-value as percentage of randomizations with average z value
        % larger than original average z value. We do 1-pvalue to get positive zscores.
        z_temp = norminv(1-pvalues_orig);
        z_orig = sum(z_temp)/sqrt(nSpec);

        z_temp = norminv(1-pvalues_rand);
        z_rand = sum(z_temp,1)./sqrt(nSpec);

        t = sum(z_rand >= z_orig);

    case 'both'
        % Two-sided test
        % ==============
        % calculate p-value as percentage of randomizations with average z value
        % more extreme than original average z value on both ends of the distribution of z
        % values rand. Pvalues are diveded by two because the original t-tests were 
        % two-sided.  We do 1-pvalue to get positive zscores.
        z_temp = norminv(1-(pvalues_orig./2));
        z_orig = sum(z_temp)/sqrt(nSpec);

        z_temp = norminv(1-(pvalues_rand./2));
        z_rand = sum(z_temp,1)./sqrt(nSpec);

        t = sum(z_rand >= z_orig);

end
p_aggregate = t/nRand;


end