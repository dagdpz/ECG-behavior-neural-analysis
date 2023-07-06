function [out_norm] = ecg_bna_compute_shufflePredictor_normalization_general(variable,ecg_bna_cfg)
% ecg_bna_compute_shufflePredictor_normalization_general - normalizing the real tfs
% and evoked data based on the shuffle predictor results
%
% USAGE:
%	[out_norm] = ecg_bna_compute_shufflePredictor_normalization(variable,ecg_bna_cfg)
%
% INPUTS:
%		variable  	    - 1xN struct containing the real and shuffled data
%       ecg_bna_cfg     - struct containing configuration settings
%		
% OUTPUTS: 
%       out_norm        - 1xN struct containing the normalized results of
%       pow, itpc, lfp, or itpcbp data
% ======================================================================= %


method = ecg_bna_cfg.shuffle_normalization_method;

if strcmp(method , 'subtraction')
    out_norm    = variable.real.mean - variable.shuffled_mean.mean;
elseif strcmp(method , 'division')
    out_norm    = variable.real.mean / variable.shuffled_mean.mean;
elseif strcmp(method , 'zscore')
    out_norm    = (variable.real.mean - variable.shuffled_mean.mean) / variable.shuffled_mean.std;
elseif strcmp(method , 'not normalized')
    out_norm    = variable.real.mean ;
end

