function [normalized] = ecg_bna_compute_shufflePredictor_normalization_general(real,shuffled,ecg_bna_cfg)
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
parameters={'pow','itpc','lfp','itpcbp'};
for p=1:numel(parameters)
    parameter=parameters{p};
    if strcmp(method , 'subtraction')
        normalized.(parameter)    = real.(parameter).mean-shuffled.(parameter).mean;
    elseif strcmp(method , 'division')
        normalized.(parameter)    = real.(parameter).mean./shuffled.(parameter).mean;
    elseif strcmp(method , 'zscore')
        normalized.(parameter)    = real.(parameter).mean-shuffled.(parameter).mean./shuffled.(parameter).std;
    elseif strcmp(method , 'not normalized')
        normalized.(parameter)    = real.(parameter).mean;
    end
end

