function [out_norm] = ecg_bna_compute_shufflePredictor_normalization(real_tfs,real_evoked,shuffled_mean,ecg_bna_cfg)
% ecg_bna_compute_shufflePredictor_normalization - normalizing the real tfs
% and evoked data based on the shuffle predictor results
%
% USAGE:
%	[out_norm] = ecg_bna_compute_shufflePredictor_normalization(real_tfs,real_evoked,shuffled_mean,ecg_bna_cfg)
%
% INPUTS:
%		real_tfs  	    - 1xN struct containing the tfs (ITPC and POWER) data
%		real_evoked  	- 1xN struct containing the LFP and ITPCbp data
%		shuffled_mean  	- 1xN struct containing the average of all shuffled data
%       ecg_bna_cfg     - struct containing configuration settings
%		
% OUTPUTS: 
%       out_norm        - 1xN struct containing the normalized results of
%       pow, itpc, lfp, and itpcbp data
% ======================================================================= %


method = ecg_bna_cfg.shuffle_normalization_method;
% name = fieldnames(shuffled_mean);
% for vr = 1:numel(fieldnames(shuffled_mean))

if strcmp(method , 'subtraction')
    out_norm.pow    = real_tfs.pow.mean - shuffled_mean.pow.mean;
    out_norm.itpc   = real_tfs.itpc.mean - shuffled_mean.itpc.mean;
    out_norm.lfp    = real_evoked.lfp.mean - shuffled_mean.lfp.mean;
    out_norm.pow    = real_evoked.itpcbp.mean - shuffled_mean.itpcbp.mean;
elseif strcmp(method , 'division')
    out_norm.pow    = real_tfs.pow.mean / shuffled_mean.pow.mean;
    out_norm.itpc   = real_tfs.itpc.mean / shuffled_mean.itpc.mean;
    out_norm.lfp    = real_evoked.lfp.mean / shuffled_mean.lfp.mean;
    out_norm.pow    = real_evoked.itpcbp.mean / shuffled_mean.itpcbp.mean;
elseif strcmp(method , 'zscore')
    out_norm.pow    = (real_tfs.pow.mean - shuffled_mean.pow.mean) / shuffled_mean.pow.std;
    out_norm.itpc   = (real_tfs.itpc.mean - shuffled_mean.itpc.mean) / shuffled_mean.itpc.std;
    out_norm.lfp    = (real_evoked.lfp.mean - shuffled_mean.lfp.mean) / shuffled_mean.lfp.std;
    out_norm.pow    = (real_evoked.itpcbp.mean - shuffled_mean.itpcbp.mean) / shuffled_mean.itpcbp.std;
elseif strcmp(method , 'not normalized')
    out_norm.pow    = real_tfs.pow.mean ;
    out_norm.itpc   = real_tfs.itpc.mean ;
    out_norm.lfp    = real_evoked.lfp.mean ;
    out_norm.pow    = real_evoked.itpcbp.mean ;
end

