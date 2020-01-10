function [ site_lfp ] = ecg_bna_get_ECG_peaks( site_lfp, block_ECG )
% ecg_bna_get_ECG_peaks - get the ECG Rpeak times for all trials of a
% session and combine it with LFP information in the incoming struct
%
% USAGE:
%	[ site_lfp ] = ecg_bna_get_ECG_peaks( site_lfp, block_ECG )
%
% INPUTS:
%		site_lfp      	- struct containing LFP data for all trials of a
%		session from a single site
%       block_ECG       - struct containing Rpeak information for a block
%       (this usually comes from
%       Y:\Projects\PhysiologicalRecording\Data\[Monkey]\[Date]\[Date]_ecg.mat
%
% OUTPUTS:
%		site_lfp     	- same as input structure but with additional
%		fields to store the ECG Rpeaks for all trials in a block
% 
% See also ecg_bna_process_combined_LFP_ECG, ecg_bna_get_ECG_raw 

%loop through each run
for b = (unique([site_lfp.trials.block]))
    fprintf('Extracting ECG for block %g\n-----------------------\n', b);
    % concatenate all trials for this run
    concat_time = []; % to concatenate sample time
    block_LFP = [];
    trials_time = [];
    trials_idx = find([site_lfp.trials.block] == b);
    for t = 1:length(trials_idx) % ignore first trial
        
        ts = site_lfp.trials(trials_idx(t)).tsample;
        if ~isempty(concat_time) 
            trial_timestamps = ...
                linspace(ts, ts*length(site_lfp.trials(trials_idx(t)).time), ...
                length(site_lfp.trials(trials_idx(t)).time)) + concat_time(end);
        else
            trial_timestamps = ...
                linspace(0, ts*(length(site_lfp.trials(trials_idx(t)).time)-1), ...
                length(site_lfp.trials(trials_idx(t)).time));
        end

        concat_time = [concat_time trial_timestamps];
        trials_time = [trials_time; [trial_timestamps(1), trial_timestamps(end)]];
        
    end
    
    % get ECG timestamps for this block
    ECG_timestamps = block_ECG.out(b).Rpeak_t;    
    ECG_peaksamples = round(ECG_timestamps/ts) + 1;
    trials_samples = round(trials_time / ts) + 1;
    
    % ECG spikes based on ECG timestamps
    ECG_spikes = false(size(concat_time));
    ECG_spikes(ECG_peaksamples) = true;
               
    % now divide into trials
    for t = 1:length(trials_idx)
        
        trial_ECG_spikes = ECG_spikes(...
            concat_time >= trials_time(trials_idx == trials_idx(t), 1) & ...
            concat_time <= trials_time(trials_idx == trials_idx(t), 2));
        site_lfp.trials(trials_idx(t)).ECG_spikes = trial_ECG_spikes;

    end
        
end



end
