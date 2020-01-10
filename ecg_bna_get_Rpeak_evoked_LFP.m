function ecg_triggered_evoked = ecg_bna_get_Rpeak_evoked_LFP( trials_lfp, state )
% ecg_bna_get_Rpeak_evoked_LFP - computes the Rpeak evoked LFP for a specified 
% time window around the Rpeak onset for given trials (usually trials
% belonging to a condition) in a session 
%
% USAGE:
%	ecg_triggered_evoked = ecg_bna_get_Rpeak_evoked_LFP( trials_lfp, state )
%
% INPUTS:
%       trials          - 1xN struct containing LFP and Rpeak data of N
%       trials
%       state           - a cell array specifying time window around Rpeak
%       during which evoked response should be obtained
% OUTPUTS:
%		ecg_triggered_evoked  - struct containing Rpeak onset triggered
%		evoked LFP from the given trials 
%
% See also ecg_bna_get_Rpeak_based_STA, ecg_bna_get_shuffled_Rpeak_evoked_LFP

state_name = state{2};
width = state{4} - state{3};

ecg_triggered_evoked.time = {}; % timestamps
ecg_triggered_evoked.lfp = {}; % evoked LFP response
ecg_triggered_evoked.state = [];
ecg_triggered_evoked.state_name = state_name;

for t = 1:length(trials_lfp)

    trialperiod           = trials_lfp(t).trialperiod;
    % get the LFP samples and timestamps for the trial
    lfp_data = trials_lfp(t).lfp_data(...
        (trials_lfp(t).time >= trialperiod(1) & ...
        trials_lfp(t).time <= trialperiod(2)));
    lfp_time = trials_lfp(t).time(...
        (trials_lfp(t).time >= trialperiod(1) & ...
        trials_lfp(t).time <= trialperiod(2)));
    ecg_peaks = trials_lfp(t).ECG_spikes(...
        (trials_lfp(t).time >= trialperiod(1) & ...
        trials_lfp(t).time <= trialperiod(2)));
    % lfp sample indices
    lfp_sample_idx = 1:length(lfp_time);
    % sample time
    lfp_ts = 1/trials_lfp(t).fsample;
    
    % number of samples in each window
    nsamples_window = round(width/lfp_ts);
      
    % now get the windows to combine
    window_mid_idx = find(ecg_peaks);

    % loop through each window
    for w = 1:length(window_mid_idx)
        if window_mid_idx(w) - round(nsamples_window/2) < 1 || ...
                window_mid_idx(w) + round(nsamples_window/2) > length(lfp_time)
            continue;
        end
        % evoked LFP for this state
        ecg_triggered_evoked.lfp = [ecg_triggered_evoked.lfp, ...
            lfp_data(window_mid_idx(w) - round(nsamples_window/2):...
            window_mid_idx(w) + round(nsamples_window/2))];
        % timestamps
        ecg_triggered_evoked.lfp_time = ...
            lfp_time(window_mid_idx(w) - round(nsamples_window/2):...
            window_mid_idx(w) + round(nsamples_window/2)) ;
        % set mid-timestamp to zero
        ecg_triggered_evoked.lfp_time = ecg_triggered_evoked.lfp_time - ...
            ecg_triggered_evoked.lfp_time(round(length(ecg_triggered_evoked.lfp_time)/2));
        ecg_triggered_evoked.time = [ecg_triggered_evoked.time, ...
            ecg_triggered_evoked.lfp_time];
    end       
    
    
end


if ~isempty(ecg_triggered_evoked.lfp)

    % crop each lfp to same number of samples
    nsamples = min(cellfun('length', ecg_triggered_evoked.lfp));
    for k = 1:length(ecg_triggered_evoked.lfp)
        ecg_triggered_evoked.lfp{k} = ecg_triggered_evoked.lfp{k}(1:nsamples);
        ecg_triggered_evoked.lfp_time = ecg_triggered_evoked.time{k}(1:nsamples);
    end
    ecg_triggered_evoked.lfp_time = ecg_triggered_evoked.lfp_time(1:nsamples);
    
    % evoked LFP average
    arr_state_lfp = vertcat(ecg_triggered_evoked.lfp{:});
    ecg_triggered_evoked.mean = nanmean(arr_state_lfp, 1);
    ecg_triggered_evoked.std = nanstd(arr_state_lfp, 0, 1);
    
end

end

