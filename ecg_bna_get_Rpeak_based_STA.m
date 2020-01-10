function ecg_triggered_evoked = ecg_bna_get_Rpeak_based_STA( trials, label, cfg_ecg, signal_type )
% ecg_bna_get_Rpeak_based_STA - computes the Rpeak evoked LFP/ECG for a specified 
% time window around the Rpeak onset for given trials (usually trials
% belonging to a condition) in a session 
%
% USAGE:
%	ecg_triggered_evoked = ecg_bna_get_Rpeak_based_STA( trials, site_ID, cfg_ecg, 'lfp' ) 
%   ecg_triggered_evoked = ecg_bna_get_Rpeak_based_STA( trials, site_ID, cfg_ecg, 'ecg' ) 
%
% INPUTS:
%       trials          - 1xN struct containing LFP and ECG data of N
%       trials
%       label           - a string specifying a label of the signal whose
%       evoked response is computed
%       cfg_ecg         - a cell array specifying time window around Rpeak
%       during which evoked response should be obtained
%       signal_type     - type of signal whose Rpeak triggered evoked
%       response should be obtained. Can be 'lfp' or 'tfs'
% OUTPUTS:
%		ecg_triggered_evoked  - struct containing Rpeak onset triggered
%		evoked LFP/ECG from the given trials
% 
% REQUIRES:	ft_spiketriggeredaverage (FieldTrip toolbox)
%
% See also ecg_bna_get_Rpeak_evoked_LFP, ecg_bna_get_shuffled_Rpeak_evoked_LFP

ecg_triggered_evoked.trial = {};
ecg_triggered_evoked.time = {};
ecg_triggered_evoked.state = cfg_ecg{1};
ecg_triggered_evoked.state_name = cfg_ecg{2};

% create a FT datatype
ft_data_all = struct();
ft_data_all.label = {label, 'ECG'};
ft_data_all.time = {}; % timestamps
ft_data_all.trial = {}; % evoked LFP response

for t = 1:length(trials)

    trialperiod           = trials(t).trialperiod;
    
    % check if ECG spike data exists for this trial
    if isempty(trials(t).ECG_spikes)
        continue;
    end
    
    % get the LFP samples and timestamps for the trial
    if strcmp(signal_type, 'lfp')
        signal = trials(t).lfp_data;
    elseif strcmp(signal_type, 'ecg')
        signal = trials(t).ecg_data;
    end
    
    ecg_peaks = trials(t).ECG_spikes;
    signal_time = trials(t).time;
    ft_data_all.trial = [ft_data_all.trial, ...
        [signal; ecg_peaks]];
    ft_data_all.time = [ft_data_all.time, signal_time];

end

if ~isempty(ft_data_all.trial)

    % evoked LFP average
    cfg                 = [];
    cfg.keeptrials      = 'yes';
    cfg.timwin          = [cfg_ecg{3} cfg_ecg{4}];
    cfg.channel         = ft_data_all.label(1); % LFP chan
    cfg.spikechannel    = ft_data_all.label(2); % ECG peak
    
    ecg_based_sta       = ft_spiketriggeredaverage(cfg, ft_data_all);
    
    ecg_triggered_evoked.trial = ecg_based_sta.trial;
    ecg_triggered_evoked.dimord = 'npeaks_time';
    ecg_triggered_evoked.time = ecg_based_sta.time;
    ecg_triggered_evoked.mean = ecg_based_sta.avg;
    ecg_triggered_evoked.std = permute(nanstd(ecg_based_sta.trial, 0, 1), [2 3 1]);
%     ecg_based_sta.std = nanstd(arr_state_lfp, 0, 1);

    clear ft_data_all ecg_based_sta;
    
end

end

