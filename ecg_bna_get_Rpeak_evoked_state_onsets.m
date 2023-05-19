function Rpeak_evoked_state = ecg_bna_get_Rpeak_evoked_state_onsets(trials_ecg, state)
% ecg_bna_get_Rpeak_evoked_state_onsets - computes the event onset
% probability for a specified time window around the Rpeak onset for given
% trials (usually trials belonging to a condition) in a session 
%
% USAGE:
%	Rpeak_evoked_state = ecg_bna_get_Rpeak_evoked_state_onsets(trials_ecg,
%	state) 
%
% INPUTS:
%       trials          - 1xN struct containing LFP and Rpeak data of N
%       trials
%       state           - a cell array specifying time window around Rpeak
%       during which evoked response should be obtained
% OUTPUTS:
%		Rpeak_evoked_state  - struct containing event onset probability in
%		a given time window around Rpeak onset triggered for the given
%		trials 
%
% See also ecg_bna_compute_session_Rpeak_evoked_state_onsets

state_id = state{1};
state_name = state{2};
Rpeak_ref_abs_time = [state{3} state{4}];
Rpeak_ref_rel_time = state{5};%state_reftend = state{4};

Rpeak_evoked_state.abs_onset_times = [];
Rpeak_evoked_state.rel_onset_times = [];
Rpeak_evoked_state.valid_trials = [];
Rpeak_evoked_state.state_id = state_id;
Rpeak_evoked_state.state_name = state_name;

for t = 1:length(trials_ecg)
    
    if isempty(trials_ecg(t).ECG_spikes) || isempty(trials_ecg(t).states) 
        continue;
    end

    states          = trials_ecg(t).states;
    s=[states(:).id] == state_id;
    if any(s)
        state_onset_t   = states(s).onset_t;
    else
        continue;
    end
    
%     if ~all(trials_ecg(t).ECG_valid)%(state_onset_sample)
%         continue;
%     end
    
    % Rpeaks
    trial_Rpeaks = trials_ecg(t).ECG_spikes;
    start_time=trials_ecg(t).tstart;
    ts=trials_ecg(t).tsample;
    nsamples = numel(trials_ecg(t).ecg_data);
    end_time = start_time + (ts*(nsamples-1));
    trial_timestamps = linspace(start_time, end_time, nsamples);

    %trial_timestamps = trials_ecg(t).time;
    Rpeak_times = trial_timestamps(trial_Rpeaks);
    Rpeak_ref_onset_times = state_onset_t - Rpeak_times;
        
    Rpeak_rel_onset_time = nan;
    % check if state onset occured after first Rpeak of this trial
    abs_onset_time_after_Rpeak = min(Rpeak_ref_onset_times(Rpeak_ref_onset_times > 0));
    if ~isempty(abs_onset_time_after_Rpeak)
        % find the id of Rpeak after which state onset occured
        Rpeak_idx = find(Rpeak_ref_onset_times == abs_onset_time_after_Rpeak);
        if Rpeak_idx < length(Rpeak_ref_onset_times)
            Rpeak_rel_onset_time = abs_onset_time_after_Rpeak / (Rpeak_times(Rpeak_idx + 1) - Rpeak_times(Rpeak_idx));
        end
    end
    
    abs_onset_time_around_Rpeak = Rpeak_ref_onset_times(Rpeak_ref_onset_times < Rpeak_ref_abs_time(2) & Rpeak_ref_onset_times > Rpeak_ref_abs_time(1));
%     abs_onset_time_before_Rpeak = ...
%         max(Rpeak_ref_onset_times(Rpeak_ref_onset_times < 0 & ...
%         Rpeak_ref_onset_times > Rpeak_ref_abs_time(1)));
%     abs_onset_time_around_Rpeak = [abs_onset_time_before_Rpeak, ...
%         abs_onset_time_after_Rpeak];   
        
%     elseif strcmp(Rpeak_ref_rel_time, 'beforeRpeak')
%         Rpeak_ref_onset_times = Rpeak_times - state_onset_t;
%         abs_onset_time_after_Rpeak = ...
%             max(abs_onset_time_after_Rpeak(abs_onset_time_after_Rpeak < 0));
%         % check if state onset occured before last Rpeak of this trial
%         if ~isempty(abs_onset_time_after_Rpeak)
%             % find the id of Rpeak before which state onset occured
%             Rpeak_idx = find(Rpeak_ref_onset_times == ...
%                 abs_onset_time_after_Rpeak);
%             if Rpeak_idx > 1
%                 Rpeak_rel_onset_time = abs_onset_time_after_Rpeak / ...
%                     (Rpeak_times(Rpeak_idx) - Rpeak_times(Rpeak_idx - 1));
%             end
%         end
%     elseif strcmp(Rpeak_ref_rel_time, 'mindist')
%         abs_onset_time_after_Rpeak = Rpeak_times - state_onset_t;
%         abs_onset_time_after_Rpeak = abs_onset_time_after_Rpeak(abs(abs_onset_time_after_Rpeak) == ...
%             min(abs(abs_onset_time_after_Rpeak))); 
%     end
    
    if ~isempty(abs_onset_time_after_Rpeak)
        Rpeak_evoked_state.rel_onset_times = [Rpeak_evoked_state.rel_onset_times Rpeak_rel_onset_time];
        Rpeak_evoked_state.valid_trials = [Rpeak_evoked_state.valid_trials, t];
    end
    
    if ~isempty(abs_onset_time_around_Rpeak)
        Rpeak_evoked_state.abs_onset_times = [Rpeak_evoked_state.abs_onset_times abs_onset_time_around_Rpeak];
    end
end

% get histogram counts
% absolute time from Rpeak
abs_tbinwidth = 0.025; % move this to settings
abs_tbinedges = Rpeak_ref_abs_time(1) : abs_tbinwidth : Rpeak_ref_abs_time(2);
[Rpeak_evoked_state.abs_histcounts.prob, Rpeak_evoked_state.abs_histcounts.timebins] = ...
    histcounts(Rpeak_evoked_state.abs_onset_times, abs_tbinedges, 'Normalization', 'probability');
% relative time from Rpeak
rel_tbinwidth = 1/36; % move this to settings
rel_tbinedges = 0:rel_tbinwidth:1;
[Rpeak_evoked_state.rel_histcounts.prob, Rpeak_evoked_state.rel_histcounts.timebins] = ...
    histcounts(Rpeak_evoked_state.rel_onset_times, rel_tbinedges, 'Normalization', 'probability');
Rpeak_evoked_state.ntrials = length(Rpeak_evoked_state.rel_onset_times);

end