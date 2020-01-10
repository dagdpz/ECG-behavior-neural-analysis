function shuffled_LFP_evoked = ecg_bna_get_shuffled_Rpeak_evoked_LFP( trials_lfp, cfg_state, nshuffles )
% ecg_bna_get_shuffled_Rpeak_evoked_ECG - computes the evoked LFP for
% a specified  time window around randomly shuffled triggers for given
% trials (usually trials belonging to a condition) in a session. The random 
% triggers are obtained by randomly shuffling the Rpeaks based on their 
% peak-to-peak intervals 
%
% USAGE:
%	shuffled_ecg_evoked = ecg_bna_get_shuffled_Rpeak_evoked_ECG(
%	trials_ecg, cfg_ecg, nshuffles ) 
%
% INPUTS:
%       trials_lfp      - 1xN struct containing LFP and Rpeak data of N
%       trials
%       cfg_state       - a cell array specifying time window around the
%       trigger during which evoked response should be obtained
%       nshuffles       - number of times the Rpeak triggers should be 
%       randomly shuffled
% OUTPUTS:
%		shuffled_LFP_evoked   - struct containing shuffled Rpeak triggered
%		evoked LFP from the given trials (mean and standard deviation
%		across shuffles)
%
% REQUIRES:	ft_spiketriggeredaverage (FieldTrip toolbox)
%
% See also ecg_bna_compute_session_Rpeak_evoked_LFP, ecg_bna_get_Rpeak_evoked_LFP, ecg_bna_get_shuffled_Rpeak_evoked_ECG

% whether to plot the distribution of ECG peak-to-peak interval
distplot = false;

shuffled_LFP_evoked.lfp = {};
shuffled_LFP_evoked.lfp_time = {};
shuffled_LFP_evoked.state = cfg_state{1};
shuffled_LFP_evoked.state_name = cfg_state{2};

LFP_raw = [];
%timestamps = [];
ecg_peaks = [];
ecg_R2Rt = [];

for t = 1:length(trials_lfp)
    if trials_lfp(t).completed == 0
        continue;
    end
    % check if ECG spike data exists for this trial
    if isempty(trials_lfp(t).ECG_spikes)
        continue;
    end
    % get the ECG raw data for the trial
    if isempty(trials_lfp(t).lfp_data)
        continue;
    else
        trial_LFP_raw = trials_lfp(t).lfp_data;
        %trial_LFP_raw(~trials_lfp(t).ECG_valid) = [];
        LFP_raw = [LFP_raw, trial_LFP_raw];
    end
    
    %timestamps = [timestamps, trials_ecg(t).time];
    ecg_peaks = [ecg_peaks, ...
        logical(trials_lfp(t).ECG_spikes)];
    ecg_R2Rt = [ecg_R2Rt, ...
        trials_lfp(t).ECG_b2btime(trials_lfp(t).ECG_spikes & ...
        ~isnan(trials_lfp(t).ECG_b2btime))];
    
end

mean_ecg_R2Rt = mean(ecg_R2Rt);

ts = trials_lfp(1).tsample;
timestamps = linspace(0, ts*(length(LFP_raw) - 1), length(LFP_raw));

if distplot
    h = figure;
    histogram(ecg_R2Rt);

end
    
Rpeak_evoked_LFP = struct();
for i = 1:nshuffles
    fprintf('Shuffle %g\n', i);
    shuffled_ecg_p2pt = ecg_R2Rt(randperm(length(ecg_R2Rt)));
    peak_idx = 1;
    shuffled_ecg_peaks = false(size(ecg_peaks));
    npeaks = numel(shuffled_ecg_p2pt);
    shuffled_ecg_peak_times = cumsum(shuffled_ecg_p2pt...
        (peak_idx:peak_idx + npeaks - 1));
    shuffled_ecg_peak_times = shuffled_ecg_peak_times(...
        shuffled_ecg_peak_times < timestamps(end));
    shuffled_ECG_peak_samples = round(shuffled_ecg_peak_times / ts);
    % randomly shift shuffled ECG peaks
    shuffled_ECG_peak_samples = ...
        shuffled_ECG_peak_samples + ...
        round((2*rand(size(shuffled_ECG_peak_samples)) - 1) * (mean_ecg_R2Rt/ts));  
    shuffled_ECG_peak_samples(shuffled_ECG_peak_samples > ...
        length(ecg_peaks) | shuffled_ECG_peak_samples < 1) = [];
    % shuffled ECG peaks
    shuffled_ecg_peaks(shuffled_ECG_peak_samples) = 1;
    
    ft_data_all = struct();
    ft_data_all.trial = {[LFP_raw; shuffled_ecg_peaks]};
    ft_data_all.time = {timestamps};
    ft_data_all.label = {'LFP', 'Rpeak'};
    
    % now find the evoked LFP
    % evoked LFP average
    cfg                 = [];
    cfg.timwin          = [cfg_state{3} cfg_state{4}];
    cfg.channel         = ft_data_all.label(1); % LFP channel
    cfg.spikechannel    = ft_data_all.label(2); % ECG peak

    ecg_based_sta       = ft_spiketriggeredaverage(cfg, ft_data_all);
    
    Rpeak_evoked_LFP(i).lfp_time = ecg_based_sta.time;
    Rpeak_evoked_LFP(i).mean = ecg_based_sta.avg;    

end



shuffled_LFP_evoked.lfp = cat(1, Rpeak_evoked_LFP.mean);
shuffled_LFP_evoked.dimord = 'nshuffles_time';
shuffled_LFP_evoked.lfp_time = Rpeak_evoked_LFP(1).lfp_time;
shuffled_LFP_evoked.mean = nanmean(cat(1, Rpeak_evoked_LFP.mean), 1);
shuffled_LFP_evoked.std = nanstd(cat(1, Rpeak_evoked_LFP.mean), 0, 1);

end

