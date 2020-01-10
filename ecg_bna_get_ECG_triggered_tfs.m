function ecg_triggered_tfs = ecg_bna_get_ECG_triggered_tfs( site_lfp, cond_trials, state, ecg_bna_cfg )
% ecg_bna_get_ECG_triggered_tfs - gets the LFP time frequency spectrogram
% for a specified time window around the Rpeak onset for a single site for
% given trials (usually trials belonging to a condition) in a session
%
% USAGE:
%	ecg_triggered_tfs = ecg_bna_get_ECG_triggered_tfs( site_lfp,
%	cond_trials, state, ecg_bna_cfg ) 
%
% INPUTS:
%		site_lfp            - struct containing LFP data for all trials of 
%		a session from a single site
%       cond_trials         - indices of trials whose time freq spectrogram
%       are to be obtaied
%       state               - a cell array specifying time window during
%       which LFP tfs should be obtained
%       ecg_bna_cfg         - struct containing settings
%       Required fields:
%           tfr.timestep            - width of a timebin in LFP TFS in
%           number of LFP samples
%           baseline_method         - method used for baseline
%           normalization of LFP TFS
%           baseline_perturbation   - whether to use
%           control(0)/inactivation(1) trials for baseline
%           baseline_use_choice_trial - whether to use
%           choice(1)/instructed(0) trials for baseline
%           baseline_use_type       - type number mof trials to be used for
%           baseline 
%           baseline_use_effector   - effector value of trials to be used
%           for baseline 
% OUTPUTS:
%		ecg_triggered_tfs   - struct containing LFP tfs for all given
%		trials
% 
% REQUIRES:	lfp_tfa_baseline_normalization
%
% See also ecg_bna_get_Rpeak_based_STA, ecg_bna_get_Rpeak_evoked_LFP, 
% ecg_bna_get_shuffled_Rpeak_evoked_LFP

state_name = state{2};
width = state{4} - state{3};

ecg_triggered_tfs.powspctrm = {}; % power spectrogram
ecg_triggered_tfs.time = {}; % timebins fo spectrogram
ecg_triggered_tfs.freq = {}; % freq bins
ecg_triggered_tfs.state = [];
ecg_triggered_tfs.state_name = state_name;

for t = find(cond_trials)

    trialperiod           = site_lfp.trials(t).trialperiod;
    
    % get the LFP samples and timestamps between start and end states
    lfp_tfs_pow = site_lfp.trials(t).tfs.powspctrm;
    lfp_tfs_time = site_lfp.trials(t).tfs.time;
    lfp_tfs_freq = site_lfp.trials(t).tfs.freq;
    
    % ecg peak times
    lfp_time = site_lfp.trials(t).time; 
    ecg_peaks = site_lfp.trials(t).ECG_spikes; %...
%         ((site_lfp.trials(t).time >= trialperiod(1) & ...
%         site_lfp.trials(t).time <= trialperiod(2))); 
    ecg_peak_times = lfp_time(ecg_peaks);
    
    % lfp sample indices
    lfp_timebin_idx = 1:length(lfp_tfs_time);
    % width of each time bin
    lfp_ts = 1/site_lfp.trials(t).fsample;
    tbin_width = ecg_bna_cfg.tfr.timestep * lfp_ts;
    
    % number fo time bins in the window
    ntbins_window = round(width/tbin_width);
    
    ecg_triggered_tfs.cfg = site_lfp.trials(t).tfs.cfg;
    
    % now get the midpoint of windows around ECG peak
    window_mid_idx = [];
    for p = ecg_peak_times
        if min(abs(lfp_tfs_time - p)) > tbin_width
            continue;
        end
        window_mid_idx = [window_mid_idx, find(abs(lfp_tfs_time - p) == ...
            min(abs(lfp_tfs_time - p)))];
    end

    % loop through each window
    for w = 1:length(window_mid_idx)
        if window_mid_idx(w) - round(ntbins_window/2) < 1 || ...
                window_mid_idx(w) + round(ntbins_window/2) > length(lfp_tfs_time)
            continue;
        end
        % power spectrum for this window
        ecg_triggered_tfs.powspctrm = [ecg_triggered_tfs.powspctrm, ...
            lfp_tfs_pow(:, :, window_mid_idx(w) - round(ntbins_window/2):...
             window_mid_idx(w) + round(ntbins_window/2))];
        % timestamps
        ecg_triggered_tfs.time = ...
            lfp_tfs_time(window_mid_idx(w) - round(ntbins_window/2):...
             window_mid_idx(w) + round(ntbins_window/2));
        % set mid-timestamp to zero
        ecg_triggered_tfs.time = ecg_triggered_tfs.time - ...
            ecg_triggered_tfs.time(round(length(ecg_triggered_tfs.time)/2));

        ecg_triggered_tfs.freq = lfp_tfs_freq;
    end
       
end


% find number of time bins in power
% spectrogram
ntimebins = min(cellfun('size', ecg_triggered_tfs.powspctrm, 3));
nfreqbins = numel(ecg_triggered_tfs.freq);
% crop each tfs to the ntimebins
for k = 1:length(ecg_triggered_tfs.powspctrm)
    ecg_triggered_tfs.powspctrm{k} = ecg_triggered_tfs.powspctrm{k}(1,:,1:ntimebins);
    
end
ecg_triggered_tfs.time = ecg_triggered_tfs.time(1:ntimebins);

% average power spectrum for each state
arr_state_pow = zeros(1, nfreqbins, ntimebins);

if ~isempty(ecg_triggered_tfs.powspctrm)

    % find the average TFS for each state
    arr_state_pow = cat(1, ecg_triggered_tfs.powspctrm{:});
    ecg_triggered_tfs.powspctrm_rawmean = nanmean(arr_state_pow, 1);

    % baseline normalization
    cfg_baseline.method = ecg_bna_cfg.baseline_method;
    if ~strcmp(cfg_baseline.method, 'none')
        if isinf(ecg_bna_cfg.baseline_perturbation)
            baseline_cnd_idx = isinf([site_lfp.baseline.perturbation]);
        end
        if isinf(ecg_bna_cfg.baseline_use_choice_trial)
            baseline_cnd_idx = baseline_cnd_idx & isinf([site_lfp.baseline.choice]);
        end
        if isinf(ecg_bna_cfg.baseline_use_type)
            baseline_cnd_idx = baseline_cnd_idx & isinf([site_lfp.baseline.type]);
        end
        if isinf(ecg_bna_cfg.baseline_use_effector)
            baseline_cnd_idx = baseline_cnd_idx & isinf([site_lfp.baseline.effector]);
        end
        if ~isinf(ecg_bna_cfg.baseline_perturbation) && ~isinf(ecg_bna_cfg.baseline_use_choice_trial)
            baseline_cnd_idx = ...
                [site_lfp.baseline.perturbation] == ecg_bna_cfg.baseline_perturbation & ...
                [site_lfp.baseline.choice] == ecg_bna_cfg.baseline_use_choice_trial & ...
                [site_lfp.baseline.type] == ecg_bna_cfg.baseline_use_type & ...
                [site_lfp.baseline.effector] == ecg_bna_cfg.baseline_use_effector;
        end
        cfg_baseline.mean = site_lfp.baseline(baseline_cnd_idx).pow_mean;
        cfg_baseline.std = site_lfp.baseline(baseline_cnd_idx).pow_std;
        ecg_triggered_tfs.powspctrm_normmean = lfp_tfa_baseline_normalization(...
            ecg_triggered_tfs.powspctrm_rawmean, cfg_baseline); 
    else
        ecg_triggered_tfs.powspctrm_normmean = ecg_triggered_tfs.powspctrm_rawmean;
    end
    
    
end

end

