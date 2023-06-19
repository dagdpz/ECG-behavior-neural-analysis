function ecg_triggered_tfs = ecg_bna_get_ECG_triggered_tfs_split_shuffled( site_lfp, cond_ecg, state, ecg_bna_cfg )
% ecg_bna_get_ECG_triggered_tfs_split_shuffled - gets the time-frequency
% pow,phs spectrogram and phsBP for a specified time window around all shuffled Rpeak onset for a single site for
% given trials (usually trials belonging to a condition) in a session
%
% USAGE:
%	ecg_triggered_tfs = ecg_bna_get_ECG_triggered_tfs_split_shuffled( site_lfp,
%	cond_ecg, state, ecg_bna_cfg )
%
% INPUTS:
%		site_lfp            - struct containing LFP data for all trials of
%		a session from a single site
%       cond_ecg         - ecg data of trials whose time freq spectrogram
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
% REQUIRES:	lfp_tfa_baseline_normalization ???
%
% See also ecg_bna_get_Rpeak_based_STA, ecg_bna_get_Rpeak_evoked_LFP,
% ecg_bna_get_shuffled_Rpeak_evoked_LFP
ecg_triggered_tfs = struct;

state_name = state{2};
width = state{4} - state{3};


n_shuffles = size(cond_ecg(1).ECG_spikes_shuffled,1);

[ecg_triggered_tfs(1:n_shuffles).powspctrm] = deal({}); % power spectrogram
[ecg_triggered_tfs(1:n_shuffles).phasespctrm] = deal({}); % phase spectrogram
[ecg_triggered_tfs(1:n_shuffles).phaseBP] = deal({}); % bandpassed phase

[ecg_triggered_tfs(1:n_shuffles).time] = deal({}); % timebins fo spectrogram
[ecg_triggered_tfs(1:n_shuffles).freq] = deal({}); % freq bins
[ecg_triggered_tfs(1:n_shuffles).state] = deal(NaN);%[];
[ecg_triggered_tfs(1:n_shuffles).state_name] = deal(state_name);

for t = 1:numel(site_lfp.trials)
    trialperiod           = site_lfp.trials(t).trialperiod;
    
    % get the LFP samples and timestamps between start and end states
    lfp_tfs_pow = site_lfp.trials(t).tfs.powspctrm;
    lfp_tfs_phs = site_lfp.trials(t).tfs.phase;
    lfp_tfs_pBP = site_lfp.trials(t).tfs.phase_bandpassed;
    lfp_tfs_time = site_lfp.trials(t).tfs.time;
    lfp_tfs_freq = site_lfp.trials(t).tfs.freq;
    
    % ecg peak times
    lfp_time = site_lfp.trials(t).time;
    
    for sh = 1: n_shuffles
        
        ecg_peaks = cond_ecg(t).ECG_spikes_shuffled(sh,:);
        ecg_peak_times = lfp_time(ecg_peaks(1:min(numel(ecg_peaks),numel(lfp_time))));
        
        % lfp sample indices
        % lfp_timebin_idx = 1:length(lfp_tfs_time);
        % width of each time bin
        lfp_ts = 1/site_lfp.trials(t).fsample;
        tbin_width = ecg_bna_cfg.tfr.timestep * lfp_ts;
        
        % number for time bins in the window
        ntbins_window = round(width/tbin_width);
        
        %ecg_triggered_tfs.cfg = site_lfp.trials(t).tfs.cfg;
        
        % now get the midpoint of windows around ECG peak
        window_mid_idx = [];
        for p = ecg_peak_times
            if min(abs(lfp_tfs_time - p)) > tbin_width
                continue;
            end
            window_mid_idx = [window_mid_idx, find(abs(lfp_tfs_time - p) ==  min(abs(lfp_tfs_time - p)))];
        end
        
        % loop through each window
        for w = 1:length(window_mid_idx)
            if window_mid_idx(w) - round(ntbins_window/2) < 1 || window_mid_idx(w) + round(ntbins_window/2) > length(lfp_tfs_time)
                continue;
            end
            window=window_mid_idx(w) - round(ntbins_window/2):window_mid_idx(w) + round(ntbins_window/2);
            % power spectrum for this window
            ecg_triggered_tfs(sh).powspctrm = [ecg_triggered_tfs(sh).powspctrm, lfp_tfs_pow(:,:,window)];
            % phase spectrum for this window
            ecg_triggered_tfs(sh).phasespctrm = [ecg_triggered_tfs(sh).phasespctrm, lfp_tfs_phs(:,:,window)];
            % bandpassed phase spectrum for this window
            ecg_triggered_tfs(sh).phaseBP = [ecg_triggered_tfs(sh).phaseBP, lfp_tfs_pBP(:,:,window)];
            
            
            % timestamps
            ecg_triggered_tfs(sh).time = lfp_tfs_time(window);
            % set mid-timestamp to zero
            ecg_triggered_tfs(sh).time = ecg_triggered_tfs(sh).time - ecg_triggered_tfs(sh).time(round(length(ecg_triggered_tfs(sh).time)/2));
            ecg_triggered_tfs(sh).freq = lfp_tfs_freq;
        end
    end
end


for sh = 1: size(cond_ecg(t).ECG_spikes_shuffled , 1)
    % find number of time bins in power spectrogram
    ntimebins = min(cellfun('size', ecg_triggered_tfs(sh).powspctrm, 3));
    % crop each tfs to the ntimebins
    for k = 1:length(ecg_triggered_tfs(sh).powspctrm)
        ecg_triggered_tfs(sh).powspctrm{k} = ecg_triggered_tfs(sh).powspctrm{k}(1,:,1:ntimebins);
        ecg_triggered_tfs(sh).phasespctrm{k} = ecg_triggered_tfs(sh).phasespctrm{k}(1,:,1:ntimebins);
    end
    ecg_triggered_tfs(sh).time = ecg_triggered_tfs(sh).time(1:ntimebins);
    
    % average power spectrum for each state
    % arr_state_pow = zeros(1, nfreqbins, ntimebins);
    
    if ~isempty(ecg_triggered_tfs(sh).powspctrm)
        % first clear out the parameters to avoid overwrites:
        arr_state_pow       = [];
        arr_state_phase     = [];
        arr_state_phaseBP   = [];
        
        % find the average TFS for each state in each shuffled data
        arr_state_pow = cat(1, ecg_triggered_tfs(sh).powspctrm{:});
        arr_state_phase = cat(1,ecg_triggered_tfs(sh).phasespctrm{:});
        arr_state_phaseBP = cat(1,ecg_triggered_tfs(sh).phaseBP{:});
        ecg_triggered_tfs(sh).powspctrm_rawmean = nanmean(arr_state_pow, 1);
        ecg_triggered_tfs(sh).powspctrm_rawstd  = nanstd(arr_state_pow, 1);
        
        
        ecg_triggered_tfs(sh).phasespctrm_rawmean   = abs(nanmean(exp(1i*arr_state_phase), 1));
        ecg_triggered_tfs(sh).phasespctrm_rawstd    = nanstd(abs(mean(exp(1i*arr_state_phase), 1)), 0, 3); 
        % should we caclculate std over all time points in each freq (?)
        ecg_triggered_tfs(sh).phaseBP_rawmean = abs(nanmean(exp(1i*arr_state_phaseBP), 1));
        ecg_triggered_tfs(sh).phaseBP_rawstd  = nanstd(abs(mean(exp(1i*arr_state_phaseBP), 1)), 0, 3);
        % should we caclculate std over all time points in each freq (?)
        
        % baseline normalization
        cfg_baseline.method = ecg_bna_cfg.baseline_method;
        
        %% BASELINE
        if ~strcmp(cfg_baseline.method, 'none')
            
            baseline_cnd_idx=...
                find_conditon(site_lfp,ecg_bna_cfg,'perturbation','baseline_perturbation') & ...
                find_conditon(site_lfp,ecg_bna_cfg,'choice','baseline_use_choice_trial') & ...
                find_conditon(site_lfp,ecg_bna_cfg,'type','baseline_use_type') & ...
                find_conditon(site_lfp,ecg_bna_cfg,'effector','baseline_use_effector');
            
            %% We can not pre-calculate baseline (for ITPC!), we need to calculate it AFTER triggering, so somewhere here.
            cfg_baseline.mean = nanmean(vertcat(site_lfp.baseline(baseline_cnd_idx).pow_mean),1);
            cfg_baseline.std  = nanmean(vertcat(site_lfp.baseline(baseline_cnd_idx).pow_std),1);
            ecg_triggered_tfs(sh).powspctrm_normmean = lfp_tfa_baseline_normalization(ecg_triggered_tfs(sh).powspctrm_rawmean, cfg_baseline);
        else
            ecg_triggered_tfs(sh).powspctrm_normmean = ecg_triggered_tfs(sh).powspctrm_rawmean;
        end
    end
end

end

function baseline_cnd_idx=find_conditon(site_lfp,ecg_bna_cfg,site_field,cfg_field)
if isinf(ecg_bna_cfg.(cfg_field))
    baseline_cnd_idx = isinf([site_lfp.baseline.(site_field)]);
elseif isnan(ecg_bna_cfg.(cfg_field))
    baseline_cnd_idx = isnan([site_lfp.baseline.(site_field)]);
else
    baseline_cnd_idx = ismember([site_lfp.baseline.(site_field)],ecg_bna_cfg.(cfg_field));
end
end

