function triggered = ecg_bna_get_ECG_triggered_tfr( site_lfp, cond_ecg, state, ecg_bna_cfg )
% ecg_bna_get_triggered_split_shuffled - gets the time-frequency
% pow,phs spectrogram and phsBP for a specified time window around all shuffled Rpeak onset for a single site for
% given trials (usually trials belonging to a condition) in a session
%
% USAGE:
%	triggered = ecg_bna_get_triggered_split_shuffled( site_lfp,
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
%		triggered   - struct containing LFP tfs for all given
%		trials
%
% REQUIRES:	lfp_tfa_baseline_normalization ???
%
% See also ecg_bna_get_Rpeak_based_STA, ecg_bna_get_Rpeak_evoked_LFP,
% ecg_bna_get_shuffled_Rpeak_evoked_LFP
triggered = struct;

state_name = state{2};
width = state{4} - state{3};


trials_lfp=site_lfp.trials;
n_shuffles=size(trials_lfp(1).ECG_spikes,1);

% [triggered(1:n_shuffles).powspctrm] = deal({}); % power spectrogram
% [triggered(1:n_shuffles).phasespctrm] = deal({}); % phase spectrogram
[spctrm(1:n_shuffles).pow] = deal({}); % power spectrogram
[spctrm(1:n_shuffles).phase] = deal({}); % phase spectrogram

[triggered(1:n_shuffles).time] = deal({}); % timebins fo spectrogram
[triggered(1:n_shuffles).freq] = deal({}); % freq bins
[triggered(1:n_shuffles).state] = deal(NaN);%[];
[triggered(1:n_shuffles).state_name] = deal(state_name);

for t = 1:numel(site_lfp.trials)
    trialperiod           = site_lfp.trials(t).trialperiod;
    
    % get the LFP samples and timestamps between start and end states
    lfp_tfs_pow  = site_lfp.trials(t).tfs.powspctrm;
    lfp_tfs_phs  = site_lfp.trials(t).tfs.phase;
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
        w_nsamples = round(width/tbin_width);
        
        %triggered.cfg = site_lfp.trials(t).tfs.cfg;
        
        % now get the midpoint of windows around ECG peak
        w_mid = [];
        for p = ecg_peak_times
            if min(abs(lfp_tfs_time - p)) > tbin_width
                continue;
            end
            w_mid = [w_mid, find(abs(lfp_tfs_time - p) ==  min(abs(lfp_tfs_time - p)))];
        end
        
        % loop through each window
        for w = 1:length(w_mid)
            if w_mid(w) - round(w_nsamples/2) < 1 || w_mid(w) + round(w_nsamples/2) > min(length(lfp_tfs_time),size(lfp_tfs_pow,3))
                continue;
            end
            window=w_mid(w) - round(w_nsamples/2):w_mid(w) + round(w_nsamples/2);
            % power spectrum for this window
            spctrm(sh).pow = [spctrm(sh).pow, lfp_tfs_pow(:,:,window)];
            % phase spectrum for this window
            spctrm(sh).phase = [spctrm(sh).phase, lfp_tfs_phs(:,:,window)];
            
            
            % timestamps
            triggered(sh).time = lfp_tfs_time(window);
            % set mid-timestamp to zero
            triggered(sh).time = triggered(sh).time - triggered(sh).time(round(length(triggered(sh).time)/2));
            triggered(sh).freq = lfp_tfs_freq;
        end
    end
end


for sh = 1: n_shuffles
    % find number of time bins in power spectrogram
    ntimebins = min(cellfun('size', spctrm(sh).pow, 3));
    % crop each tfs to the ntimebins
    for k = 1:length(spctrm(sh).pow)
        spctrm(sh).pow{k} = spctrm(sh).pow{k}(1,:,1:ntimebins);
        spctrm(sh).phase{k} = spctrm(sh).phase{k}(1,:,1:ntimebins);
    end
    triggered(sh).time = triggered(sh).time(1:ntimebins);
    
    % average power spectrum for each state
    % arr_state_pow = zeros(1, nfreqbins, ntimebins);
    
    if ~isempty(spctrm(sh).pow)
        
        % find the average TFS for each state in each shuffled data
        cat_pow = cat(1, spctrm(sh).pow{:});
        cat_phase = cat(1,spctrm(sh).phase{:});
        triggered(sh).pow.mean = nanmean(cat_pow, 1);
        triggered(sh).pow.std  = nanstd(cat_pow, 1);
        
        
        triggered(sh).itpc.mean   = abs(nanmean(exp(1i*cat_phase), 1));
        triggered(sh).itpc.std    = nanstd(abs(mean(exp(1i*cat_phase), 1)), 0, 3); 
        % should we caclculate std over all time points in each freq (?)
        
%         % baseline normalization
%         cfg_baseline.method = ecg_bna_cfg.baseline_method;
%         
%         %% BASELINE
%         if ~strcmp(cfg_baseline.method, 'none')
%             
%             baseline_cnd_idx=...
%                 find_conditon(site_lfp,ecg_bna_cfg,'perturbation','baseline_perturbation') & ...
%                 find_conditon(site_lfp,ecg_bna_cfg,'choice','baseline_use_choice_trial') & ...
%                 find_conditon(site_lfp,ecg_bna_cfg,'type','baseline_use_type') & ...
%                 find_conditon(site_lfp,ecg_bna_cfg,'effector','baseline_use_effector');
%             
%             %% We can not pre-calculate baseline (for ITPC!), we need to calculate it AFTER triggering, so somewhere here.
%             cfg_baseline.mean = nanmean(vertcat(site_lfp.baseline(baseline_cnd_idx).pow_mean),1);
%             cfg_baseline.std  = nanmean(vertcat(site_lfp.baseline(baseline_cnd_idx).pow_std),1);
%             triggered(sh).powspctrm_normmean = lfp_tfa_baseline_normalization(triggered(sh).powspctrm_rawmean, cfg_baseline);
%         else
%             triggered(sh).powspctrm_normmean = triggered(sh).powspctrm_rawmean;
%         end
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

