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

trials=site_lfp.trials;
n_shuffles=size(trials(1).ECG_spikes,1);
lfp_ts = 1/site_lfp.trials(1).fsample;
tbin_width = ecg_bna_cfg.tfr.timestep * lfp_ts;

% number for time bins in the window
w_nsamples = round(width/tbin_width);
freq = site_lfp.trials(1).tfs.freq;

nanpreallocator=NaN(n_shuffles, numel(freq), round(w_nsamples/2)*2+1);

pow.mean    = nanpreallocator;
pow.std     = nanpreallocator;
pow.conf95_high     = nanpreallocator;
pow.conf95_low     = nanpreallocator;
itpc.mean   = nanpreallocator;
itpc.std    = nanpreallocator;
itpc.conf95_high     = nanpreallocator;
itpc.conf95_low     = nanpreallocator;

for sh = 1: n_shuffles
    [spctrm.pow] = {}; % power spectrogram
    [spctrm.phase] = {}; % phase spectrogram
    for t = 1:numel(trials)
        % get the LFP samples and timestamps between start and end states
        lfp_tfs_pow  = trials(t).tfs.powspctrm;
        lfp_tfs_phs  = trials(t).tfs.phase;
        lfp_tfs_time = trials(t).tfs.time;
        %lfp_tfs_freq = trials(t).tfs.freq;
        
        % ecg peak times
        trial_time = trials(t).time;
        ecg_peaks = cond_ecg(t).ECG_spikes_shuffled(sh,:);
        ecg_peak_times = trial_time(ecg_peaks(1:min(numel(ecg_peaks),numel(trial_time))));
        
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
            spctrm.pow = [spctrm.pow, lfp_tfs_pow(:,:,window)];
            % phase spectrum for this window
            spctrm.phase = [spctrm.phase, lfp_tfs_phs(:,:,window)];
            % timestamps
            time = lfp_tfs_time(window);
            % set mid-timestamp to zero
            time = time - time(round(length(time)/2));
        end
    end
    
    if ~isempty(spctrm.pow)
        % find the average TFS for each state in each shuffled data
        cat_pow = cat(1, spctrm.pow{:});
        cat_phase = cat(1,spctrm.phase{:});
        pow.mean(sh,:,:)    = nanmean(cat_pow, 1);
        pow.std(sh,:,:)     = nanstd(cat_pow, 1);
        pow.conf95_high(sh,:,:)  = prctile(cat_pow,97.5,1);
        pow.conf95_low(sh,:,:)  = prctile(cat_pow,2.5,1);
        itpc.mean(sh,:,:)   = abs(nanmean(exp(1i*cat_phase), 1));
        
        %% these still make no sense, think about using phase information for std
        % problem a) wrapping for circular phase to be taken as closest to the average pahse
        % problem b) scaling of std and confidence interval
        
        itpc.std(sh,:,:)    = repmat(nanstd(abs(mean(exp(1i*cat_phase), 1)), 0, 3),1,1,length(time));
        itpc.conf95_high(sh,:,:)  = prctile(cat_phase,97.5,1);
        itpc.conf95_low(sh,:,:)   = prctile(cat_phase,2.5,1);
    end
end
triggered.time=time;
triggered.freq=freq;
triggered.state_name=state_name;
triggered.state=state;

%% differentiate between
%%      a) mean of std and std of (shuffle) means
%%      a) mean confidence interval and confidence interval of (shuffle) means
triggered.pow.mean=mean(cat(1,pow.mean),1);
triggered.pow.std=std(cat(1,pow.mean),1);
triggered.pow.std_mean=mean(cat(1,pow.std),1);
triggered.pow.conf95=prctile(cat(1,pow.mean),[97.5, 2.5],1);

triggered.itpc.std=std(cat(1,itpc.mean),1);
triggered.itpc.conf95=prctile(cat(1,itpc.mean),[97.5, 2.5],1);
triggered.itpc.mean=mean(cat(1,itpc.mean),1);
triggered.itpc.std_mean=mean(cat(1,itpc.std),1);
end

% function baseline_cnd_idx=find_conditon(site_lfp,ecg_bna_cfg,site_field,cfg_field)
% if isinf(ecg_bna_cfg.(cfg_field))
%     baseline_cnd_idx = isinf([site_lfp.baseline.(site_field)]);
% elseif isnan(ecg_bna_cfg.(cfg_field))
%     baseline_cnd_idx = isnan([site_lfp.baseline.(site_field)]);
% else
%     baseline_cnd_idx = ismember([site_lfp.baseline.(site_field)],ecg_bna_cfg.(cfg_field));
% end
% end

