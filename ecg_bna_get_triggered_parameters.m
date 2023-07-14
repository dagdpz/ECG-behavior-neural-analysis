function triggered = ecg_bna_get_triggered_parameters( site_lfp, state, ecg_bna_cfg )
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
w_nsamples = round(width/lfp_ts);
freq_bp=size(trials(1).phase_bandpassed,2);
tbin_width = ecg_bna_cfg.tfr.timestep * lfp_ts;
w_nsamples_tfr = round(width/tbin_width);
freq = site_lfp.trials(1).tfs.freq;

nanpreallocator=NaN(n_shuffles, numel(freq), round(w_nsamples_tfr/2)*2+1);

pow.mean    = nanpreallocator;
pow.std     = nanpreallocator;
pow.conf95_high     = nanpreallocator;
pow.conf95_low     = nanpreallocator;
itpc.mean   = nanpreallocator;
itpc.std    = nanpreallocator;
itpc.conf95_high     = nanpreallocator;
itpc.conf95_low     = nanpreallocator;


nanpreallocator=NaN(n_shuffles, freq_bp, round(w_nsamples/2)*2+1);
itpcbp.mean    = nanpreallocator;
itpcbp.std     = nanpreallocator;
itpcbp.conf95_high     = nanpreallocator;
itpcbp.conf95_low     = nanpreallocator;

nanpreallocator=NaN(n_shuffles, round(w_nsamples/2)*2+1);
lfp.mean   = nanpreallocator;
lfp.std    = nanpreallocator;
lfp.conf95_high     = nanpreallocator;
lfp.conf95_low     = nanpreallocator;

for sh = 1: n_shuffles
    trigg.pow = {}; % power spectrogram
    trigg.phase = {}; % phase spectrogram
    trigg.lfp       = {}; % evoked LFP response
    trigg.phaseBP   = {}; % bandpassed phase
    for t = 1:numel(trials)
        % get the LFP samples and timestamps between start and end states
        trial_pow       = trials(t).tfs.powspctrm;
        trial_phs       = trials(t).tfs.phase; 
        trial_lfp       = trials(t).lfp_data;
        trial_pBP       = trials(t).phase_bandpassed;
        trial_tfs_time  = trials(t).tfs.time;
        trial_lfp_time  = trials(t).time;     
        
        % ecg peak times
        ecg_peaks  = trials(t).ECG_spikes(sh,:);
        ecg_peak_times = trial_lfp_time(ecg_peaks(1:min(numel(ecg_peaks),numel(trial_lfp_time))));
        w_idx=find(ecg_peaks);
        
        % now get the midpoint of windows around ECG peak
        w_tfr_mid = [];
        w_lfp_mid = [];
        for idx=1:numel(ecg_peak_times)
            p = ecg_peak_times(idx);
            if min(abs(trial_tfs_time - p)) > tbin_width
                continue;
            end
            c_w_mid_tfr=find(abs(trial_tfs_time - p) ==  min(abs(trial_tfs_time - p)));
            c_w_mid_lfp=w_idx(idx);
            %if w_center(w) - round(w_nsamples/2) < 1 || 
            if c_w_mid_lfp - round(w_nsamples/2) < 1 ||...
               c_w_mid_lfp + round(w_nsamples/2) > length(trial_lfp_time) || ...
               c_w_mid_tfr - round(w_nsamples_tfr/2) < 1 ||...
               c_w_mid_tfr + round(w_nsamples_tfr/2) > min(length(trial_tfs_time),size(trial_pow,3))
                continue;
            end
            w_tfr_mid = [w_tfr_mid, c_w_mid_tfr];
            w_lfp_mid = [w_lfp_mid, c_w_mid_lfp];
        end
        
        % loop through each window
        for w = 1:length(w_tfr_mid)
            w_tfr=w_tfr_mid(w) - round(w_nsamples_tfr/2):w_tfr_mid(w) + round(w_nsamples_tfr/2);
            % power spectrum for this window
            trigg.pow = [trigg.pow, trial_pow(:,:,w_tfr)];
            % phase spectrum for this window
            trigg.phase = [trigg.phase, trial_phs(:,:,w_tfr)];
            % timestamps
            tfr_time = trial_tfs_time(w_tfr);
            % set mid-timestamp to zero
            tfr_time = tfr_time - tfr_time(round(length(tfr_time)/2));
            
            w_lfp=w_lfp_mid(w) - round(w_nsamples/2):w_lfp_mid(w) + round(w_nsamples/2);
            % evoked LFP for this state
            trigg.lfp = [trigg.lfp, trial_lfp(:,w_lfp)];
            % bandpassed phase spectrum for this window
            trigg.phaseBP = [trigg.phaseBP, trial_pBP(:,:,w_lfp)];
            % timestamps, set mid-timestamp to zero
            lfp_time =  trial_lfp_time(w_lfp) ;
            lfp_time = lfp_time - lfp_time(round(length(lfp_time)/2));
        end
    end
    
    if ~isempty(trigg.pow)
        % find the average TFS for each state in each shuffled data
        cat_pow = cat(1, trigg.pow{:});
        cat_phase = cat(1,trigg.phase{:});
        pow.mean(sh,:,:)    = nanmean(cat_pow, 1);
        pow.std(sh,:,:)     = nanstd(cat_pow, 1);
        pow.conf95_high(sh,:,:)  = prctile(cat_pow,97.5,1);
        pow.conf95_low(sh,:,:)  = prctile(cat_pow,2.5,1);
        itpc.mean(sh,:,:)   = abs(nanmean(exp(1i*cat_phase), 1));
        
        %% these still make no sense, think about using phase information for std
        % problem a) wrapping for circular phase to be taken as closest to the average pahse
        % problem b) scaling of std and confidence interval
        
        itpc.std(sh,:,:)    = repmat(nanstd(abs(mean(exp(1i*cat_phase), 1)), 0, 3),1,1,length(tfr_time));
        itpc.conf95_high(sh,:,:)  = prctile(cat_phase,97.5,1);
        itpc.conf95_low(sh,:,:)   = prctile(cat_phase,2.5,1);
        
        cat_phaseBP = cat(1,trigg.phaseBP{:});
        % should we caclculate std over all time points in each freq (?)
        itpcbp.mean(sh,:,:) = abs(nanmean(exp(1i*cat_phaseBP), 1));
        itpcbp.std(sh,:,:) = repmat(nanstd(abs(mean(exp(1i*cat_phaseBP), 1)), 0, 3),1,1,length(lfp_time));
        
        % evoked LFP average
        cat_lfp = vertcat(trigg.lfp{:});
        lfp.mean(sh,:) = nanmean(cat_lfp, 1);
        lfp.std(sh,:) = nanstd(cat_lfp, 0, 1);
    end
end
triggered.time=lfp_time;
triggered.freq=freq;
triggered.state_name=state_name;
triggered.state=state;
triggered.tfr_time=tfr_time;

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

triggered.lfp.mean=mean(cat(1,lfp.mean),1);
triggered.lfp.std=std(cat(1,lfp.mean),1);
triggered.lfp.conf95 = prctile(cat(1,lfp.mean),[97.5, 2.5],1);
triggered.lfp.std_mean=mean(cat(1,lfp.std),1);

triggered.itpcbp.mean=mean(cat(1,itpcbp.mean),1);
triggered.itpcbp.std=std(cat(1,itpcbp.mean),1);
triggered.itpcbp.conf95 = prctile(cat(1,itpcbp.mean),[97.5, 2.5],1);
triggered.itpcbp.std_mean=mean(cat(1,itpcbp.std),1);
end

