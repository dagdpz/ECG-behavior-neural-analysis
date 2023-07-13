function triggered = ecg_bna_get_Rpeak_triggered_LFP( site_lfp, state )
% ecg_bna_get_Rpeak_evoked_LFP - computes the Rpeak evoked LFP for a specified
% time window around the Rpeak onset for given trials (usually trials
% belonging to a condition) in a session
%
% USAGE:
%	triggered_out = ecg_bna_get_Rpeak_evoked_LFP( trials_lfp, state )
%
% INPUTS:
%       trials          - 1xN struct containing LFP and Rpeak data of N
%       trials
%       state           - a cell array specifying time window around Rpeak
%       during which evoked response should be obtained
% OUTPUTS:
%		triggered_out  - struct containing Rpeak onset triggered
%		evoked LFP from the given trials
%
% See also ecg_bna_get_Rpeak_based_STA, ecg_bna_get_shuffled_Rpeak_evoked_LFP

state_name = state{2};
width = state{4} - state{3};
trials=site_lfp.trials;

% sample time
lfp_ts = 1/trials(1).fsample;
% number of samples in each window
w_nsamples = round(width/lfp_ts);

n_shuffles=size(trials(1).ECG_spikes,1);
freq_bp=size(trials(1).phase_bandpassed,2);

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

for sh=1:n_shuffles
    trigg_all.lfp       = {}; % evoked LFP response
    trigg_all.phaseBP   = {}; % bandpassed phase
    for t = 1:length(trials)
        trialperiod           = trials(t).trialperiod;
        
        % get the LFP samples and timestamps for the trial
        idx=trials(t).time >= trialperiod(1) & trials(t).time <= trialperiod(2);
        lfp_data = trials(t).lfp_data(idx);
        lfp_tfs_pBP = site_lfp.trials(t).phase_bandpassed;
        time = trials(t).time(idx);
        ecg_peaks = trials(t).ECG_spikes(:,idx);
        % now get the windows to combine
        w_center = find(ecg_peaks(sh,:));        
        % loop through each window
        for w = 1:length(w_center)
            if w_center(w) - round(w_nsamples/2) < 1 || w_center(w) + round(w_nsamples/2) > length(time)
                continue;
            end
            window=w_center(w) - round(w_nsamples/2):w_center(w) + round(w_nsamples/2);
            % evoked LFP for this state
            trigg_all.lfp = [trigg_all.lfp, lfp_data(window)];
            % bandpassed phase spectrum for this window
            trigg_all.phaseBP = [trigg_all.phaseBP, lfp_tfs_pBP(:,:,window)];
            % timestamps, set mid-timestamp to zero
            time =  time(window) ;
            time = time - time(round(length(time)/2));
        end
    end
    %triggered(sh).time   = time;
    cat_phaseBP = cat(1,trigg_all.phaseBP{:});
    % should we caclculate std over all time points in each freq (?)
    itpcbp.mean(sh,:,:) = abs(nanmean(exp(1i*cat_phaseBP), 1));
    itpcbp.std(sh,:,:) = repmat(nanstd(abs(mean(exp(1i*cat_phaseBP), 1)), 0, 3),1,1,length(time));
    
    % evoked LFP average
    cat_lfp = vertcat(trigg_all.lfp{:});
    lfp.mean(sh,:) = nanmean(cat_lfp, 1);
    lfp.std(sh,:) = nanstd(cat_lfp, 0, 1);
end

triggered.time=time;
triggered.state_name=state_name;
triggered.state=state;

%% differentiate between 
%%      a) mean of std and std of (shuffle) means
%%      a) mean confidence interval and confidence interval of (shuffle) means 

triggered.lfp.mean=mean(cat(1,lfp.mean),1);
triggered.lfp.std=std(cat(1,lfp.mean),1);
triggered.lfp.conf95 = prctile(cat(1,lfp.mean),[97.5, 2.5],1);
triggered.lfp.std_mean=mean(cat(1,lfp.std),1);
triggered.itpcbp.mean=mean(cat(1,itpcbp.mean),1);
triggered.itpcbp.std=std(cat(1,itpcbp.mean),1);
triggered.itpcbp.conf95 = prctile(cat(1,itpcbp.mean),[97.5, 2.5],1);
triggered.itpcbp.std_mean=mean(cat(1,itpcbp.std),1);
end

