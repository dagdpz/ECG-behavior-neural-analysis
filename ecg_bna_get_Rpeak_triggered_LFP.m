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

n_shuffles=size(trials(1).ECG_spikes,1);
[trigg_all(1:n_shuffles).time]      = deal({}); % timestamps
[trigg_all(1:n_shuffles).lfp]       = deal({}); % evoked LFP response
[trigg_all(1:n_shuffles).phaseBP]   = deal({}); % bandpassed phase
%[triggered(1:n_shuffles).state]     = deal({}); % not sure if this is needed (?)
[triggered(1:n_shuffles).state_name] = deal(state_name); % not sure if this is needed (?)
isvalid=0;
for t = 1:length(trials)
    trialperiod           = trials(t).trialperiod;
    
    % get the LFP samples and timestamps for the trial
    idx=trials(t).time >= trialperiod(1) & trials(t).time <= trialperiod(2);
    lfp_data = trials(t).lfp_data(idx);
    lfp_tfs_pBP = site_lfp.trials(t).phase_bandpassed;
    time = trials(t).time(idx);
    ecg_peaks = trials(t).ECG_spikes(:,idx);
    
    % sample time
    lfp_ts = 1/trials(t).fsample;
    % number of samples in each window
    w_nsamples = round(width/lfp_ts);
    
    for sh=1:n_shuffles
        % now get the windows to combine
        w_center = find(ecg_peaks(sh,:));
        
        % loop through each window
        for w = 1:length(w_center)
            if w_center(w) - round(w_nsamples/2) < 1 || w_center(w) + round(w_nsamples/2) > length(time)
                continue;
            end
            window=w_center(w) - round(w_nsamples/2):w_center(w) + round(w_nsamples/2);
            % evoked LFP for this state
            trigg_all(sh).lfp = [trigg_all(sh).lfp, lfp_data(window)];
            % bandpassed phase spectrum for this window
            trigg_all(sh).phaseBP = [trigg_all(sh).phaseBP, lfp_tfs_pBP(:,:,window)];
            % timestamps, set mid-timestamp to zero
            temp_time =  time(window) ;
            temp_time = temp_time - temp_time(round(length(temp_time)/2));
            trigg_all(sh).time     = [trigg_all(sh).time, temp_time];
            isvalid=1;
        end
    end
end

if isvalid
    % crop each lfp to same number of samples
    nsamples = min(cellfun('length', trigg_all(1).lfp));%% question here really is what happens if some shuffles do not contain an RPeak in one of hte trials
   
    for sh=1:n_shuffles
        for k = 1:length(trigg_all(sh).lfp)
            trigg_all(sh).lfp{k} = trigg_all(sh).lfp{k}(1:nsamples);
            triggered(sh).time   = trigg_all(sh).time{k}(1:nsamples);
        end
        %triggered(sh).time = trigg_all(sh).time(1:nsamples);
        cat_phaseBP = cat(1,trigg_all(sh).phaseBP{:});
        % should we caclculate std over all time points in each freq (?)
        triggered(sh).itpcbp.mean = abs(nanmean(exp(1i*cat_phaseBP), 1));
        triggered(sh).itpcbp.std  = nanstd(abs(mean(exp(1i*cat_phaseBP), 1)), 0, 3);
        
        % evoked LFP average
        cat_lfp = vertcat(trigg_all(sh).lfp{:});
        triggered(sh).lfp.mean = nanmean(cat_lfp, 1);
        triggered(sh).lfp.std = nanstd(cat_lfp, 0, 1);
    end
else
    triggered.time = [];
    triggered.lfp.mean = [];
    triggered.lfp.std = [];
    triggered.itpcbp.mean = [];
    triggered.itpcbp.std = [];
end
end

