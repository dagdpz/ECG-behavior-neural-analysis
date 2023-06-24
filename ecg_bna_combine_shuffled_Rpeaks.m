function session_ecg = ecg_bna_combine_shuffled_Rpeaks(session_ecg, Rpeaks,ts,cfg)

n_shuffles=cfg.n_permutations;

Rpeak_blocks=[Rpeaks.block];
session_blocks=unique([session_ecg.trials.block]);
blocks=intersect(session_blocks,Rpeak_blocks);

%% fill up nonexisting blocks...
%     if isempty(ECG_timestamps)
%         fprintf('No ECG data found for block %g\n', b);
%         session_ecg.trials(trials_idx) = [];
%         return;
%     end
emptyblocks=session_blocks(~ismember(session_blocks,Rpeak_blocks));
for b=emptyblocks
    trials_idx = find([session_ecg.trials.block] == b);
    [session_ecg.trials(trials_idx).block]=deal(NaN);
end


for b=blocks
    block_Rpeak=Rpeaks(Rpeak_blocks==b);
    
    % find trials for this block
    trials_idx = find([session_ecg.trials.block] == b);
    
    % get ECG timestamps for this block
    ECG_timestamps = block_Rpeak.RPEAK_ts-block_Rpeak.offset;
    ECG_timestamps_shuffled = block_Rpeak.shuffled_ts(1:n_shuffles,:)-block_Rpeak.offset;
    %ECG_R2Rt = block_Rpeak.R2R_t;
    
    
    % ECG_R2Rsamples = block_Rpeak.R2R_sample;
    % ECG_R2Rvalid = block_Rpeak.R2R_valid;
    % ECG_R2Rvalid_bpm = block_Rpeak.R2R_valid_bpm;
    
    trials_time = vertcat(session_ecg.trials(trials_idx).trialperiod); %% uhmmmmmmmm doublecheck this
    if nargin<3
        ts = session_ecg.trials(trials_idx(1)).tsample;
    end
    session_ecg.tsample=ts;
    block_ecg_timestamps = (0:ts:round(trials_time(end)/ts)*ts);
    
    
    ECG_peaksamples = round(ECG_timestamps/ts) + 1;
    ECG_shuffled_peaksamples = round(ECG_timestamps_shuffled/ts) + 1;
    trials_samples = round(trials_time / ts) + 1;
    
    % ECG spikes based on ECG timestamps
    ECG_spikes = false(size(block_ecg_timestamps));
    ECG_spikes(ECG_peaksamples(ECG_peaksamples<numel(ECG_spikes))) = true;
    
    ECG_spikes_shuffled = false(size(ECG_shuffled_peaksamples,1),size(block_ecg_timestamps,2));
    for p=1:size(ECG_shuffled_peaksamples,1) %% this is rather stupid
        ECG_spikes_shuffled(p,ECG_shuffled_peaksamples(p,~isnan(ECG_shuffled_peaksamples(p,:)))) = true;
    end
    
    % ECG_b2bt = single(nan(size(block_ecg_timestamps)));
    % ECG_b2bt(ECG_R2Rsamples) = ECG_R2Rvalid;
    % ECG_bpm = single(nan(size(block_ecg_timestamps)));
    % ECG_bpm(ECG_R2Rsamples) = ECG_R2Rvalid_bpm;
    %
    % % fill missing values
    % nanx = isnan(ECG_b2bt);
    % t    = 1:numel(ECG_b2bt);
    % ECG_b2bt(nanx) = interp1(t(~nanx), ECG_b2bt(~nanx), t(nanx));
    % nanx = isnan(ECG_bpm);
    % t    = 1:numel(ECG_bpm);
    % ECG_bpm(nanx) = interp1(t(~nanx), ECG_bpm(~nanx), t(nanx));
    %
    % % remove invalid intervals
    % ECG_b2bt(~block_ecg_validsamples) = nan;
    % ECG_bpm(~block_ecg_validsamples) = nan;
    
    % now divide into trials
    for t = 1:length(trials_idx)
        
        trial_ECG_spikes    = ECG_spikes(trials_samples(t, 1):trials_samples(t, 2));
        trial_ECG_spikes_shuffled    = ECG_spikes_shuffled(:,trials_samples(t, 1):trials_samples(t, 2));
        session_ecg.trials(trials_idx(t)).ECG_spikes    = trial_ECG_spikes;
        session_ecg.trials(trials_idx(t)).ECG_spikes_shuffled    = trial_ECG_spikes_shuffled;
        session_ecg.trials(trials_idx(t)).nRpeaks       = sum(trial_ECG_spikes);
        
        %
        %     trial_ECG_bpm       = ECG_bpm(trials_samples(t, 1):trials_samples(t, 2));
        %     trial_ECG_b2btime   = ECG_b2bt(trials_samples(t, 1):trials_samples(t, 2));
        %     trial_ECG_valid     = block_ecg_validsamples(trials_samples(t, 1):trials_samples(t, 2));
        %     session_ecg.trials(trials_idx(t)).ECG_bpm       = trial_ECG_bpm;
        %     session_ecg.trials(trials_idx(t)).ECG_b2btime   = trial_ECG_b2btime;
        %     session_ecg.trials(trials_idx(t)).ECG_valid     = trial_ECG_valid;
    end
    
end
end