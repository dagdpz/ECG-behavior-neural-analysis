function out=ecg_bna_compute_PSTH(Trigger,SD,PSTH_time,during_trial_index,bins,cfg)

% RPEAK_ts,RPEAK_dur,RAST,
compute_PSTH(RPEAK_ts,RPEAK_dur,RAST,SD,PSTH_time,during_trial_index,bins,cfg)
RPEAK_samples=round((RPEAK_ts-PSTH_time(1))/cfg.spk.PSTH_binwidth);

%% preallocate "out" structure
n=size(RPEAK_samples,1);
m=numel(bins);

out.raster    = NaN;
out.mean      = nan([n m]);
out.std       = nan([n m]);
out.SEM       = nan([n m]);
out.n_events  = NaN([n 1]);
out.RTs       = {};
out.RDs       = {};

%% loop through rows of RPEAK_samples: 1 row for real, nReshuffles rows of reshuffled data
for p=1:n
% %% remove samples that would land outside
    rpeaks2_valid = RPEAK_samples(p,:)<=numel(SD) & RPEAK_samples(p,:)>0;
    
    RT=RPEAK_ts(p,rpeaks2_valid);      % take time-stamps
    RS=RPEAK_samples(p,rpeaks2_valid); % take samples
    RD=RPEAK_dur(p,rpeaks2_valid);     % take durations
    
    within_trial_idx = during_trial_index(RS);
    
    RT=RT(within_trial_idx);
    RS=RS(within_trial_idx);
    RD=RD(within_trial_idx);
    
    idx_by_event = bsxfun(@plus,RS',bins); % bsxfun produces a RPEAK-(nBins-1) matrix with samples taken from SDF for each R-peak
    if n == 1
        out.raster = RAST(idx_by_event);
    end
    PSTH = SD(idx_by_event); 
    out.mean(p,:)=mean(PSTH, 1);
    out.std(p,:) = std(PSTH, [], 1);
    out.SEM(p,:)=sterr(PSTH);
    out.n_events(p)=numel(RS);
    out.RTs{p}=RT;
    out.RDs{p}=RD;
end
end
