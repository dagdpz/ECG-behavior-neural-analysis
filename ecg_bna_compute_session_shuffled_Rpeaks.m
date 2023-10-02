function Rpeaks=ecg_bna_compute_session_shuffled_Rpeaks(session_info,ecg_bna_cfg)
load(session_info.Input_ECG);

% condition for using CAP data instead of ECG 
if (ecg_bna_cfg.outNameCap)  
    oldnames = fieldnames(out_cap);
    newnames = strrep(oldnames,'B2B','R2R');
    idx = find(~strcmpi(newnames,oldnames));
    for n = 1: numel(newnames)
        [out_cap.(newnames{n})] = deal(out_cap.(oldnames{n}));
    end
    for n = 1:numel(idx)
        out_cap = rmfield(out_cap,oldnames{idx(n)});
    end
    out = out_cap;
    nrblockscell=num2cell(ses.nrblock_combinedFiles);
    [out.nrblock_combinedFiles] = deal(nrblockscell{:});
end

offset_blocks_Rpeak=0;
for b=1:numel(out)
    Rpeaks(b).block=out(b).nrblock_combinedFiles;
    Rpeaks(b).offset=offset_blocks_Rpeak(b);
    if isempty(out(b).nrblock_combinedFiles) || isempty(out(b).Rpeak_t) || isempty(out(b).R2R_t)
        Rpeaks(b).block=NaN;
        Rpeaks(b).RPEAK_ts=[];
        Rpeaks(b).shuffled_ts=[];
        offset_blocks_Rpeak(b+1)=offset_blocks_Rpeak(b);
        continue
    end
    
    RPEAK_ts=out(b).Rpeak_t;    
    [~,consecutive_idx]=ismember(out(b).R2R_t(out(b).idx_valid_R2R_consec),out(b).Rpeak_t);
    valid_idx=consecutive_idx-1;
    %valid_idx=out(b).idx_valid_R2R_consec-1;                       % idx of R_peaks surrounded by valid R2R intervals
    next_invalid=diff(valid_idx)~=1;                               % is the next Rpeak invalid (i.e. followed by invalid R2R interval)
    iv_starts  =[0  RPEAK_ts(valid_idx([next_invalid true]))];     % start of invalid intervals: Timestamps of valid Rpeaks followed by invalid ones
                                                                   % First Segment (for 0 to first valid Rpeak) and last segment
                                                                   % (everything after last valid Rpeak) are always invalid
    iv_ends    =[RPEAK_ts(valid_idx([true next_invalid]))   inf];  % end of invalid intervals: Timestamps of valid Rpeaks PRECEDED by invalid ones  
                                                                   % we can get Rpeak-ts preceded by invalid R2R by shifting next_invalid  
    grace_window=mean(out(b).R2R_valid)/2;                         % +/- Range for shuffled Rpeaks to be allowed inside invalid segments
    allowed_jitter_range=std(out(b).R2R_valid);                    % one-sided std of of normally distributed jitter values
    RPEAK_ts=RPEAK_ts(out(b).idx_valid_R2R_consec);                % take only Rpeaks surrounded by valid R2R
    RPEAKS_intervals=diff(RPEAK_ts);                               % These are the intervals of our valid Rpeaks, the ones we are going to jitter
    
    %% shuffling the VALID R2R intervals to append at the end of our jittered intervals, 
    %  just in case we randomly end up not covering the entire time window
    %  Don't worry, most (typically all) of it is going to be removed 
    N=ecg_bna_cfg.n_permutations;
    [~,ix]=sort(rand(N,length(out(b).R2R_valid)),2);
    perm=out(b).R2R_valid(ix);
    
    %% here we jitter all Rpeak intervals
    RPEAKS_intervals_p = [repmat(RPEAK_ts(1),N,1) repmat(RPEAKS_intervals,N,1)+randn(N,length(RPEAKS_intervals))*allowed_jitter_range perm]; % jittering every interval!
    RPEAK_ts_p         = cumsum(RPEAKS_intervals_p,2);
    
    %% remove jittered Rpeaks that fell into invalid segments
    for iv=1:numel(iv_starts)
        RPEAK_ts_p(RPEAK_ts_p>iv_starts(iv)+grace_window & RPEAK_ts_p<iv_ends(iv)-grace_window)=NaN;
    end
    RPEAK_ts_p(RPEAK_ts_p>max(RPEAK_ts)+allowed_jitter_range)=NaN;
    RPEAK_ts_p(:,all(isnan(RPEAK_ts_p),1))=[];
    Rpeaks(b).RPEAK_ts=RPEAK_ts+offset_blocks_Rpeak(b);         % this offset is just a trick to be able to append Rpeaks across blocks easily
    Rpeaks(b).shuffled_ts=RPEAK_ts_p+offset_blocks_Rpeak(b);
    offset_blocks_Rpeak(b+1)=offset_blocks_Rpeak(b)+max(RPEAK_ts)+allowed_jitter_range*2;
end

end