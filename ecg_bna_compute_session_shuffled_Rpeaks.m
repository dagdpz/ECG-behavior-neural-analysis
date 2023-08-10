function Rpeaks=ecg_bna_compute_session_shuffled_Rpeaks(session_info,ecg_bna_cfg)
load(session_info.Input_ECG);

Rpeaks = struct('block', {out.nrblock_combinedFiles});

offset_blocks_Rpeak=0;
for b=1:numel(out)
    Rpeaks(b).offset=offset_blocks_Rpeak(b);
    if isempty(out(b).nrblock_combinedFiles) || isempty(out(b).Rpeak_t) || isempty(out(b).R2R_t)
        Rpeaks(b).RPEAK_ts=[];
        Rpeaks(b).shuffled_ts=[];
        offset_blocks_Rpeak(b+1)=offset_blocks_Rpeak(b);
        continue
    end
    RPEAK_ts=[out(b).Rpeak_t(1) intersect(out(b).Rpeak_t,out(b).R2R_t)];
    RPEAKS_intervals=diff(RPEAK_ts);
    ecg_R2Rt_mean=mean(RPEAKS_intervals);
    idx_valid = RPEAKS_intervals<1.5*ecg_R2Rt_mean; %use mode or mean ?
    nonval_idx=find([0, ~idx_valid]);
    invalid_intervals(:,1)=RPEAK_ts(nonval_idx-1);
    invalid_intervals(:,2)=RPEAK_ts(nonval_idx);
    RPEAKS_intervals=RPEAKS_intervals(idx_valid);
    ecg_R2Rt_std=std(RPEAKS_intervals);
    N=ecg_bna_cfg.n_permutations;
    [~,ix]=sort(rand(N,length(RPEAKS_intervals)),2);
    perm=RPEAKS_intervals(ix);
    
    %% here we jitter all Rpeak intervals
    RPEAKS_intervals_p = [repmat(RPEAK_ts(1),N,1) repmat(RPEAKS_intervals,N,1)+randn(N,length(RPEAKS_intervals))*ecg_R2Rt_std perm]; % jittering every interval!
    RPEAK_ts_p=cumsum(RPEAKS_intervals_p,2);
    
    %% re-evaluating valid intervals... this is important to fix surrogates being higher due to periods of increased spiking that correlate with invalid Rpeaks
    for iv=1:size(invalid_intervals,1)
        RPEAK_ts_p(RPEAK_ts_p>invalid_intervals(iv,1)+ecg_R2Rt_mean/2 & RPEAK_ts_p<invalid_intervals(iv,2)-ecg_R2Rt_mean/2)=NaN;
    end
    RPEAK_ts_p(RPEAK_ts_p>max(RPEAK_ts)+max(RPEAKS_intervals))=NaN;
    RPEAK_ts_p(:,all(isnan(RPEAK_ts_p),1))=[];
    Rpeaks(b).RPEAK_ts=RPEAK_ts+offset_blocks_Rpeak(b);
    Rpeaks(b).shuffled_ts=RPEAK_ts_p+offset_blocks_Rpeak(b);
    offset_blocks_Rpeak(b+1)=offset_blocks_Rpeak(b)+max(RPEAK_ts)+max(RPEAKS_intervals)*2;
end

end