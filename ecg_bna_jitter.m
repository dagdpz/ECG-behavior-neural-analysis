function Rpeaks=ecg_bna_jitter(session_info,cfg)
load(session_info.Input_ECG);

% get the number of permutations from settings
N=cfg.n_permutations;

% % condition for using CAP data instead of ECG 
% if (ecg_bna_cfg.outNameCap)  
%     oldnames = fieldnames(out_cap);
%     newnames = strrep(oldnames,'B2B','R2R');
%     idx = find(~strcmpi(newnames,oldnames));
%     for n = 1: numel(newnames)
%         [out_cap.(newnames{n})] = deal(out_cap.(oldnames{n}));
%     end
%     for n = 1:numel(idx)
%         out_cap = rmfield(out_cap,oldnames{idx(n)});
%     end
%     out = out_cap;
%     nrblockscell=num2cell(ses.nrblock_combinedFiles);
%     [out.nrblock_combinedFiles] = deal(nrblockscell{:});
% end

offset_blocks_Rpeak=0;
for b=1:numel(out)
    Rpeaks(b).block=out(b).nrblock_combinedFiles;
    Rpeaks(b).offset=offset_blocks_Rpeak(b);
    if isempty(out(b).nrblock_combinedFiles) || isempty(out(b).Rpeak_t) || isempty(out(b).R2R_t)
        Rpeaks(b).block=NaN;
        Rpeaks(b).RPEAK_ts=[];
        Rpeaks(b).RR_durations = [];
        Rpeaks(b).shuffled_ts=[];
        offset_blocks_Rpeak(b+1)=offset_blocks_Rpeak(b);
        continue
    end
    
    RPEAK_ts=out(b).Rpeak_t; 
        
	%% here we jitter all Rpeak intervals
	%allowed_jitter_range=min(max(out(b).R2R_valid)-min(out(b).R2R_valid),min(out(b).R2R_valid));                    % one-sided std of normally distributed jitter values
    allowed_jitter_range=max(max(out(b).R2R_valid)-min(out(b).R2R_valid),min(out(b).R2R_valid));                    % one-sided std of normally distributed jitter values

    Nrandoms=numel(RPEAK_ts);
    jitterthing=repmat(rand(N,1)-0.5,1,Nrandoms)*allowed_jitter_range;
    RPEAK_ts_p         = repmat(RPEAK_ts,N,1)+jitterthing(:,1:numel(RPEAK_ts));
    RPEAK_ts_dur       = [zeros(N,1) diff(RPEAK_ts_p,1,2)];
    
    %% figure out consecutive RR-intervals - we do it only after jittering because we need to know durations of jittered intervals corresponding to the consecutive ones
    [~,consecutive_idx]=ismember(out(b).R2R_t(out(b).idx_valid_R2R_consec),out(b).Rpeak_t); % indexing of consecutive_idx corrsponds to out(b).Rpeak_t
    valid_idx=consecutive_idx-1;
    
    %% take data corresponding to consecutive R-peaks
    RPEAK_ts     = RPEAK_ts(valid_idx);                                  % take only Rpeaks surrounded by valid R2R
    RPEAK_dur    = out(b).R2R_valid(out(b).idx_valid_R2R_consec-1);      % 
    
    idx2exclude_2 = RPEAK_ts_p>max(RPEAK_ts)+allowed_jitter_range | RPEAK_ts_p<0;
    RPEAK_ts_p    (idx2exclude_2) = 0;
    RPEAK_ts_dur  (idx2exclude_2) = NaN;
    
    RPEAK_ts_p   (:,all(RPEAK_ts_p==0,1)) = [];
    RPEAK_ts_dur (:,all(isnan(RPEAK_ts_dur),1)) = [];
    
    Rpeaks(b).RPEAK_ts          =RPEAK_ts+offset_blocks_Rpeak(b);         % this offset is just a trick to be able to append Rpeaks across blocks easily
    Rpeaks(b).shuffled_ts       =RPEAK_ts_p+offset_blocks_Rpeak(b);
    Rpeaks(b).RPEAK_ts(Rpeaks(b).RPEAK_ts==offset_blocks_Rpeak(b))=0;    
    Rpeaks(b).shuffled_ts(Rpeaks(b).shuffled_ts==offset_blocks_Rpeak(b))=0;
    
    Rpeaks(b).RPEAK_dur=RPEAK_dur; % durations of RR-intervals (the corresponding ends of those intervals are in Rpeaks(b).RPEAK_ts)
    Rpeaks(b).shuffled_dur = RPEAK_ts_dur; % durations of reshuffled RR-intervals (the corresponding ends of those intervals are in Rpeaks(b).shuffled_ts)
    offset_blocks_Rpeak(b+1)=offset_blocks_Rpeak(b)+max(RPEAK_ts)+allowed_jitter_range*2;
end

end