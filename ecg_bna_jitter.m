function Triggers=ecg_bna_jitter(IN,cfg)

N=cfg.n_permutations;

ts=IN.ts;
intervals=IN.intervals;
valid_idx=IN.valid_ts;
blocks=IN.blocks;
td=diff(ts); % These are the ACTUAL distances of our Rpeaks

%% shuffling the VALID R2R intervals to append at the end of our jittered intervals,
switch cfg.jitter_method
    case 'interval_jitter'
        %  just in case we randomly end up not covering the entire time window
        %  Don't worry, most (typically all) of it is going to be removed
        [~,ix]=sort(rand(N,length(intervals)),2);
        perm=intervals(ix);
        %% here we jitter all Rpeak intervals
        jitter_range = 2*std(intervals);                    % two-sided std of normally distributed jitter values
        td_jit   = [repmat(ts(1),N,1) repmat(td,N,1)+(randn(N,length(td))-0.5)*jitter_range perm]; % jittering every interval!
        ts_jit   = cumsum(td_jit,2);
        td_jit(1,:)=[];
        
        % now we define
        next_invalid=diff(valid_idx)~=1;                               % is the next Rpeak invalid (i.e. followed by invalid R2R interval)
        iv_starts  =[0  ts(valid_idx([next_invalid true]))];     % start of invalid intervals: Timestamps of valid Rpeaks followed by invalid ones
        % First Segment (for 0 to first valid Rpeak) and last segment
        % (everything after last valid Rpeak) are always invalid
        iv_ends    =[ts(valid_idx([true next_invalid]))   inf];  % end of invalid intervals: Timestamps of valid Rpeaks PRECEDED by invalid ones
        % we can get Rpeak-ts preceded by invalid R2R by shifting next_invalid
        grace_window=mean(intervals)/2;                         % +/- Range for shuffled Rpeaks to be allowed inside invalid segments
        
        
        %% remove jittered Rpeaks and corresponding durations that fell into invalid segments
        for iv=1:numel(iv_starts)
            idx2exclude_1 = ts_jit>iv_starts(iv)+grace_window & ts_jit<iv_ends(iv)-grace_window;
            ts_jit   (idx2exclude_1) = 0;
            td_jit (idx2exclude_1) = NaN;
        end
        idx2exclude_2 = ts_jit>max(ts)+jitter_range;
        ts_jit (idx2exclude_2) = 0;
        td_jit (idx2exclude_2) = NaN;
        
        ts_jit (:,all(ts_jit==0,1)) = [];
        td_jit (:,all(isnan(td_jit),1)) = [];
        
    case 'train_jitter'
        %% here we jitter the entire train of all Rpeaks instead
        jitter_range    = mean(intervals);                    %
        jitterthing     = repmat(rand(N,1)-0.5,1,numel(ts))*jitter_range;
        ts_jit          = repmat(ts,N,1)+jitterthing;    % jittered timestamps
        td_jit          = [zeros(N,1) diff(ts_jit,1,2)]; % durations of intervals (before corresponding trigger n)
    case 'trigger_jitter'
        %% here we jitter all Rpeaks individually instead - within a box that spans half into both intervals
        % ignore those with long intervals around them - not valid because not consecutive
                interval_min    = [0 td/2];
                interval_ran    = interval_min+[td/2 0];
%         interval_min    = repmat(mean(intervals)/2,size(ts));
%         interval_ran    = repmat(mean(intervals),size(ts));        
        
        jitterthing     = rand([N,numel(ts)]).*repmat(interval_ran,N,1)-repmat(interval_min,N,1);
        ts_jit          = repmat(ts,N,1)+jitterthing;
        td_jit          = [zeros(N,1) diff(ts_jit,1,2)];
        
end

%% take data corresponding to valid events
ts  = ts(valid_idx);                                          % take only Rpeaks surrounded by valid R2R
V   = repmat(ismember(1:size(ts_jit,2),valid_idx),N,1);       % logical index to reduce


%% put the data together
Triggers.blocks       =blocks;
Triggers.ts           =ts;
Triggers.intervals    =intervals; %% all intervals deemed valid (see outside function)
Triggers.shuffled_ts  =reshape(ts_jit(V),N,numel(ts));
%Triggers.shuffled_ts(Triggers.shuffled_ts==offset_blocks)=0;
Triggers.shuffled_intervals = reshape(td_jit(V),N,numel(ts)); % durations of reshuffled RR-intervals (the corresponding ends of those intervals are in Rpeaks(b).shuffled_ts)







%% This figure plots all R-peaks and consecutive R-peaks to make sure our selection procedure does the right thing
%     figure
%     hold on
%     for this_index_wont_be_used = 1:length(iv_starts)
%         f= fill([iv_starts(this_index_wont_be_used) iv_starts(this_index_wont_be_used) ...
%             iv_ends(this_index_wont_be_used) iv_ends(this_index_wont_be_used) ...
%             iv_starts(this_index_wont_be_used)], ...
%             [0 1 1 0 0], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
%     end
%     st1 = stem(out(b).Rpeak_t, ones(length(out(b).Rpeak_t),1)); % all the R-peaks from out.Rpeak_t
%     st2 = stem(out(b).R2R_t, 0.75*ones(length(out(b).R2R_t),1));
%     st3 = stem(RPEAK_ts, 0.5*ones(length(RPEAK_ts), 1),'g'); % plot consecutive R-peaks
%     legend([f st1 st2 st3], {'Invalid Time Window', 'All R-peaks', 'Ends of valid RRs', 'Consecutive R-peaks: 1 valid RR before & 1 valid RR after'})
%     xlabel('Time, s')
end

