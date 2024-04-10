function Event=ecg_bna_get_events(session_info,trial,event,block_anchors)

%ecg_bna_get_triggers(sessions_info(i)

if strcmp(event,'CAP') %% adjust output file
    load(sessions_info(i).Input_ECG);
    nrblockscell={out.nrblock_combinedFiles};
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
    [out.nrblock_combinedFiles] = deal(nrblockscell{:});
    
elseif strcmp(event,'Rpeak')
    load(sessions_info(i).Input_ECG);
    
end

switch event
    case {'Rpeak','Rpeak_insp','Rpeak_exp','CAP'}
        %% treat block_merging as invalid?
        Event.valid_ts     = [];
        Event.ts           = [];
        Event.intervals    = [];
        
        offset=0;
        B=0;
        for b=1:numel(out)
            
            %% adjustment for isnp/exp
            valid_interval=out(b).idx_valid_R2R_consec;
            if strcmp(event,'Rpeak_insp')
                valid_interval=intersect(find(out(b).is_RR_insp),valid_interval);
            elseif strcmp(event,'Rpeak_exp')
                valid_interval=intersect(find(out(b).is_RR_exp),valid_interval);
            end
            
            %% is this needed ?
            %             if idx_valid_R2R_consec(1) == 1
            %                 % drop the first Rpeak if it's consecutive as the jittering
            %                 % function complains about such a situation
            %                 idx_valid_R2R_consec(1) = [];
            %             end
            
            [~,consecutive_idx]=ismember(out(b).R2R_t(valid_interval),out(b).Rpeak_t); % indexing of consecutive_idx corrsponds to out(b).R2R_t
            valid_idx=consecutive_idx-1;
            intervals    = out(b).R2R_valid(valid_interval);
            
            if isempty(valid_idx) %isempty(block) || isempty(out(b).Rpeak_t) || isempty(out(b).R2R_t) || isempty(out(b).idx_valid_R2R_consec)
                continue
            end
            
            
            %% add to offset
            
            %% use trial information to get INI length
            %% for a given block, this is what it is (?) -> do we have a global offset (?)
            %trials
            block=out(b).nrblock_combinedFiles;
            %offset=offset+max(max(Event.ts))+max(intervals);%onset of recording relative to state 2
            offset=offset+blockstart(block); %%%%
            
            Event.valid_ts        = [Event.valid_ts valid_idx+numel(Event.ts) ]; % +numel(ts)
            Event.ts              = [Event.ts out(b).Rpeak_t+offset ]; %+offset
            Event.intervals       = [Event.intervals intervals];
            
            
            B=B+1;
            valid_offsets(B)      = offset; % how/do we need to signify which block it is ?
            valid_blocks(B)       = out(b).nrblock_combinedFiles;
        end
    case 'microstim'
    case 'state'
        %create one continous string of state_onsets,
        %trials
        
end