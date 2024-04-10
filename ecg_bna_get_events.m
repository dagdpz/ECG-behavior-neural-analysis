function Event=ecg_bna_get_events(session_info,trials,event,block_anchors)

%ecg_bna_get_triggers(session_info
eventtype=strrep(event{1,2}, ' ','');
%                     cfg.event_types=cfg.analyse_states(:,2);
%                     cfg.events=cfg.analyse_states(:,1);
%                         eventtype=strrep(cfg.event_types{e}, ' ','');
%% load respective file
switch eventtype
    case {'CAP'}
    load(session_info.Input_ECG);
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
    
    case {'Rpeak','Rpeak_insp','Rpeak_exp'} %% why load it several times ... not that it matters right now
    load(session_info.Input_ECG);
    
end


        Event.valid_ts     = [];
        Event.ts           = [];
        Event.intervals    = [];
        Event.blocks       = [];
switch eventtype
    case {'Rpeak','Rpeak_insp','Rpeak_exp','CAP'}
        %% treat block_merging as invalid?
        
        for b=1:numel(out)
            
            %% adjustment for isnp/exp
            valid_interval=out(b).idx_valid_R2R_consec;
            if strcmp(eventtype,'Rpeak_insp')
                valid_interval=intersect(find(out(b).is_RR_insp),valid_interval);
            elseif strcmp(eventtype,'Rpeak_exp')
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
            offset=block_anchors(block); %%%%
            
            Event.valid_ts        = [Event.valid_ts valid_idx+numel(Event.ts) ]; % +numel(ts)
            Event.ts              = [Event.ts out(b).Rpeak_t+offset ]; %+offset
            Event.intervals       = [Event.intervals intervals];
            Event.blocks          = [Event.blocks block];
%             B=B+1;
%             valid_offsets(B)      = offset; % how/do we need to signify which block it is ?
%             valid_blocks(B)       = out(b).nrblock_combinedFiles;
        end
    case 'microstim'
    case 'state'
        state=event{3};
        %find state in each trial
        [O]=ecg_bna_synchronize_with_trigger(block_anchors,trials);
        Event.ts=O.state_onsets(O.states==state);
        Event.valid_ts     = 1:numel(Event.ts);
        Event.blocks       = unique([trials.block]);
        Event.intervals    = O.trial_ends-O.trial_starts;
        
end