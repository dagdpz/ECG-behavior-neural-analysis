function bpm_Change_all_Blocks_LS

ecg_preprocess_folder='Y:\Data\BodySignals\ECG_CAP\';
monkeys={'Magnus','Bacchus'};

%% these are the ones we want to run - they have both task and rest condition in at least one of the nuclei
clear sessions Task Rest
sessions{1}=sort([20221115, 20221118, 20221122, 20221206, ...
    20221222, 20221229, 20230104, 20230106, 20230112, ...
    20230126, 20230511, 20230518, 20230519, 20230524, ...
    20230525, 20230526, 20230601, 20230602, 20230607, ...
    20230609, 20230614, 20230615, 20230616, 20230621, 20230623]);

sessions{2} = unique([20210716, 20210720, 20210723, ...
    20210730, 20210805, 20210826, 20210827, 20210905, ...
    20210906, 20210930, 20211007, 20211012, 20211013, ...
    20211014, 20211019, 20211028, 20211103, 20211116, ...
    20211117, 20211207, 20211214, 20220105, 20220106, ...
    20220203, 20220211, 20220221, 20220222, 20220224, ...
    20220225, 20220309, 20220310, 20220315, 20220318, 20220322]);

%%
for m = 1: length(monkeys)
    %% check data directory
    Monkey = monkeys{m};
    monkeydir=[ecg_preprocess_folder filesep Monkey filesep];
    all_ECG_out_data = dir([monkeydir filesep '*.mat']);
    for s=1:numel(all_ECG_out_data)
        all_ECG_out_data(s).session=str2double(all_ECG_out_data(s).name(1:8));
    end
    
    %% loop through desired sessions
    sessions_to_loop=sessions{1,m};    
    sn = 0;
    stoskip=[];
    for s = 1:numel(sessions_to_loop)
        session=sessions_to_loop(s);
        si=[all_ECG_out_data.session]==session;
        if ~any(si)
            stoskip= [stoskip, s];
            continue;
        else        
            load([monkeydir all_ECG_out_data(si).name],'out')
            load(['Y:\Projects\Pulv_bodysignal\ephys\ECG_TaskRest_',Monkey,...
                '_merged\trials_',Monkey,'_',num2str(session),'.mat']) 
            %fprintf('\n Number of Block in session %s of Monkey %s = %i \n',num2str(session),Monkey,length(out));
            
            type = [];
            btokeep=[];
            R2R_valid=[];
            for b = 1:numel(out)
                if  isempty(out(b).nrblock_combinedFiles) || isnan(out(b).nrblock_combinedFiles)
                    continue;
                end
                b_idx = find([trials.block] == out(b).nrblock_combinedFiles);
                trial_nBlocks_type = unique([trials(b_idx).type]);
                if isempty(b_idx) || isempty(out(b).Rpeak_t) || isempty(out(b).R2R_t) || isempty(out(b).idx_valid_R2R_consec)
                    continue
                end
                type{b} = num2str(trial_nBlocks_type);
                R2R_valid{b} = out(b).R2R_valid;
                btokeep= [btokeep, b];
            end
            type=type(btokeep);
            out=out(btokeep);
            R2R_valid=R2R_valid(btokeep);
            if  (numel(type)<1) || (numel(unique(type)) < 2) % skip session if not both task types are present
                stoskip= [stoskip, s];
                continue;
            end
            
            type = strrep(type,'1','Rest');
            type = strrep(type,'2','Task');
            
        end
        if ~ismember(s,stoskip)
            sn=sn+1;
            r=ismember(type,'Rest');
            t=ismember(type,'Task');
            BPMtask=60./[R2R_valid{t}];
            BPMrest=60./[R2R_valid{r}];
            Task.(Monkey).mean_bpm(sn)   = mean(BPMtask);
            Task.(Monkey).median_bpm(sn) = median(BPMtask);
            Task.(Monkey).std_bpm(sn)    = median(BPMtask);
            Rest.(Monkey).mean_bpm(sn)   = mean(BPMrest);
            Rest.(Monkey).median_bpm(sn) = median(BPMrest);
            Rest.(Monkey).std_bpm(sn)    = median(BPMrest);
        end
    end
end

%% store
results_file = fullfile(['Y:\Projects\Pulv_bodysignal\HeartRate_Changes_all_Blocks',filesep,'Bpm_both_monkeys_taskrest']);
save(results_file,'Task','Rest');
end
