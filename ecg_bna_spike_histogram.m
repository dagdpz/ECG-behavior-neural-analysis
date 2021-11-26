function ecg_bna_spike_histogram(basepath_ecg,basepath_spikes,basepath_to_save,session_info)

%20210706

session_info{1}={'Bacchus','20211001',[1 2 3 4 5 6 7]};
session_info{2}={'Bacchus','20210720',[4 5 6 7 8]};

basepath_ecg='Y:\Projects\Pulv_distractor_spatial_choice\Data\';
basepath_spikes='Y:\Projects\Pulv_distractor_spatial_choice\ephys\StimulusType_Difficulty_Position_LS\';
basepath_to_save='Y:\Projects\Pulv_distractor_spatial_choice\Data\ECG_triggered_PSTH';


Sanity_check=0;


ECG_event=-1;
keys.PSTH_WINDOWS={'ECG',ECG_event,-0.5,0.5};
keys.PSTH_binwidth=0.01;
keys.kernel_type='gaussian';
keys.gaussian_kernel=0.02;



if ~exist(basepath_to_save,'dir')
    mkdir(basepath_ecg, 'ECG_triggered_PSTH');
end


condition=struct('unit',{});

u=0; % unit counter across sessions

for s=1:numel(session_info)
    monkey=session_info{s}{1};
    session=session_info{s}{2};
    blocks=session_info{s}{3};
    
    
    load([basepath_ecg monkey filesep 'ECG' filesep session filesep session  '_ecg.mat']);
    load([basepath_spikes 'population_' monkey '_' session '.mat']);
    
    
    for U=1:numel(population)
        pop=population(U);
        unit_ID=population(U).unit_ID;
        u=U+u;
        
        for b=1:numel(blocks)
            block=blocks(b);
            
            o=find([out.nrblock_combinedFiles]==block); %% this would be the easy version, but somehow this number can also be empty...
            for oo=1:numel(out)
                if out(oo).nrblock_combinedFiles==block
                    o=oo;
                end
            end
            
            %% here we could potentially further reduce trials
            tr=[pop.trial.block]==block;
            
            
            trcell=num2cell(pop.trial(tr));
            % add trial onset time to each spike so its basically one stream again
            % also, make sure spikes aren't counted twice (because previous
            % trial is appended in beginning;
            arrival_times=cellfun(@(x) [x.arrival_times(x.arrival_times>0)]+x.TDT_ECG1_t0_from_rec_start+x.TDT_ECG1_tStart,trcell,'Uniformoutput',false);
            trial_onsets=cellfun(@(x) x.TDT_ECG1_t0_from_rec_start+x.TDT_ECG1_tStart,trcell); % add trial onset time to each spike so its basically one stream again
            trial_ends=cellfun(@(x) x.states_onset(x.states==98)+x.TDT_ECG1_t0_from_rec_start+x.TDT_ECG1_tStart,trcell,'Uniformoutput',false); % no clue why this needs to be nonuniformoutput, it did work earlier so this is confusing...
            trial_ends=[trial_ends{:}];
            
            AT=vertcat(arrival_times{:});
            if ~numel(trial_onsets)>0
%                 AT=AT-trial_onsets(1);                      %% subtract trial onset time of first trial per block so all arrival times are aligned to fixation acquisition of first trial!
%                 trial_ends=trial_ends-trial_onsets(1);      %% subtract trial onset time of first trial per block so all arrival times are aligned to fixation acquisition of first trial!
%                 trial_onsets=trial_onsets-trial_onsets(1);  %% subtract trial onset time of first trial per block so all arrival times are aligned to fixation acquisition of first trial!
%             else
                continue;
            end
            RPEAK_ts=[out(o).Rpeak_t];
            
            %% reduce RPEAK_ts potentially ? (f.e.: longer than recorded ephys, inter-trial-interval?)
            during_trial_index = arrayfun(@(x) any(trial_onsets<=x+keys.PSTH_WINDOWS{1,3} & trial_ends>=x+keys.PSTH_WINDOWS{1,4}),RPEAK_ts);
            RPEAK_ts=RPEAK_ts(during_trial_index);
            
            %% find if to append it to rest or task condition
            %not ideal, but should work for  now
            if trcell{1}.type==1
                tasktype=1; % rest
            elseif trcell{1}.type==2
                tasktype=2; % task
            else
                continue
            end
            
            %% check how many Rpeaks ("trials") we have already for that condition, so we can append across blocks
            if numel(condition)>=tasktype && numel(condition(tasktype).unit)>=u  && isfield (condition(tasktype).unit(u),'trial')
                T=numel(condition(tasktype).unit(u).trial);
            else
                T=0;
            end
            
            %% now the tricky part: sort by ECG peaks ...
            for t=1:numel(RPEAK_ts)
                condition(tasktype).unit(u).trial(t+T).states=ECG_event;
                condition(tasktype).unit(u).trial(t+T).states_onset=0;
                AT_temp=AT-RPEAK_ts(t);
                condition(tasktype).unit(u).trial(t+T).arrival_times=AT_temp(AT_temp>keys.PSTH_WINDOWS{1,3}-0.2 & AT_temp<keys.PSTH_WINDOWS{1,4}+0.2);
            end
            
            %% The part following here is internal sanity check and should be turned off in general since there typically is no ECG data in the spike format
            if Sanity_check
                ECG_data=[pop.trial(tr).TDT_ECG1];
                ECG_time=cellfun(@(x) [1/x.TDT_ECG1_SR:1/x.TDT_ECG1_SR:numel(x.TDT_ECG1)/x.TDT_ECG1_SR]+x.TDT_ECG1_t0_from_rec_start+x.TDT_ECG1_tStart,trcell,'Uniformoutput',false);
                CG_time=[ECG_time{:}];
                tt=0;
                clear ECG_to_plot
                for t=1:numel(RPEAK_ts)
                    ECg_t_idx=ECG_time>RPEAK_ts(t)-0.5 & ECG_time<RPEAK_ts(t)+0.5;
                    if round(sum(ECg_t_idx)-trcell{1}.TDT_ECG1_SR)~=0
                        continue
                    end
                    tt=tt+1;
                    ECG_to_plot(tt,:)=ECG_data(ECg_t_idx);
                end
                
                figure
                filename=[unit_ID '_block_' num2str(block) '_ECG_average'];
                lineProps={'color','b','linewidth',1};
                shadedErrorBar(1:size(ECG_to_plot,2),mean(ECG_to_plot,1),sterr(ECG_to_plot,1),lineProps,1);
                export_fig([basepath_to_save, filesep, filename], '-pdf','-transparent') % pdf by run
                close(gcf);
                
                figure
                filename=[unit_ID '_block_' num2str(block) '_ECG_per5trials'];
                plot(ECG_to_plot(1:5:end,:)');
                export_fig([basepath_to_save, filesep, filename], '-pdf','-transparent') % pdf by run
                close(gcf);
            end
            
        end
        
        figure;
        title(unit_ID,'interpreter','none');
        hold on
        if numel(condition(1).unit) >= u
            trial=condition(1).unit(u).trial;
            [SD  bins SD_VAR SD_SEM]=ph_spike_density(trial,1,keys,zeros(size(trial)),ones(size(trial)));
            shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,SD,SD_SEM,lineProps,1);
        end
        
        if numel(condition(2).unit) >= u
            trial=condition(2).unit(u).trial;
            [SD  bins SD_VAR SD_SEM]=ph_spike_density(trial,1,keys,zeros(size(trial)),ones(size(trial)));
            lineProps={'color','r','linewidth',1};
            shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,SD,SD_SEM,lineProps,1);
        end
        
        filename=unit_ID;
        export_fig([basepath_to_save, filesep, filename], '-pdf','-transparent') % pdf by run
        close(gcf);
        
        
    end
    
    
    
end

%% Here comes some sort of across population plot i assume?
aaa=12;


end





