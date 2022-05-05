function Output=ecg_bna_compute_session_spike_histogram(session_info)
savePlot = 1;
Sanity_check=0; % ECG triggered ECG, turn off since typically there is no ECG data in the spike format
remove_ini=1;   % to remove inter-trial intervals from ECG peaks (useful if ITI spikes were excluded during waveclus preprocessing)
keys.n_permutations=100; %100;
    keys.significance_window=[-0.25 0.25];

ECG_event=-1;
keys.PSTH_WINDOWS={'ECG',ECG_event,-0.5,0.5};
keys.PSTH_binwidth=0.01;
keys.kernel_type='gaussian';
keys.gaussian_kernel=0.02;

basepath_to_save=[session_info.SPK_fldr filesep 'per_unit'];
if ~exist(basepath_to_save,'dir')
    mkdir(basepath_to_save);
end

condition(1).unit=[];
condition(2).unit=[];

% Rpeaks derived from concatenated ECG data [First_trial_INI.ECG1 trial.TDT_ECG1]
load(session_info.Input_ECG);
load(session_info.Input_spikes);
Output = [];
    
for u=1:numel(population)
    pop=population(u);
    unit_ID=population(u).unit_ID;
    target =population(u).target;
    blocks=unique([pop.trial.block]);
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
        % also, make sure spikes aren't counted twice (because previous trial is appended in beginning;
        % removing overlapping spikes here            % add trial onset time
        arrival_times=cellfun(@(x) [x.arrival_times(x.arrival_times>x.states_onset(x.states==1))]+x.TDT_ECG1_t0_from_rec_start,trcell,'uniformoutput',false);
        % trial_onsets and ends only relevant if removing ITI
        trial_onsets=cellfun(@(x) x.TDT_ECG1_t0_from_rec_start+x.TDT_ECG1_tStart,trcell);
        trial_ends=cellfun(@(x) x.states_onset(x.states==98)+x.TDT_ECG1_t0_from_rec_start+x.TDT_ECG1_tStart,trcell,'uniformoutput',false); % no clue why this needs to be nonuniformoutput, it did work earlier so this is confusing...
        trial_ends=[trial_ends{:}];
        
        AT=vertcat(arrival_times{:});
        if ~numel(trial_onsets)>0 || isempty(o) || isempty(out(o).Rpeak_t(1)) || isempty(out(o).R2R_t)
            continue;
        end
        
        %% re-evaluating valid intervals... this is important to fix surrogates being higher due to periods of increased spiking that correlate with invalid Rpeaks
        RPEAK_ts=[out(o).Rpeak_t(1) intersect(out(o).Rpeak_t,out(o).R2R_t)];
        RPEAK_ts(RPEAK_ts> trial_ends(end))=[];
        RPEAKS_intervals=diff(RPEAK_ts);
        idx_valid = RPEAKS_intervals<1.5*mode(RPEAKS_intervals);
        invalid_intervals=[NaN,NaN];
        nonval_idx=find([0, ~idx_valid]);
        for iv=1:numel(nonval_idx)
            invalid_intervals(iv,1)=RPEAK_ts(nonval_idx(iv)-1);
            invalid_intervals(iv,2)=RPEAK_ts(nonval_idx(iv));
        end
        RPEAKS_intervals=RPEAKS_intervals(idx_valid);
        ecg_R2Rt_std=std(RPEAKS_intervals);
        ecg_R2Rt_mean=mean(RPEAKS_intervals);
        
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
        
        %% now the tricky part: sort by ECG peaks first, compute spike density across all trials!
        SD_full_block=ecg_bna_spike_density(AT,trial_onsets,trial_ends,keys);
        sorted=sort_by_rpeaks_temp(RPEAK_ts,SD_full_block,trial_onsets,trial_ends,keys,remove_ini);
        condition(tasktype).unit(u).trial(T+1:T+numel(sorted))=sorted;
        
        %% make surrogates
        for p=1:keys.n_permutations
            RPEAKS_intervals_p = [RPEAK_ts(1)+randn(1)*ecg_R2Rt_std RPEAKS_intervals(randperm(length(RPEAKS_intervals))) RPEAKS_intervals(randperm(length(RPEAKS_intervals)))]; %% slightly jitter first Rpeak
            RPEAK_ts_perm=cumsum(RPEAKS_intervals_p);
            RPEAK_ts_perm(RPEAK_ts_perm> trial_ends(end))=[];
            RPEAK_ts_perm(RPEAK_ts_perm<=0)=[];
            %% remove Rpeaks that landed in invalid intervals
            for iv=1:size(invalid_intervals,1)
                RPEAK_ts_perm(RPEAK_ts_perm>invalid_intervals(iv,1)+ecg_R2Rt_mean/2 & RPEAK_ts_perm<invalid_intervals(iv,2)-ecg_R2Rt_mean/2)=[];
            end
            sorted=sort_by_rpeaks_temp(RPEAK_ts_perm,SD_full_block,trial_onsets,trial_ends,keys,remove_ini);
            condition(tasktype).unit(u).permuations(p).trial(T+1:T+numel(sorted))=sorted;
        end
        
        %% The part following here is internal sanity check and should be turned off in general since there typically is no ECG data in the spike format
        if Sanity_check
            ECG_data=[pop.trial(tr).TDT_ECG1];
            ECG_time=cellfun(@(x) [1/x.TDT_ECG1_SR:1/x.TDT_ECG1_SR:numel(x.TDT_ECG1)/x.TDT_ECG1_SR]+x.TDT_ECG1_t0_from_rec_start+x.TDT_ECG1_tStart,trcell,'uniformoutput',false);
            CG_time=[ECG_time{:}];
            tt=0;
            clear ECG_to_plot
            for t=1:numel(RPEAK_ts)
                ECg_t_idx=CG_time>RPEAK_ts(t)-0.5 & CG_time<RPEAK_ts(t)+0.5;
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
    
    
%     
%     %% NANs if task type not present
%     SD=NaN(size(BINS));SD_SEM=NaN(size(BINS));SDPmean=NaN(size(BINS));SDPconf=[NaN(size(BINS));NaN(size(BINS))];
%     if numel(condition(1).unit) >= u
%         trial=condition(1).unit(u).trial;
%         SD=mean(vertcat(trial.PSTH),1);
%         SD_SEM=sterr(vertcat(trial.PSTH),1);
%         
%         % get mean and confidence intervals of shuffle predictor
%         for p=1:n_permutations
%             trial=condition(1).unit(u).permuations(p).trial;
%             SDP(p,:)=mean(vertcat(trial.PSTH),1);
%         end
%         SDPmean=nanmean(SDP,1);
%         SDPconf(1,:)=abs(prctile(SDP,2.5,1)-SDPmean);
%         SDPconf(2,:)=abs(prctile(SDP,97.5,1)-SDPmean);
%     end
    
    % separate for Rest and Task, group for Target
    [SD,BINS]=ecg_bna_condition_spike_density(condition(1),u,keys);
    Output.(target).Rest.SD(u,:)      = SD.SD_mean ;
    Output.(target).Rest.SD_SEM(u,:)  = SD.SD_SEM ;
    Output.(target).Rest.SDP(u,:)     = SD.SDPmean ;
    Output.(target).Rest.SDPCL(u,:)   = SD.SDPconf(1,:) ;
    Output.(target).Rest.SDPCu(u,:)   = SD.SDPconf(2,:) ;
    
    
    Output.(target).Rest.sig_all(u,:)       = SD.sig_all;
    Output.(target).Rest.sig(u,:)           = SD.sig;
    Output.(target).Rest.sig_FR_diff(u,:)   = SD.sig_FR_diff;
    Output.(target).Rest.sig_time(u,:)      = SD.sig_time;
    Output.(target).Rest.sig_n_bins(u,:)    = SD.sig_n_bins;
    Output.(target).Rest.sig_sign(u,:)      = SD.sig_sign;
    
    [SD,BINS]=ecg_bna_condition_spike_density(condition(2),u,keys);
    Output.(target).Task.SD(u,:)      = SD.SD_mean ; %% not good because
    Output.(target).Task.SD_SEM(u,:)  = SD.SD_SEM ; %% not good because
    Output.(target).Task.SDP(u,:)     = SD.SDPmean ;
    Output.(target).Task.SDPCL(u,:)   = SD.SDPconf(1,:) ;
    Output.(target).Task.SDPCu(u,:)   = SD.SDPconf(2,:) ;
    
    Output.(target).Task.sig_all(u,:)       = SD.sig_all;
    Output.(target).Task.sig(u,:)           = SD.sig;
    Output.(target).Task.sig_FR_diff(u,:)   = SD.sig_FR_diff;
    Output.(target).Task.sig_time(u,:)      = SD.sig_time;
    Output.(target).Task.sig_n_bins(u,:)    = SD.sig_n_bins;
    Output.(target).Task.sig_sign(u,:)      = SD.sig_sign;
%     
%     SD=NaN(size(BINS));SD_SEM=NaN(size(BINS));SDPmean=NaN(size(BINS));SDPconf=[NaN(size(BINS));NaN(size(BINS))];
%     if numel(condition(2).unit) >= u
%         trial=condition(2).unit(u).trial;
%         SD=mean(vertcat(trial.PSTH),1);
%         SD_SEM=sterr(vertcat(trial.PSTH),1);
%         
%         % get mean and confidence intervals of shuffle predictor
%         for p=1:n_permutations
%             trial=condition(2).unit(u).permuations(p).trial;
%             SDP(p,:)=mean(vertcat(trial.PSTH),1);
%         end
%         SDPmean=nanmean(SDP,1);
%         SDPconf(1,:)=abs(prctile(SDP,2.5,1)-SDPmean);
%         SDPconf(2,:)=abs(prctile(SDP,97.5,1)-SDPmean);
%     end
        
    if savePlot
        figure;
        title([unit_ID,'__',target ],'interpreter','none');
        hold on
        
        lineProps={'color','b','linewidth',1};
        shadedErrorBar(BINS,Output.(target).Rest.SD(u,:),Output.(target).Rest.SD_SEM(u,:),lineProps,1);
        lineProps={'color','b','linewidth',1,'linestyle',':'};
        shadedErrorBar(BINS,Output.(target).Rest.SDP(u,:),[Output.(target).Rest.SDPCL(u,:);Output.(target).Rest.SDPCu(u,:)],lineProps,1);
        lineProps={'color','r','linewidth',1};
        
        
        shadedErrorBar(BINS,Output.(target).Task.SD(u,:),Output.(target).Task.SD_SEM(u,:),lineProps,1);
        lineProps={'color','r','linewidth',1,'linestyle',':'};
        shadedErrorBar(BINS,Output.(target).Task.SDP(u,:),[Output.(target).Task.SDPCL(u,:);Output.(target).Task.SDPCu(u,:)],lineProps,1);
        
        y_lims=get(gca,'ylim');
        offset=NaN;factor=NaN;
        if Output.(target).Rest.sig_sign(u,:);
            offset=mean(y_lims);
            factor=diff(y_lims)/3;
        end
        to_plot=Output.(target).Rest.sig(u,:);to_plot(to_plot==0)=NaN;
        plot(BINS,to_plot*factor+offset,'b','linewidth',5);
        offset=NaN;factor=NaN;
        if Output.(target).Task.sig_sign(u,:);
            offset=mean(y_lims);
            factor=diff(y_lims)/3;
        end
        to_plot=Output.(target).Task.sig(u,:);to_plot(to_plot==0)=NaN;
        plot(BINS,to_plot*factor+offset,'r','linewidth',5);
        
        
        filename= [unit_ID, '__' target]; export_fig([basepath_to_save, filesep, filename], '-pdf','-transparent');
        close(gcf);
    end % pdf by run
    
    
    
end 
%% save output
save([basepath_to_save, filesep, session_info.session],'Output')
end

function out=sort_by_rpeaks_temp(RPEAK_ts,SD,trial_onsets,trial_ends,keys,remove_ini)
%% now the tricky part: sort by ECG peaks ...
if remove_ini
    %% reduce RPEAK_ts potentially ? (f.e.: longer than recorded ephys, inter-trial-interval?)
    during_trial_index = arrayfun(@(x) any(trial_onsets<=x+keys.PSTH_WINDOWS{1,3} & trial_ends>=x+keys.PSTH_WINDOWS{1,4}),RPEAK_ts);
    RPEAK_ts=RPEAK_ts(during_trial_index);
end
        
RPEAK_samples=round(RPEAK_ts/keys.PSTH_binwidth);
bins_before=round(keys.PSTH_WINDOWS{3}/keys.PSTH_binwidth);
bins_after=round(keys.PSTH_WINDOWS{4}/keys.PSTH_binwidth);

%% remove samples that would land outside
RPEAK_ts(RPEAK_samples<=-bins_before)=[];
RPEAK_samples(RPEAK_samples<=-bins_before)=[];
RPEAK_ts(RPEAK_samples>=numel(SD)-bins_after)=[];
RPEAK_samples(RPEAK_samples>=numel(SD)-bins_after)=[];

PSTH=NaN(numel(RPEAK_samples),bins_after-bins_before+1);

for s=1:numel(RPEAK_samples)
    PSTH(s,:)=SD(bins_before+RPEAK_samples(s):bins_after+RPEAK_samples(s));
end

%% get back trial structure
for t=1:numel(trial_onsets)
    rpeak_idx=RPEAK_ts>trial_onsets(t) & RPEAK_ts<trial_ends(t);
    out(t).PSTH=PSTH(rpeak_idx,:);
end
end

function SD=ecg_bna_spike_density(AT,trial_onsets,trial_ends,keys)
%% make SD across all trials appended (no average)!
switch keys.kernel_type
    case 'gaussian'
        Kernel=normpdf(-5*keys.gaussian_kernel:keys.PSTH_binwidth:5*keys.gaussian_kernel,0,keys.gaussian_kernel);
    case 'box'
        n_bins=round(2*keys.gaussian_kernel/keys.PSTH_binwidth);
        Kernel=ones(1,n_bins)/n_bins/keys.PSTH_binwidth; %%*1000 cause a one full spike in one 1ms bin means 1000sp/s locally
end
PSTH_time=trial_onsets(1):keys.PSTH_binwidth:trial_ends(end);
SD= conv(hist(AT,PSTH_time),Kernel,'same');
end

function [SD,BINS]=ecg_bna_condition_spike_density(condition,u,keys)
    BINS=(keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000;
    %% NANs if task type not present
    SD_mean=NaN(size(BINS));SD_SEM=NaN(size(BINS));SDPmean=NaN(size(BINS));SDPconf=[NaN(size(BINS));NaN(size(BINS))];
    if numel(condition.unit) >= u
        trial=condition.unit(u).trial;
        SD_mean=mean(vertcat(trial.PSTH),1);
        SD_SEM=sterr(vertcat(trial.PSTH),1);
        
        % get mean and confidence intervals of shuffle predictor
        for p=1:keys.n_permutations
            trial=condition.unit(u).permuations(p).trial;
            SDP(p,:)=mean(vertcat(trial.PSTH),1);
        end
        SDPmean=nanmean(SDP,1);
        SDPconf(1,:)=abs(prctile(SDP,2.5,1)-SDPmean);
        SDPconf(2,:)=abs(prctile(SDP,97.5,1)-SDPmean);
    end
    SD.SD_mean=SD_mean;
    SD.SD_SEM=SD_SEM;
    SD.SDPmean=SDPmean;
    SD.SDPconf=SDPconf;
    
    %% signficance
    sig_to_check=BINS>keys.significance_window(1)*1000 & BINS<keys.significance_window(2)*1000;
    sig_above=SD_mean>(SDPmean+SDPconf(2,:))&sig_to_check;
    sig_below=SD_mean<(SDPmean-SDPconf(1,:))&sig_to_check;
    % find longest period of significance
    sig_idx_above_start=find(diff([0 sig_above 0])>0);
    sig_idx_above_end=find(diff([0 sig_above 0])<0);
    [ma,m_ia]=max(sig_idx_above_end-sig_idx_above_start);
    sig_idx_below_start=find(diff([0 sig_below 0])>0);
    sig_idx_below_end=find(diff([0 sig_below 0])<0);
    [mb,m_ib]=max(sig_idx_below_end-sig_idx_below_start);
    
    if ma>mb
        sig_sign=1;
        m=ma;
        m_i=m_ia;
        sig_start_end=[sig_idx_above_start(m_i):sig_idx_above_end(m_i)-1];
    else
        sig_sign=-1;
        m=mb;
        m_i=m_ib;
        sig_start_end=[sig_idx_below_start(m_i):sig_idx_below_end(m_i)-1];
    end
    
    if isempty(m)
        sig_sign=0;
        m=NaN;
        t_start_end=[NaN NaN];
        maxidx=1;
        max_FR_diff=NaN;
    else
        [max_FR_diff, maxidx]=max(abs(SD_mean(sig_start_end)-SDPmean(sig_start_end)));
        t_start_end=BINS(sig_start_end);
    end
    sig_period=BINS>=t_start_end(1) & BINS<=t_start_end(end) ;
    max_time=t_start_end(maxidx);
    SD.sig_all=sig_above-sig_below;
    SD.sig=(sig_above&sig_period)-(sig_below&sig_period);
    SD.sig_FR_diff=max_FR_diff;
    SD.sig_time=max_time;
    SD.sig_n_bins=m;
    SD.sig_sign=sig_sign;
end




