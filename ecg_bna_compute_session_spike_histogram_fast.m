function Output=ecg_bna_compute_session_spike_histogram_fast(session_info,Rpeaks,ecg_bna_cfg)
savePlot = 1;
Sanity_check=0; % ECG triggered ECG, turn off since typically there is no ECG data in the spike format
ecg_bna_cfg.n_permutations=1000; %100;
ecg_bna_cfg.significance_window=[-0.25 0.25];
ecg_bna_cfg.PSTH_binwidth=0.01;
ecg_bna_cfg.kernel_type='gaussian';
ecg_bna_cfg.gaussian_kernel=0.02;

basepath_to_save=[session_info.SPK_fldr filesep 'per_unit'];
if ~exist(basepath_to_save,'dir')
    mkdir(basepath_to_save);
end

BINS=(ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000;
condition_labels={'Rest','Task'};
condition_colors={'b','r'};

% Rpeaks derived from concatenated ECG data [First_trial_INI.ECG1 trial.TDT_ECG1]
%load(session_info.Input_ECG);
load(session_info.Input_spikes);

% offset_blocks_Rpeak=0;
% for b=1:numel(out)
%     Rpeaks(b).block=out(b).nrblock_combinedFiles;
%     if isempty(out(b).nrblock_combinedFiles) || isempty(out(b).Rpeak_t(1)) || isempty(out(b).R2R_t)
%         Rpeaks(b).block=NaN;
%     end
%     RPEAK_ts=[out(b).Rpeak_t(1) intersect(out(b).Rpeak_t,out(b).R2R_t)];
%     RPEAKS_intervals=diff(RPEAK_ts);
%     ecg_R2Rt_mean=mean(RPEAKS_intervals);
%     idx_valid = RPEAKS_intervals<1.5*ecg_R2Rt_mean; %use mode or mean ?
%     nonval_idx=find([0, ~idx_valid]);
%     invalid_intervals=[NaN,NaN];
%     for iv=1:numel(nonval_idx)
%         invalid_intervals(iv,1)=RPEAK_ts(nonval_idx(iv)-1);
%         invalid_intervals(iv,2)=RPEAK_ts(nonval_idx(iv));
%     end
%     RPEAKS_intervals=RPEAKS_intervals(idx_valid);
%     ecg_R2Rt_std=std(RPEAKS_intervals);
%     N=ecg_bna_cfg.n_permutations;
%     [~,ix]=sort(rand(N,length(RPEAKS_intervals)),2);
%     perm=RPEAKS_intervals(ix);
%     
%     %% here we jitter all Rpeak intervals
%     RPEAKS_intervals_p = [repmat(RPEAK_ts(1),N,1) repmat(RPEAKS_intervals,N,1)+randn(N,length(RPEAKS_intervals))*ecg_R2Rt_std perm]; % jittering every interval!
%     RPEAK_ts_p=cumsum(RPEAKS_intervals_p,2);
%     
%     %% re-evaluating valid intervals... this is important to fix surrogates being higher due to periods of increased spiking that correlate with invalid Rpeaks
%     for iv=1:size(invalid_intervals,1)
%         RPEAK_ts_p(RPEAK_ts_p>invalid_intervals(iv,1)+ecg_R2Rt_mean/2 & RPEAK_ts_p<invalid_intervals(iv,2)-ecg_R2Rt_mean/2)=NaN;
%     end
%     RPEAK_ts_p(RPEAK_ts_p>max(RPEAK_ts)+max(RPEAKS_intervals))=NaN;
%     RPEAK_ts_p(:,all(isnan(RPEAK_ts_p),1))=[];
%     Rpeaks(b).RPEAK_ts=RPEAK_ts+offset_blocks_Rpeak(b);
%     Rpeaks(b).shuffled_ts=RPEAK_ts_p+offset_blocks_Rpeak(b);
%     Rpeaks(b).offset=offset_blocks_Rpeak(b);
%     offset_blocks_Rpeak(b+1)=offset_blocks_Rpeak(b)+max(RPEAK_ts)+max(RPEAKS_intervals)*2;
% end
% offset_blocks_Rpeak(end)=[];
% Rblocks=[Rpeaks.block];

offset_blocks_Rpeak=[Rpeaks.offset];
Rblocks=[Rpeaks.block];
for u=1:numel(population)
    pop=population(u);
    unit_ID=population(u).unit_ID;
    target =population(u).target;
    %% Make sure we only take overlapping blocks
    blocks_unit=unique([pop.trial.block]);
    blocks=intersect(blocks_unit,Rblocks);
    b=ismember(Rblocks,blocks);
    
    for tasktype=1:2
        L=condition_labels{tasktype};
        Output.(L).unit_ID{u}            = unit_ID;
        Output.(L).target{u}             = target;
        Output.(L).quantSNR(u,:)         = population(u).SNR_rating;
        Output.(L).Single_rating(u,:)    = population(u).Single_rating;
        Output.(L).stability_rating(u,:) = population(u).stability_rating;
        Output.(L).SD(u,:)               = NaN(size(BINS));
        Output.(L).SD_STD(u,:)           = NaN(size(BINS));
        Output.(L).SD_SEM(u,:)           = NaN(size(BINS));
        Output.(L).SDP(u,:)              = NaN(size(BINS));
        Output.(L).SDPCL(u,:)            = NaN(size(BINS));
        Output.(L).SDPCu(u,:)            = NaN(size(BINS));
        Output.(L).sig_all(u,:)          = zeros(size(BINS));
        Output.(L).sig(u,:)              = zeros(size(BINS));
        Output.(L).sig_FR_diff(u,:)      = NaN;
        Output.(L).sig_time(u,:)         = NaN;
        Output.(L).sig_n_bins(u,:)       = 0;
        Output.(L).sig_sign(u,:)         = 0;
        Output.(L).NrTrials(u,:)         = NaN;
        Output.(L).NrEvents(u,:)         = NaN;
        Output.(L).FR(u,:)               = NaN;
        Nooutput.(L).Rts=NaN;
        Nooutput.(L).Rts_perm={NaN;NaN};
        
        %% here we could potentially further reduce trials
        tr=ismember([pop.trial.block],blocks) & [pop.trial.type]==tasktype;
        trcell=num2cell(pop.trial(tr));
        
        % add trial onset time to each spike so its basically one stream again
        % also, make sure spikes aren't counted twice (because previous trial is appended in beginning;
        % removing overlapping spikes here            % add trial onset time         % add block separator
        arrival_times=cellfun(@(x) x.arrival_times(x.arrival_times>x.states_onset(x.states==1)) + x.TDT_ECG1_t0_from_rec_start+offset_blocks_Rpeak(Rblocks==x.block),trcell,'uniformoutput',false);
        trial_onsets=cellfun(@(x) x.TDT_ECG1_t0_from_rec_start+x.TDT_ECG1_tStart+offset_blocks_Rpeak(Rblocks==x.block),trcell);
        trial_ends=cellfun(@(x) x.states_onset(x.states==98)+x.TDT_ECG1_t0_from_rec_start+x.TDT_ECG1_tStart+offset_blocks_Rpeak(Rblocks==x.block),trcell,'uniformoutput',false); % no clue why this needs to be nonuniformoutput, it did work earlier so this is confusing...
        trial_ends=[trial_ends{:}];
        
        if ~numel(trial_onsets)>0
            continue; % out(1).nrblock_combinedFiles might be empty!
        end
        
        %% compute spike density as one continuous vector across all concatenated trials (hmmm there migth be a problem with interleaved trial types here)
        AT=vertcat(arrival_times{:});
        AT(AT>trial_ends(end))=[];
        RPEAK_ts=[Rpeaks(b).RPEAK_ts];
        RPEAK_ts_perm=[Rpeaks(b).shuffled_ts];
        [SD_all_trials, PSTH_time]=ecg_bna_spike_density(AT,trial_onsets,trial_ends,ecg_bna_cfg);
        
        %% compute raster
        AT = AT';
        bins = ecg_bna_cfg.analyse_states{3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{4};
        spikeTimesByRpeak = arrayfun(@(x) AT(AT > x + ecg_bna_cfg.analyse_states{3} & AT < x + ecg_bna_cfg.analyse_states{4})-x, RPEAK_ts(2:end), 'Uniformoutput', false); % figure out arrival times 0.5s before and after the current R-peak, center around R-peak (timepoint 0)
        rasterByRpeak = cellfun(@(x) histc(x, bins), spikeTimesByRpeak, 'UniformOutput', false);
        raster = cat(1, rasterByRpeak{:});
        raster = logical(raster); % values >1 will become 1 anyway
        
        %% define which parts of the continous PSTH are during a trial
        trial_onset_samples=ceil((trial_onsets-PSTH_time(1))/ecg_bna_cfg.PSTH_binwidth);
        trial_ends_samples=floor((trial_ends-PSTH_time(1))/ecg_bna_cfg.PSTH_binwidth);
        trial_onset_samples(trial_onset_samples==0)=1;
        during_trial_index=false(size(PSTH_time));
        for t=1:numel(trial_onset_samples)
            during_trial_index(trial_onset_samples(t):trial_ends_samples(t))=true;
        end
        
        realPSTHs=compute_PSTH(RPEAK_ts(2:end),SD_all_trials,PSTH_time,during_trial_index,ecg_bna_cfg);
        shuffledPSTH=compute_PSTH(RPEAK_ts_perm(:,2:end),SD_all_trials,PSTH_time,during_trial_index,ecg_bna_cfg);
        SD=do_statistics(realPSTHs,shuffledPSTH,BINS,ecg_bna_cfg);
        
        Output.(L).SD(u,:)            = SD.SD_mean;
        Output.(L).SD_STD(u,:)        = SD.SD_STD;
        Output.(L).SD_SEM(u,:)        = SD.SD_SEM ;
        Output.(L).SDP(u,:)           = SD.SDPmean ;
        Output.(L).SDPCL(u,:)         = SD.SDPconf(1,:) ;
        Output.(L).SDPCu(u,:)         = SD.SDPconf(2,:) ;
        Output.(L).sig_all(u,:)       = SD.sig_all;
        Output.(L).sig(u,:)           = SD.sig;
        Output.(L).sig_FR_diff(u,:)   = SD.sig_FR_diff;
        Output.(L).sig_time(u,:)      = SD.sig_time;
        Output.(L).sig_n_bins(u,:)    = SD.sig_n_bins;
        Output.(L).sig_sign(u,:)      = SD.sig_sign;
        Output.(L).NrTrials(u,:)      = sum(tr);
        Output.(L).NrEvents(u,:)      = realPSTHs.n_events;
        Output.(L).FR(u,:)            = mean(SD_all_trials); %% not too sure this was the intended one...
        Output.(L).raster             = raster;
        
        Nooutput.(L).Rts=realPSTHs.RTs{1};
        Nooutput.(L).Rts_perm=shuffledPSTH.RTs;
        
        
        %% The part following here is internal sanity check and should be turned off in general since there typically is no ECG data in the spike format
        if Sanity_check %% this needs to be fixed as well, this might be incorrect after least update...
            ECG_data=[pop.trial(tr).TDT_ECG1];
            ECG_time=cellfun(@(x) [1/x.TDT_ECG1_SR:1/x.TDT_ECG1_SR:numel(x.TDT_ECG1)/x.TDT_ECG1_SR]+x.TDT_ECG1_t0_from_rec_start+x.TDT_ECG1_tStart+offset_blocks_Rpeak(Rblocks==x.block),trcell,'uniformoutput',false);
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
    
    if savePlot
        figure; % raster
%         set(gcf, 'Position', [2054 104 1436 891])
        sgtitle(['Raster_' unit_ID,'__',target ],'interpreter','none');
        for tasktype=1:numel(condition_labels)
            subplot(1,2,tasktype)
            hold on
            L=condition_labels{tasktype};
            col=condition_colors{tasktype};
            if isfield(Output.(L), 'raster')
                % figure out how many R-peaks we have and if < 100, plot
                % everything; if > 100, choose only 100
                if size(Output.(L).raster,1) <= 100
                    stepSize = 1; % just step through all the Rpeaks one by one
                else
                    % choose median intervals of linearly spaced intervals
                    % between 1 and the number of Rpeaks
                    stepSize = median(diff(round(linspace(1, size(Output.(L).raster,1), 100))));
                end
                % plot spikes only for 100 Rpeaks or everything if we have
                % fewer
                a = 1; % introduce row counter
                for RpeakNum = 1:stepSize:size(Output.(L).raster,1)
                    x = Output.(L).raster(RpeakNum,:);
                    line([bins; bins], [x; 2*x]+a-1, 'Color', col, 'LineWidth', 1)
                    a = a + 1;
                end
                title([{L ': 50-ms bins, ', '100 R-peak intervals out of ' num2str(size(Output.(L).raster,1)) ' (linearly spaced)'}])
            end
            ylim([0 100])
            xlabel('Time from R-peak, ms')
            ylabel('Number of Plotted Row')
            box on
        end

        filename= ['Raster_' unit_ID, '__' target];
        export_fig([basepath_to_save, filesep, filename], '-pdf','-transparent') % pdf by run
        close(gcf);
        
        
        figure; %% PSTH
        title(['PSTH_' unit_ID,'__',target ],'interpreter','none');
        hold on
        for tasktype=1:numel(condition_labels)
            L=condition_labels{tasktype};
            col=condition_colors{tasktype};
            lineProps={'color',col,'linewidth',1};
            shadedErrorBar(BINS,Output.(L).SD(u,:),Output.(L).SD_SEM(u,:),lineProps,1);
            lineProps={'color',col,'linewidth',1,'linestyle',':'};
            shadedErrorBar(BINS,Output.(L).SDP(u,:),[Output.(L).SDPCu(u,:);Output.(L).SDPCL(u,:)],lineProps,1);
            ypos=NaN;
            if Output.(L).sig_sign(u,:)==-1;
                ypos=min(Output.(L).SD(u,:))*-1;
            elseif Output.(L).sig_sign(u,:)==1;
                ypos=max(Output.(L).SD(u,:));
            end
            to_plot=Output.(L).sig(u,:);to_plot(to_plot==0)=NaN;
            plot(BINS,to_plot*ypos,col,'linewidth',5);
        end
        y_lims=get(gca,'ylim');
        for tasktype=1:numel(condition_labels)
            L=condition_labels{tasktype};
            text(BINS(10),y_lims(2)-diff(y_lims)*tasktype*1/20, [L ': trials = ' ,num2str(Output.(L).NrTrials(u,:)), ' events = ' ,num2str(Output.(L).NrEvents(u,:)) ],'Color',condition_colors{tasktype});
        end
        ylabel('Firing rate');
        xlabel('time to Rpeak');
        filename= ['PSTH_' unit_ID, '__' target];
        export_fig([basepath_to_save, filesep, filename], '-pdf','-transparent') % pdf by run
        close(gcf);
        
        figure; %% interval distribution
        title(['Rpeak_intervals_' unit_ID,'__',target ],'interpreter','none');
        for tasktype=1:numel(condition_labels)
            L=condition_labels{tasktype};
            col=condition_colors{tasktype};
            subplot(2,1,tasktype);
            hold on
            histbins=0.2:0.02:0.8;
            H=hist(diff(Nooutput.(L).Rts),histbins);
            plot(histbins,H,'linewidth',2,'color',col);
            for p=1:numel(Nooutput.(L).Rts_perm)
                H(p,:)=hist(diff(Nooutput.(L).Rts_perm{p}),histbins);
            end
            lineProps={'color','k','linewidth',1,'linestyle',':'};
            shadedErrorBar(histbins,mean(H,1),std(H,1),lineProps,1);
            ylabel('N');
            xlabel('Rpeak interval');
            title('grey is mean and std of surrogates');
        end
        filename= ['Rpeak_intervals_' unit_ID, '__' target];
        export_fig([basepath_to_save, filesep, filename], '-pdf','-transparent') % pdf by run
        close(gcf);
    end % pdf by run
end
%% save output
save([basepath_to_save, filesep, session_info.session],'Output')
end

function out=compute_PSTH(RPEAK_ts,SD,PSTH_time,during_trial_index,cfg)
RPEAK_samples=round((RPEAK_ts-PSTH_time(1))/cfg.PSTH_binwidth);
bins_before=round(cfg.analyse_states{3}/cfg.PSTH_binwidth);
bins_after=round(cfg.analyse_states{4}/cfg.PSTH_binwidth);
bins=bins_before:bins_after;

%% remove samples that would land outside
RPEAK_ts(RPEAK_samples>=numel(SD)-bins_after)=NaN;
RPEAK_samples(RPEAK_samples>=numel(SD)-bins_after)=NaN;
RPEAK_ts(RPEAK_samples<=-bins_before)=NaN;
RPEAK_samples(RPEAK_samples<=-bins_before)=NaN;

%% distinguish between shuffled and real data?
for p=1:size(RPEAK_samples,1)
    RT=RPEAK_ts(p,~isnan(RPEAK_samples(p,:)));
    RS=RPEAK_samples(p,~isnan(RPEAK_samples(p,:)));
    RT=RT(during_trial_index(RS));
    RS=RS(during_trial_index(RS));
    PSTH = SD(bsxfun(@plus,RS',bins));
    out.mean(p,:)=mean(PSTH, 1);
    out.std(p,:) = std(PSTH, [], 1);
    out.SEM(p,:)=sterr(PSTH);
    out.n_events(p)=numel(RS);
    out.RTs{p}=RT;
end
end

function [SD,PSTH_time]=ecg_bna_spike_density(AT,trial_onsets,trial_ends,cfg)
%% make SD across all trials appended (no average)!
switch cfg.kernel_type
    case 'gaussian'
        Kernel=normpdf(-5*cfg.gaussian_kernel:cfg.PSTH_binwidth:5*cfg.gaussian_kernel,0,cfg.gaussian_kernel);
    case 'box'
        n_bins=round(2*cfg.gaussian_kernel/cfg.PSTH_binwidth);
        Kernel=ones(1,n_bins)/n_bins/cfg.PSTH_binwidth; %%*1000 cause a one full spike in one 1ms bin means 1000sp/s locally
end
PSTH_time=trial_onsets(1):cfg.PSTH_binwidth:trial_ends(end);
SD= conv(hist(AT,PSTH_time),Kernel,'same');
end

function [SD,BINS]=do_statistics(Real,Shuffled,BINS,cfg)
definition='max';

SD_mean=Real.mean;

SDPmean=nanmean(Shuffled.mean,1);
SDPconf(1,:)=abs(prctile(Shuffled.mean,2.5,1)-SDPmean);
SDPconf(2,:)=abs(prctile(Shuffled.mean,97.5,1)-SDPmean);

SD.SD_mean=Real.mean;
SD.SD_STD=Real.std;
SD.SD_SEM=Real.SEM;
SD.SDPmean=SDPmean;
SD.SDPconf=SDPconf;

%% signficance
sig_to_check=BINS>cfg.significance_window(1)*1000 & BINS<cfg.significance_window(2)*1000;
pos_diff=SD_mean-(SDPmean+SDPconf(2,:));
neg_diff=SD_mean-(SDPmean-SDPconf(1,:));
sig_above=pos_diff>0&sig_to_check;
sig_below=neg_diff<0&sig_to_check;

sig_idx_above_start=find(diff([0 sig_above 0])>0);
sig_idx_above_end=find(diff([0 sig_above 0])<0);
sig_idx_below_start=find(diff([0 sig_below 0])>0);
sig_idx_below_end=find(diff([0 sig_below 0])<0);

% find longest period of significance
[ma,m_ia]=max(sig_idx_above_end-sig_idx_above_start);
[mb,m_ib]=max(sig_idx_below_end-sig_idx_below_start);
% find maximum deviation from surrogates
[max_pos_diff,max_idx]=max(pos_diff.*sig_to_check);
[max_neg_diff,min_idx]=min(neg_diff.*sig_to_check);

m_imax=find(sig_idx_above_start<=max_idx & sig_idx_above_end>=max_idx);
m_imin=find(sig_idx_below_start<=min_idx & sig_idx_below_end>=min_idx);

if strcmp(definition,'max')
    if abs(max_pos_diff)>abs(max_neg_diff)
        sig_sign=1;
        m_i=m_imax;
        sig_start_end=[sig_idx_above_start(m_i):sig_idx_above_end(m_i)-1];
    else
        sig_sign=-1;
        m_i=m_imin;
        sig_start_end=[sig_idx_below_start(m_i):sig_idx_below_end(m_i)-1];
    end
elseif strcmp(definition,'dur')
    if ma>mb
        sig_sign=1;
        m_i=m_ia;
        sig_start_end=[sig_idx_above_start(m_i):sig_idx_above_end(m_i)-1];
    else
        sig_sign=-1;
        m_i=m_ib;
        sig_start_end=[sig_idx_below_start(m_i):sig_idx_below_end(m_i)-1];
    end
end

m=numel(sig_start_end);

if m==0
    sig_sign=0;
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
