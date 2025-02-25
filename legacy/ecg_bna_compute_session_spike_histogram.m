function Output=ecg_bna_compute_session_spike_histogram(trials,population,Rpeaks,cfg)
Sanity_check=0; % ECG triggered ECG, turn off since typically there is no ECG data in the spike format

BINS=(cfg.analyse_states{1,3}:cfg.spk.PSTH_binwidth:cfg.analyse_states{1,4})*1000;
Rblocks=[Rpeaks.block];

if cfg.spk.apply_exclusion_criteria
    load([cfg.unit_lists filesep cfg.spk.unit_list], 'unit_ids')
    
    pop_units  = {population.unit_ID};
    inc_ids    = ismember(pop_units, unit_ids);
    population = population(inc_ids);
    
end

for u=1:numel(population)
    tic
    pop=population(u);
    disp(['Processing ' pop.unit_ID])
    T=ph_get_unit_trials(pop,trials);
    
    %% Make sure we only take overlapping blocks
    blocks_unit=unique([pop.block]);
    blocks=intersect(blocks_unit,Rblocks);
    b=ismember(Rblocks,blocks);
    
    % preallocate 'Output' structure
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        % unit id and recording quality parameters
        Output.unit_ID               = pop.unit_ID;
        Output.target                = pop.target;
        Output.quantSNR              = pop.avg_SNR;
        Output.Single_rating         = pop.avg_single_rating;
        Output.stability_rating      = pop.avg_stability;
        % ecg-related variables
        Output.(L).SD                = single(nan(1, length(BINS)));
        Output.(L).SD_STD            = single(nan(1, length(BINS)));
        Output.(L).SD_SEM            = single(nan(1, length(BINS)));
        Output.(L).SDP               = single(nan(1, length(BINS)));
        Output.(L).SDPCL             = single(nan(1, length(BINS)));
        Output.(L).SDPCu             = single(nan(1, length(BINS)));
        Output.(L).sig_all           = single(zeros(1, length(BINS)));
        Output.(L).sig               = single(zeros(1, length(BINS)));
        Output.(L).sig_FR_diff       = single(nan(1));
        Output.(L).sig_time          = single(nan(1));
        Output.(L).sig_n_bins        = single(zeros(1));
        Output.(L).sig_sign          = single(zeros(1));
        Output.(L).NrTrials          = single(nan(1));
        Output.(L).NrEvents          = single(nan(1));
        Output.(L).FR                = single(nan(1));
        Output.(L).raster            = single(nan(1));
        Output.(L).Rts               = single(nan(1)); % RR ends
        Output.(L).Rds               = single(nan(1)); % RR durations
        Output.(L).Rds_perm          = single(nan(1));
        Output.(L).SDsubstractedSDP            = single(nan(1, length(BINS)));
        Output.(L).SDsubstractedSDP_normalized = single(nan(1, length(BINS)));
        Output.(L).FR_ModIndex_SubtrSDP        = single(nan(1));
        Output.(L).FR_ModIndex_PcS             = single(nan(1));
    end
    
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        %% get condition AND valid block trials only
        CT = ecg_bna_get_condition_trials(T, cfg.condition(c));
        tr=ismember([T.block],blocks) & CT;
        if sum(tr)==0
           continue; 
        end
        
        if (~isfield(Rpeaks, 'RPEAK_ts_insp') && cfg.process_Rpeaks_inhalation_exhalation) || (~isfield(Rpeaks, 'RPEAK_ts_exp') && cfg.process_Rpeaks_inhalation_exhalation)
            continue; % out(1).nrblock_combinedFiles might be empty! and even if there is 1 trial we're not processing
        end
        
        popcell=pop.trial(tr);
        trcell=T(tr);
        [trial_onsets, trial_ends,AT]=ecg_bna_synchronize_with_trigger(trcell,popcell,Rpeaks);
        
        
        %% compute spike density as one continuous vector across all concatenated trials (hmmm there migth be a problem with interleaved trial types here)
        PSTH_time=trial_onsets(1):cfg.spk.PSTH_binwidth:trial_ends(end);
        [SD_stream, RAST]=ecg_bna_spike_density(AT,PSTH_time,cfg.spk);
        
        %% retrieve R-peaks and their reshuffles
        RPEAK_ts=[Rpeaks(b).(['RPEAK_ts' cfg.condition(c).Rpeak_field])]; % relative to first trial INI (?!)
        RPEAK_ts_perm=[Rpeaks(b).(['shuffled_ts' cfg.condition(c).Rpeak_field])];
        RPEAK_dur = [Rpeaks(b).(['RPEAK_dur' cfg.condition(c).Rpeak_field])];
        RPEAK_dur_perm = [Rpeaks(b).(['shuffled_dur' cfg.condition(c).Rpeak_field])];
       
        [during_trial_index,bins]=get_valid_PSTH_indexes(PSTH_time,trial_onsets,trial_ends,cfg);
        
        realPSTHs         = compute_PSTH(RPEAK_ts,RPEAK_dur,RAST,SD_stream,PSTH_time,during_trial_index,bins,cfg);
        shuffledPSTH      = compute_PSTH(RPEAK_ts_perm,RPEAK_dur_perm,RAST,SD_stream,PSTH_time,during_trial_index,bins,cfg);
        SD                = do_statistics(realPSTHs,shuffledPSTH,bins,cfg.spk);
        
        Output.(L).SD                           = SD.SD_mean ;
        Output.(L).SD_STD                       = SD.SD_STD;
        Output.(L).SD_SEM                       = SD.SD_SEM ;
        Output.(L).SDP                          = SD.SDPmean ;
        Output.(L).SDPCL                        = SD.SDPconf(1,:) ;
        Output.(L).SDPCu                        = SD.SDPconf(2,:) ;
        Output.(L).sig_all                      = SD.sig_all;
        Output.(L).sig                          = SD.sig;
        Output.(L).sig_FR_diff                  = SD.sig_FR_diff;
        Output.(L).sig_time                     = SD.sig_time;
        Output.(L).sig_n_bins                   = SD.sig_n_bins;
        Output.(L).sig_sign                     = SD.sig_sign;
        Output.(L).NrTrials                     = sum(tr);
        Output.(L).NrEvents                     = realPSTHs.n_events;
        Output.(L).FR                           = mean(SD_stream); %% not too sure this was the intended one...
        Output.(L).raster                       = logical(realPSTHs.raster); % logical replaces all numbers >0 with 1 and reduces memory load
        %Output.(L).Rts                         = single(realPSTHs.RTs{1}); % unless we need this, dont save it!
        Output.(L).Rds                          = hist(realPSTHs.RDs{1},cfg.spk.histbins); % put RR durations to plot those in the histograms later
        Output.(L).Rds_perm                     = hist([shuffledPSTH.RDs{:}],cfg.spk.histbins);
        Output.(L).SDsubstractedSDP             = Output.(L).SD - Output.(L).SDP; % spikes/s, difference between mean and jittered data
        Output.(L).SDsubstractedSDP_normalized  = Output.(L).SDsubstractedSDP ./ Output.(L).SDP *100; % percent signal change
        Output.(L).FR_ModIndex_SubtrSDP         = max(Output.(L).SDsubstractedSDP) - min(Output.(L).SDsubstractedSDP); % difference between max and min FR
        Output.(L).FR_ModIndex_PcS              = max(Output.(L).SDsubstractedSDP_normalized) - min(Output.(L).SDsubstractedSDP_normalized); % difference between max and min % signal change
        
        clear realPSTHs shuffledPSTH SD
        
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
            export_fig([cfg.per_session_folder, filesep, filename], '-pdf','-transparent') % pdf by run
            close(gcf);
            
            figure
            filename=[unit_ID '_block_' num2str(block) '_ECG_per5trials'];
            plot(ECG_to_plot(1:5:end,:)');
            export_fig([cfg.per_session_folder, filesep, filename], '-pdf','-transparent') % pdf by run
            close(gcf);
        end
    end
    %% save output
	save([cfg.per_session_folder, filesep, pop.unit_ID, '_', pop.target],'Output')
    clear Output
    toc
end
end

function [during_trial_index,bins]=get_valid_PSTH_indexes(PSTH_time,trial_onsets,trial_ends,cfg)


%% define which parts of the continous PSTH are during a trial
trial_onset_samples=ceil((trial_onsets-PSTH_time(1))/cfg.spk.PSTH_binwidth);
trial_ends_samples=floor((trial_ends-PSTH_time(1))/cfg.spk.PSTH_binwidth);
trial_onset_samples=trial_onset_samples+1;
trial_ends_samples=trial_ends_samples+1;

bins_before=round(cfg.analyse_states{3}/cfg.spk.PSTH_binwidth);
bins_after=round(cfg.analyse_states{4}/cfg.spk.PSTH_binwidth);
bins=bins_before:bins_after;

during_trial_index=false(size(PSTH_time));
drop_samples = trial_onset_samples < 1 | trial_ends_samples < 1; % in very rare cases samples have negative values, drop those (??)
trial_onset_samples = trial_onset_samples(~drop_samples);
trial_ends_samples = trial_ends_samples(~drop_samples);
for t=1:numel(trial_onset_samples)
    during_trial_index(trial_onset_samples(t)-(bins_before-1):trial_ends_samples(t)-bins_after)=true;
end
during_trial_index(1:-bins_before)=false;
during_trial_index(end-bins_after:end)=false;
% 
% %% remove samples that would land outside
% rpeaks2exclude = ...
%     RPEAK_samples>=numel(SD)-bins_after | RPEAK_samples<=-bins_before;
% 
% RPEAK_ts(rpeaks2exclude)       = NaN;
% RPEAK_samples(rpeaks2exclude)  = NaN;
% RPEAK_dur(rpeaks2exclude)      = NaN;
end

function out=compute_PSTH(RPEAK_ts,RPEAK_dur,RAST,SD,PSTH_time,during_trial_index,bins,cfg)

RPEAK_samples=round((RPEAK_ts-PSTH_time(1))/cfg.spk.PSTH_binwidth);

%% preallocate "out" structure
n=size(RPEAK_samples,1);
m=numel(bins);

out.raster    = NaN;
out.mean      = nan([n m]);
out.std       = nan([n m]);
out.SEM       = nan([n m]);
out.n_events  = NaN([n 1]);
out.RTs       = {};
out.RDs       = {};

%% loop through rows of RPEAK_samples: 1 row for real, nReshuffles rows of reshuffled data
for p=1:n
% %% remove samples that would land outside
    rpeaks2_valid = RPEAK_samples(p,:)<=numel(SD) & RPEAK_samples(p,:)>0;
    
    RT=RPEAK_ts(p,rpeaks2_valid);      % take time-stamps
    RS=RPEAK_samples(p,rpeaks2_valid); % take samples
    RD=RPEAK_dur(p,rpeaks2_valid);     % take durations
    
    within_trial_idx = during_trial_index(RS);
    
    RT=RT(within_trial_idx);
    RS=RS(within_trial_idx);
    RD=RD(within_trial_idx);
    
    idx_by_event = bsxfun(@plus,RS',bins); % bsxfun produces a RPEAK-(nBins-1) matrix with samples taken from SDF for each R-peak
    if n == 1
        out.raster = RAST(idx_by_event);
    end
    PSTH = SD(idx_by_event); 
    out.mean(p,:)=mean(PSTH, 1);
    out.std(p,:) = std(PSTH, [], 1);
    out.SEM(p,:)=sterr(PSTH);
    out.n_events(p)=numel(RS);
    out.RTs{p}=RT;
    out.RDs{p}=RD;
end
end

function [SD,RAST,PSTH_time]=ecg_bna_spike_density(AT,PSTH_time,cfg)
%% make SD across all trials appended (no average)!
switch cfg.kernel_type
    case 'gaussian'
        Kernel=normpdf(-5*cfg.gaussian_kernel:cfg.PSTH_binwidth:5*cfg.gaussian_kernel,0,cfg.gaussian_kernel);
    case 'box'
        n_bins=round(2*cfg.gaussian_kernel/cfg.PSTH_binwidth);
        Kernel=ones(1,n_bins)/n_bins/cfg.PSTH_binwidth; %%*1000 cause a one full spike in one 1ms bin means 1000sp/s locally
end
RAST = hist(AT,PSTH_time);
SD= conv(RAST,Kernel,'same');
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
