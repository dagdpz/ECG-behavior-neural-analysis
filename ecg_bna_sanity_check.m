function Output=ecg_bna_sanity_check(by_block,T,Triggers,Blockoffsets,cfg)
Sanity_check=1; % ECG triggered ECG, turn off since typically there is no ECG data in the spike format

BINS=(cfg.analyse_states{1,3}:cfg.spk.PSTH_binwidth:cfg.analyse_states{1,4})*1000;
%Rblocks=[Triggers.block];




% ecg-related variables

for c=1:numel(cfg.condition)
    L=cfg.condition(c).name;
    %% get condition AND valid block trials only
    tr = ecg_bna_get_condition_trials(T, cfg.condition(c));
    %tr=ismember([T.block],blocks) & CT;
    Output.(L).NrTrials=sum(tr);
    if sum(tr)==0
        continue;
    end
    
    %         trcell=T(tr);
    %         [O]=ecg_bna_synchronize_with_trigger(Blockoffsets,trcell,popcell);
    %
    %         %% compute spike density as one continuous vector across all concatenated trials (hmmm there migth be a problem with interleaved trial types here)
    %         PSTH_time=O.trial_starts(1):cfg.spk.PSTH_binwidth:O.trial_ends(end);
    %         [SD_stream, RAST]=ecg_bna_spike_density(O.AT,PSTH_time,cfg.spk);
    %
    %
    %         events=cfg.analyse_states(:,1);
    %         for e=1:numel(events)
    %             E=events{e};
    %             cfg.window=[cfg.analyse_states{e,4:5}];
    %             trig=Triggers.(E);
    %             [during_trial_index,bins]=get_valid_PSTH_indexes(PSTH_time,O.trial_starts,O.trial_ends,cfg);
    %
    %
    %             realPSTHs         = compute_PSTH(trig,RAST,SD_stream,PSTH_time,during_trial_index,bins,'real',cfg);
    %             shuffledPSTH      = compute_PSTH(trig,RAST,SD_stream,PSTH_time,during_trial_index,bins,'shuffled',cfg);
    %
    %             SD                = do_statistics(realPSTHs,shuffledPSTH,bins,cfg.spk);
    %             Output.(L).(E)=SD;
    %
    %         end
    
    %% The part following here is internal sanity check and should be turned off in general since there typically is no ECG data in the spike format
    
    %% loop through blocks
    
    
    
    
    complete_trials=[by_block.trial];
    trcell={complete_trials(tr)};
    ECG_data=[complete_trials(tr).TDT_ECG1];
    ECG_SR={complete_trials(tr).TDT_ECG1_SR};
    from_rec={complete_trials(tr).TDT_ECG1_t0_from_rec_start};
    relative={complete_trials(tr).TDT_ECG1_tStart};
    samples={complete_trials(tr).TDT_ECG1};
    blocknum={complete_trials(tr).block};
    
    ECG_time=cellfun(@(x,y,s,r,b) [1/x:1/x:numel(y)/x]+s+r+Blockoffsets(b),ECG_SR,samples,from_rec,relative,blocknum,'uniformoutput',false);
    CG_time=[ECG_time{:}];
    
    
    events=cfg.analyse_states(:,1);
    for e=1:numel(events)
        E=events{e};
        cfg.window=[cfg.analyse_states{e,4:5}];
        trig=Triggers.(E);
        [during_trial_index,bins]=get_valid_PSTH_indexes(PSTH_time,O.trial_starts,O.trial_ends,cfg);
        
        
        realPSTHs         = compute_PSTH(trig,RAST,SD_stream,PSTH_time,during_trial_index,bins,'real',cfg);
        shuffledPSTH      = compute_PSTH(trig,RAST,SD_stream,PSTH_time,during_trial_index,bins,'shuffled',cfg);
        
        SD                = do_statistics(realPSTHs,shuffledPSTH,bins,cfg.spk);
        Output.(L).(E)=SD;
        
        
        ts=Triggers.R.ts;
        
        tt=0;
        clear ECG_to_plot
        for t=1:numel(ts)
            ECg_t_idx=CG_time>ts(t)-0.5 & CG_time<ts(t)+0.5;
            if round(sum(ECg_t_idx)-ECG_SR{1})~=0
                continue
            end
            tt=tt+1;
            ECG_to_plot(tt,:)=ECG_data(ECg_t_idx);
        end
        
        figure
        %filename=[unit_ID '_block_' num2str(block) '_ECG_average'];
        lineProps={'color','b','linewidth',1};
        shadedErrorBar(1:size(ECG_to_plot,2),mean(ECG_to_plot,1),[std(ECG_to_plot,1);std(ECG_to_plot,1)],lineProps,1);
    end
    
    %     export_fig([cfg.per_session_folder, filesep, filename], '-pdf','-transparent') % pdf by run
    %     close(gcf);
    %
    %     figure
    %     filename=[unit_ID '_block_' num2str(block) '_ECG_per5trials'];
    %     plot(ECG_to_plot(1:5:end,:)');
    %     export_fig([cfg.per_session_folder, filesep, filename], '-pdf','-transparent') % pdf by run
    %     close(gcf);
end
end
%% save output
% save([cfg.per_session_folder, filesep, pop.unit_ID, '_', pop.target],'Output')
% clear Output

%end
end
end

function [during_trial_index,bins]=get_valid_PSTH_indexes(PSTH_time,trial_onsets,trial_ends,cfg)


%% define which parts of the continous PSTH are during a trial
trial_onset_samples=ceil((trial_onsets-PSTH_time(1))/cfg.spk.PSTH_binwidth);
trial_ends_samples=floor((trial_ends-PSTH_time(1))/cfg.spk.PSTH_binwidth);
trial_onset_samples=trial_onset_samples+1;
trial_ends_samples=trial_ends_samples+1;

bins_before=round(cfg.window(1)/cfg.spk.PSTH_binwidth);
bins_after=round(cfg.window(2)/cfg.spk.PSTH_binwidth);
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

function out=compute_PSTH(trig,RAST,SD,PSTH_time,during_trial_index,bins,mode,cfg)

switch mode
    case 'real'
        
        ts=trig.ts;
        iv=trig.intervals;
    case 'shuffled'
        ts=trig.shuffled_ts;
        iv=trig.shuffled_intervals;
end

ts_samples=round((ts-PSTH_time(1))/cfg.spk.PSTH_binwidth);

%% preallocate "out" structure
n=size(ts_samples,1);
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
    valid = ts_samples(p,:)<=numel(SD) & ts_samples(p,:)>0;
    
    RT=ts(p,valid);      % take time-stamps
    RS=ts_samples(p,valid); % take samples
    RD=iv(p,valid);     % take durations
    
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
SD.SDPCL=abs(prctile(Shuffled.mean,2.5,1)-SDPmean);
SD.SDPCU=abs(prctile(Shuffled.mean,97.5,1)-SDPmean);

SD.SD=Real.mean;
SD.SD_STD=Real.std;
SD.SD_SEM=Real.SEM;
SD.SDP=SDPmean;
SD.NrEvents            = Real.n_events;
SD.raster              = logical(Real.raster); % logical replaces all numbers >0 with 1 and reduces memory load
SD.SD_diff             = SD.SD - SD.SDP; % spikes/s, difference between mean and jittered data
SD.SD_diff_normalized  = SD.SD_diff ./ SD.SDP *100; % percent signal change
SD.FR_ModIndex_diff    = max(SD.SD_diff) - min(SD.SD_diff); % difference between max and min FR
SD.FR_ModIndex_PcS     = max(SD.SD_diff_normalized) - min(SD.SD_diff_normalized); % difference between max and min % signal change


%not the place to do histograms id say (?)
SD.Rds                          = hist(Real.RDs{1},cfg.histbins); % put RR durations to plot those in the histograms later
SD.Rds_perm                     = hist([Shuffled.RDs{:}],cfg.histbins);
% Output.(L).FR                           = mean(SD_stream); %% not too sure this was the intended one...

%% signficance
sig_to_check=BINS>cfg.significance_window(1)*1000 & BINS<cfg.significance_window(2)*1000;
pos_diff=SD_mean-(SDPmean+SD.SDPCU);
neg_diff=SD_mean-(SDPmean-SD.SDPCL);
sig_above=pos_diff>0&sig_to_check;
sig_below=neg_diff<0&sig_to_check;

sig_idx_above_start	=find(diff([0 sig_above 0])>0);
sig_idx_above_end   =find(diff([0 sig_above 0])<0);
sig_idx_below_start =find(diff([0 sig_below 0])>0);
sig_idx_below_end   =find(diff([0 sig_below 0])<0);

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
