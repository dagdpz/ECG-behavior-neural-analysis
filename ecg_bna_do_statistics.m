function [SD,BINS]=ecg_bna_do_statistics(Real,Shuffled,BINS)
definition='max';

SD_mean=mean(Real);

SDPmean=nanmean(Shuffled,1);
SDPconf(1,:)=abs(prctile(Shuffled,2.5,1)-SDPmean);
SDPconf(2,:)=abs(prctile(Shuffled,97.5,1)-SDPmean);

SD.SD_mean=SD_mean;
SD.SD_STD=std(Real);
SD.SD_SEM=sterr(Real);
SD.SDPmean=SDPmean;
SD.SDPconf=SDPconf;

%% signficance
pos_diff=SD_mean-(SDPmean+SDPconf(2,:));
neg_diff=SD_mean-(SDPmean-SDPconf(1,:));

sig_above=SD_mean > SDPmean+SDPconf(2,:);
sig_below=SD_mean < SDPmean-SDPconf(1,:);

sig_idx_above_start=find(diff([0 sig_above 0])>0);
sig_idx_above_end=find(diff([0 sig_above 0])<0);
sig_idx_below_start=find(diff([0 sig_below 0])>0);
sig_idx_below_end=find(diff([0 sig_below 0])<0);

consec_above = diff([0 find(diff(sig_above)) numel(sig_above)]); % number of repetitions of the same element
consec_below =  diff([0 find(diff(sig_below)) numel(sig_below)]);

% take into account the circular nature of the data - consider clusters in
% the beginning and in the end of the cycle as one, compute its duration
if sig_above(1) && sig_above(end)
    consec_above(1) = consec_above(1)+consec_above(end);
end

if sig_below(1) && sig_below(end)
    consec_below(1) = consec_below(1)+consec_below(end);
end

if sig_above(1) % if starts with one
    clust_above = consec_above(1:2:end); % then take lengths of non-zero clusters
else
    clust_above = consec_above(2:2:end);
end

if sig_below(1)
    clust_below = consec_below(1:2:end);
else
    clust_below = consec_below(2:2:end);
end

% find longest period of significance
[ma,m_ia]=max(clust_above);
[mb,m_ib]=max(clust_below);
% find maximum deviation from surrogates
[~,m_abs_max] = max(abs(SD_mean - SDPmean));
[max_pos_diff,max_idx]=max(pos_diff);
[max_neg_diff,min_idx]=min(neg_diff);

m_imax=find(sig_idx_above_start<=max_idx & sig_idx_above_end>=max_idx);
m_imin=find(sig_idx_below_start<=min_idx & sig_idx_below_end>=min_idx);

if strcmp(definition,'max')
    sig_sign = sign(SD_mean(m_abs_max) - SDPmean(m_abs_max)); % sign at max abs diff
    if abs(max_pos_diff)>abs(max_neg_diff)
        m_i=m_imax;
        sig_start_end=[sig_idx_above_start(m_i):sig_idx_above_end(m_i)-1];
    else
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

