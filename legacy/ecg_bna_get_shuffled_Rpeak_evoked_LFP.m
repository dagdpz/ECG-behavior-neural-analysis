function shuffled_LFP_evoked = ecg_bna_get_shuffled_Rpeak_evoked_LFP( trials_lfp, cfg_state, nshuffles )
% ecg_bna_get_shuffled_Rpeak_evoked_ECG - computes the evoked LFP for
% a specified  time window around randomly shuffled triggers for given
% trials (usually trials belonging to a condition) in a session. The random
% triggers are obtained by randomly shuffling the Rpeaks based on their
% peak-to-peak intervals
%
% USAGE:
%	shuffled_ecg_evoked = ecg_bna_get_shuffled_Rpeak_evoked_ECG(
%	trials_ecg, cfg_ecg, nshuffles )
%
% INPUTS:
%       trials_lfp      - 1xN struct containing LFP and Rpeak data of N
%       trials
%       cfg_state       - a cell array specifying time window around the
%       trigger during which evoked response should be obtained
%       nshuffles       - number of times the Rpeak triggers should be
%       randomly shuffled
% OUTPUTS:
%		shuffled_LFP_evoked   - struct containing shuffled Rpeak triggered
%		evoked LFP from the given trials (mean and standard deviation
%		across shuffles)
%
% REQUIRES:	ft_spiketriggeredaverage (FieldTrip toolbox)
%
% See also ecg_bna_compute_session_Rpeak_evoked_LFP, ecg_bna_get_Rpeak_evoked_LFP, ecg_bna_get_shuffled_Rpeak_evoked_ECG

% whether to plot the distribution of ECG peak-to-peak interval
%distplot = false;

shuffled_LFP_evoked.lfp = {};
shuffled_LFP_evoked.lfp_time = {};
shuffled_LFP_evoked.state = cfg_state{1};
shuffled_LFP_evoked.state_name = cfg_state{2};

trials_original=trials_lfp;



trials_lfp=trials_original;
%% sort by block first?
%% for concatinated blocks, make sure timeline is monotoneous
u_blocks=unique([trials_lfp.block]);
end_time_prev_block=0;
for b=u_blocks
    tr=find([trials_lfp.block]==b);
    for t=tr
        trials_lfp(t).time_from_rec_start=trials_lfp(t).time_from_rec_start+end_time_prev_block;
    end
    end_time_prev_block=trials_lfp(t).time_from_rec_start(end);
end

use_t=[trials_lfp.completed] & ~arrayfun(@(x) isempty(x.ECG_spikes)||isempty(x.lfp_data),trials_lfp);
time=[trials_lfp(use_t).time];
time_from_rec_start=[trials_lfp(use_t).time_from_rec_start]; 
%% WHY ARE THERE SOMETIMES HUGE JUMPS (because conditions?)
%% NOT QUITE SURE HERE!!
% what to do with concatination across blocks??


ecg_peaks=logical([trials_lfp(use_t).ECG_spikes]);
LFP_raw= [trials_lfp(use_t).lfp_data];
fsample=[trials_lfp(use_t).fsample];
%blocks= [trials_lfp(use_t).block];

ecg_R2Rt = diff(time_from_rec_start(ecg_peaks));
RPEAK_ts=time_from_rec_start(ecg_peaks);

%% find invalid intervals...!!
idx_valid = ecg_R2Rt<1.5*mode(ecg_R2Rt);
invalid_intervals=[NaN,NaN];
nonval_idx=find([0, ~idx_valid]);
for iv=1:numel(nonval_idx)
    invalid_intervals(iv,1)=RPEAK_ts(nonval_idx(iv)-1);
    invalid_intervals(iv,2)=RPEAK_ts(nonval_idx(iv));
end

% for t = 1:length(trials_lfp)
%     if trials_lfp(t).completed == 0
%         continue;
%     end
%     % check if ECG spike data exists for this trial
%     if isempty(trials_lfp(t).ECG_spikes)
%         continue;
%     end
%     % get the ECG raw data for the trial
%     if isempty(trials_lfp(t).lfp_data)
%         continue;
%     else
%         trial_LFP_raw = trials_lfp(t).lfp_data;
%         %trial_LFP_raw(~trials_lfp(t).ECG_valid) = [];
%         LFP_raw = [LFP_raw, trial_LFP_raw];
%     end
%
%     %timestamps = [timestamps, trials_ecg(t).time];
%     ecg_peaks = [ecg_peaks, ...
%         logical(trials_lfp(t).ECG_spikes)];
% %     ecg_R2Rt = [ecg_R2Rt, diff(trials_lfp(t).time(trials_lfp(t).ECG_spikes))];
%
%
% %% this part wasnt working, as  ECG_b2btime wasnt even defined
% %         trials_lfp(t).ECG_b2btime(trials_lfp(t).ECG_spikes & ...
% %         ~isnan(trials_lfp(t).ECG_b2btime))];
%
%
% end

%mean_ecg_R2Rt = mean(ecg_R2Rt);
ecg_R2Rt=ecg_R2Rt(idx_valid);
ecg_R2Rt_std=std(ecg_R2Rt);
total_ecgR2R=sum(ecg_R2Rt);

R2R_repetitions=ceil((time_from_rec_start(end)-RPEAK_ts(1))/total_ecgR2R);
ecg_R2Rt=repmat(ecg_R2Rt,1,R2R_repetitions);


ts = trials_lfp(1).tsample;
timestamps = linspace(0, ts*(length(LFP_raw) - 1), length(LFP_raw)); 

%preparing input for ecg_bna_get_Rpeak_evoked_LFP

trials_LFP_shuffled.time=timestamps;
trials_LFP_shuffled.lfp_data=LFP_raw;
trials_LFP_shuffled.trialperiod=[timestamps(1) timestamps(end)];
trials_LFP_shuffled.fsample=fsample(1);

Rpeak_evoked_LFP = struct();
fprintf('Shuffling \n');
for i = 1:nshuffles
    %fprintf('Shuffle %g\n', i);
    shuffled_ecg_p2pt = [RPEAK_ts(1)+randn(1)*ecg_R2Rt_std ecg_R2Rt(randperm(length(ecg_R2Rt)))]; %% add one (NOT SO!) random in the beginning
    shuffled_ecg_peaks = false(size(ecg_peaks));
    shuffled_ecg_peak_times = cumsum(shuffled_ecg_p2pt);
    shuffled_ecg_peak_times = shuffled_ecg_peak_times(shuffled_ecg_peak_times < RPEAK_ts(end));
    
    %% remove Rpeaks that landed in invalid intervals
    for iv=1:size(invalid_intervals,1)
        shuffled_ecg_peak_times(shuffled_ecg_peak_times>invalid_intervals(iv,1) & shuffled_ecg_peak_times<invalid_intervals(iv,2))=[];
    end
    for p=1:numel(shuffled_ecg_peak_times)
       [~,shuffled_ECG_peak_samples(p)]=min(abs(time_from_rec_start-shuffled_ecg_peak_times(p)));
    end
%     shuffled_ECG_peak_samples = round((shuffled_ecg_peak_times-time_from_rec_start(1)) / ts);
%     shuffled_ECG_peak_samples(shuffled_ECG_peak_samples==0)=[];
    % randomly shift shuffled ECG peaks why ??
    %     shuffled_ECG_peak_samples = shuffled_ECG_peak_samples + round((2*rand(size(shuffled_ECG_peak_samples)) - 1) * (mean_ecg_R2Rt/ts));
    %     shuffled_ECG_peak_samples(shuffled_ECG_peak_samples > length(ecg_peaks) | shuffled_ECG_peak_samples < 1) = [];
    % shuffled ECG peaks
    shuffled_ecg_peaks(shuffled_ECG_peak_samples) = 1;
    
    
    
    trials_LFP_shuffled.ECG_spikes=shuffled_ecg_peaks;
    ecg_triggered_evoked = ecg_bna_get_Rpeak_evoked_LFP( trials_LFP_shuffled, cfg_state );
    Rpeak_evoked_LFP(i).lfp_time = ecg_triggered_evoked.lfp_time;
    Rpeak_evoked_LFP(i).mean = ecg_triggered_evoked.mean;
    
    
    
    %% this part here does not work, becuase ft_spiketriggeredaverage normalized (subtracting mean) within every segment!
    %     ft_data_all = struct();
    %     ft_data_all.trial = {[LFP_raw; shuffled_ecg_peaks]};
    %     ft_data_all.time = {timestamps}; % dont htink this makes a difference ??
    %     %ft_data_all.time = {timestamps};
    %     ft_data_all.label = {'LFP', 'Rpeak'};
    %
    %     % now find the evoked LFP
    %     % evoked LFP average
    %     cfg                 = [];
    %     cfg.timwin          = [cfg_state{3} cfg_state{4}];
    %     cfg.channel         = ft_data_all.label(1); % LFP channel
    %     cfg.spikechannel    = ft_data_all.label(2); % ECG peak
    %
    %     ecg_based_sta       = ft_spiketriggeredaverage(cfg, ft_data_all);
    %     Rpeak_evoked_LFP(i).lfp_time = ecg_based_sta.time;
    %     Rpeak_evoked_LFP(i).mean = ecg_based_sta.avg;
end

shuffled_LFP_evoked.lfp = cat(1, Rpeak_evoked_LFP.mean);
shuffled_LFP_evoked.dimord = 'nshuffles_time';
shuffled_LFP_evoked.lfp_time = Rpeak_evoked_LFP(1).lfp_time;
shuffled_LFP_evoked.mean = nanmean(cat(1, Rpeak_evoked_LFP.mean), 1);
shuffled_LFP_evoked.std = nanstd(cat(1, Rpeak_evoked_LFP.mean), 0, 1);
end

