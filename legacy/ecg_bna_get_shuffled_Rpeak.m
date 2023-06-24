function shuffled_Rpeaks = ecg_bna_get_shuffled_Rpeak(trials_lfp, nshuffles )
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

%% BE AWARE:
% this function only works properly if no trials have been excluded before,
% it assumes continous data (per block) is present!


T=0;
u_blocks=unique([trials_lfp.block]);
ts = trials_lfp(1).tsample;
for b=u_blocks
    block_trials=find([trials_lfp.block]==b);
    
    time_from_rec_start =[trials_lfp(block_trials).time_from_rec_start];
    ecg_peaks           =logical([trials_lfp(block_trials).ECG_spikes]);
    ecg_R2Rt            = diff(time_from_rec_start(ecg_peaks));
    RPEAK_ts            =time_from_rec_start(ecg_peaks);
    
    %% find invalid intervals...!!
    idx_valid = ecg_R2Rt<1.5*mode(ecg_R2Rt);
    invalid_intervals=[NaN,NaN];
    nonval_idx=find([0, ~idx_valid]);
    for iv=1:numel(nonval_idx)
        invalid_intervals(iv,1)=RPEAK_ts(nonval_idx(iv)-1);
        invalid_intervals(iv,2)=RPEAK_ts(nonval_idx(iv));
    end
    
    ecg_R2Rt=ecg_R2Rt(idx_valid);
    ecg_R2Rt_std=std(ecg_R2Rt);
    total_ecgR2R=sum(ecg_R2Rt);
    
    shuffled_ecg_peaks = false(nshuffles,length(ecg_peaks));
    
%     
%     
%     shuffled_ecg_peaks = cell(nshuffles,1);
%     [shuffled_ecg_peaks{:}]=deal(false(1,length(ecg_peaks)));
    if ~isempty(ecg_R2Rt) % no ecg recorded?
        
        R2R_repetitions=ceil((time_from_rec_start(end)-RPEAK_ts(1))/total_ecgR2R);
        ecg_R2Rt=repmat(ecg_R2Rt,1,R2R_repetitions);
        
        fprintf('Shuffling \n');
        for i = 1:nshuffles
            shuffled_ecg_p2pt = [RPEAK_ts(1)+randn(1)*ecg_R2Rt_std ecg_R2Rt(randperm(length(ecg_R2Rt)))]; %% add one (NOT SO!) random in the beginning
            shuffled_ecg_peak_times = cumsum(shuffled_ecg_p2pt);
            shuffled_ecg_peak_times = shuffled_ecg_peak_times(shuffled_ecg_peak_times < RPEAK_ts(end));
            
            %% remove Rpeaks that landed in invalid intervals
            for iv=1:size(invalid_intervals,1)
                shuffled_ecg_peak_times(shuffled_ecg_peak_times>invalid_intervals(iv,1) & shuffled_ecg_peak_times<invalid_intervals(iv,2))=[];
            end
            
            shuffled_ECG_peak_samples = round((shuffled_ecg_peak_times-time_from_rec_start(1)) / ts);
            shuffled_ECG_peak_samples(shuffled_ECG_peak_samples<=0)=[];
            shuffled_ecg_peaks(i,shuffled_ECG_peak_samples) = 1;
            %shuffled_ecg_peaks{i}(shuffled_ECG_peak_samples) = 1;
        end
    end
    
    
   %nsampless=arrayfun(@(x) numel(x.ECG_spikes),trials_lfp(block_trials));
    
    
    start_sample=1;
    for t=block_trials
        T=T+1;
        nsamples=numel(trials_lfp(t).ECG_spikes);
%        [shuffled_Rpeaks(:,T).ECG_spikes]=deal(shuffled_ecg_peaks{:}(start_sample:start_sample+nsamples-1));
        for i = 1:nshuffles %% this part here takes a bit too long...
            shuffled_Rpeaks(i,T).ECG_spikes=shuffled_ecg_peaks(i,start_sample:start_sample+nsamples-1);
        end
        start_sample=start_sample+nsamples;
    end
end

end

