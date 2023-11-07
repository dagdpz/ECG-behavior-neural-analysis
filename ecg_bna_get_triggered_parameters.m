function triggered = ecg_bna_get_triggered_parameters(site,triggers, width_in_samples)
% ecg_bna_get_triggered_split_shuffled - gets the time-frequency
% pow,phs spectrogram and phsBP for a specified time window around all shuffled Rpeak onset for a single site for
% given trials (usually trials belonging to a condition) in a session
%
% USAGE:
%	triggered = ecg_bna_get_triggered_split_shuffled( site_lfp,
%	cond_ecg, state, ecg_bna_cfg )
%
% INPUTS:
%		site_lfp            - struct containing LFP data for all trials of
%		a session from a single site
%       cond_ecg         - ecg data of trials whose time freq spectrogram
%       are to be obtaied
%       state               - a cell array specifying time window during
%       which LFP tfs should be obtained
%       ecg_bna_cfg         - struct containing settings
%       Required fields:
%           tfr.timestep            - width of a timebin in LFP TFS in
%           number of LFP samples
%           baseline_method         - method used for baseline
%           normalization of LFP TFS
%           baseline_perturbation   - whether to use
%           control(0)/inactivation(1) trials for baseline
%           baseline_use_choice_trial - whether to use
%           choice(1)/instructed(0) trials for baseline
%           baseline_use_type       - type number mof trials to be used for
%           baseline
%           baseline_use_effector   - effector value of trials to be used
%           for baseline
% OUTPUTS:
%		triggered   - struct containing LFP tfs for all given
%		trials
%
% REQUIRES:	lfp_tfa_baseline_normalization ???
%
% See also ecg_bna_get_Rpeak_based_STA, ecg_bna_get_Rpeak_evoked_LFP,
% ecg_bna_get_shuffled_Rpeak_evoked_LFP




triggered = struct;
w_samples=width_in_samples(1):width_in_samples(2);

%% get rid of some events to get a matrix (!!!!)
n_shuffles=size(triggers,1);
NNN=triggers==0;
AAA=all(NNN,1);
TT=triggers(:,~AAA);
NNN=TT==0;
nRpeaks=sum(~NNN,2);
minvalid=min(nRpeaks);
nrows=size(NNN,2);
todrawfrom=1:nrows;

for sh = 1: n_shuffles
    F=find(NNN(sh,:));
    todraw=todrawfrom;
    todraw(F)=[];
    randindex=rand(numel(todraw),1);
    [~, randindex]=sort(randindex);
    ndrawn=nrows-numel(F)-minvalid;
    toremove(sh,:)=   [F todraw(randindex(1:ndrawn))]+(sh-1)*nrows;
end

trigs=TT';
trigs(toremove)=[];
trigs=reshape(trigs,[n_shuffles,minvalid]); %% check dimensions

%% NOW: loop through window samples to compute shuffled parameters all at once!

tfs=[site.tfs];
for s = 1: numel(w_samples)
    w=w_samples(s);
    t=trigs+w;
    triggered.ntriggers=minvalid;
    
    FN_in={'pow','powbp','lfp','pha','phabp'};
    FN={'pow','powbp','lfp','itpc','itpcbp'};
    for f=1:numel(FN_in)
        fn=FN{f};
        fni=FN_in{f};
        
        AA=      reshape(tfs.(fni)(:,t),size(tfs.(fni),1),size(t,1),size(t,2));
        if ismember(fni,{'pha','phabp'})
            BB=abs(mean(AA,3));
        else
            BB=mean(AA,3);
        end
        
        triggered.(fn).mean(1,:,s)        = mean(BB,2);
        triggered.(fn).std(1,:,s)         = std(BB,0,2);
        %% this one not sure what to put
        %triggered.(fn).std_mean(:,s)    = mean(abs(mean(AA,3)),2);
        %% check dimensions
        triggered.(fn).conf95(1,:,s)      = prctile(BB,97.5,2);
        triggered.(fn).conf95(2,:,s)      = prctile(BB,2.5,2);
    end
    
    %% differentiate between
    %%      a) mean of std and std of (shuffle) means
    %%      a) mean confidence interval and confidence interval of (shuffle) means
    
end

end

