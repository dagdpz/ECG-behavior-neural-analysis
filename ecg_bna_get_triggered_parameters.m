function [triggered, significance, cluster_z_values] = ecg_bna_get_triggered_parameters(site,triggers, width_in_samples,cfg,observedD)
significance=struct;

switch cfg.to_trigger
    case 'LFP';
        FN_in={'pow','powbp','evoked','pha','pha','phabp'};
        FN={'pow','powbp','evoked','pha','itpc','itpcbp'};
    case 'MUA';
        FN_in={'evoked'};
        FN={'evoked'};
    case 'ECG';
        FN_in={'evoked'};
        FN={'evoked'};
end
datatype=lower(cfg.to_trigger); % mua or lfp

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
n_iterations=size(triggers,1);
NNN=triggers==0;
AAA=all(NNN,1);
TT=triggers(:,~AAA);
NNN=TT==0;
nRpeaks=sum(~NNN,2);

%% reducing rpeaks - this can lead to 0 remaining shuffled Rpeaks quite easily...
nRpeaks=sum(~NNN,2);
minvalid=min(nRpeaks);
nrows=size(NNN,2);
todrawfrom=1:nrows;
for sh = 1: n_iterations
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
%trigs=reshape(trigs,[n_shuffles,minvalid]); %% check dimensions
trigs=reshape(trigs,[minvalid,n_iterations]); %% check dimensions
trigs=trigs';

% %% alternative: remove indexes where more than half are invalid
%
% trigs=zeros(n_shuffles,max(nRpeaks));
%
% [~,sortix]=sort(-nRpeaks);
% for sh = sortix'
%     trigs(sh,1:nRpeaks(sh))=TT(sh,~TT(sh,:)==0);
% end
%
% toremove=sum(trigs==0,1)>size(trigs,1)/2;
% TT(:,toremove)=[];
% for tr = 1:size(TT,2)
%     F=TT(:,tr)==0;
%     F
%
%     todraw=todrawfrom;
%     todraw(F)=[];
%     randindex=rand(numel(todraw),1);
%     [~, randindex]=sort(randindex);
%     ndrawn=nrows-numel(F)-minvalid;
%     toremove(sh,:)=   [F todraw(randindex(1:ndrawn))]+(sh-1)*nrows;
% end
%
% %%%

tfs=[site.(datatype)];
%% check pre-allocation time saving (?)

for f=1:numel(FN_in)
    fn=FN{f};
    fni=FN_in{f};
    S2=size(tfs.(fni),1);
    S3=numel(w_samples);
    
    triggered.(fn).mean        = zeros(1,S2,S3);
    triggered.(fn).std         = zeros(1,S2,S3);
    triggered.(fn).conf95      = zeros(2,S2,S3);
%     if n_iterations==1
%         S1 = size(trigs,2);
%         triggered.(fn).complete        = zeros(S1,S2,S3);
%     end
        triggered.(fn).complete        = zeros(n_iterations,S2,S3);
end


%% NOW: loop through window samples to compute shuffled parameters all at once!

for s = 1: numel(w_samples)
    w=w_samples(s);
    t=trigs+w;
    triggered.ntriggers=minvalid;
    
    for f=1:numel(FN_in)
        fn=FN{f};
        fni=FN_in{f};
        AA=single(reshape(tfs.(fni)(:,t),size(tfs.(fni),1),size(t,1),size(t,2)));
        if ismember(fn,{'itpc','itpcbp'})
            BB=single(abs(mean(AA,3)));
        elseif ismember(fn,{'pha'})
            BB=single(angle(mean(AA,3)));
        else
            BB=single(mean(AA,3));
        end
        
        if n_iterations==1 && ismember(fn,{'evoked'})% && ismember(fni,{'ecg'})
            triggered.(fn).mean(1,:,s)        = single(mean(AA,3));
            triggered.(fn).std(1,:,s)         = single(std(AA,0,3));
            triggered.(fn).sterr(1,:,s)       = single(sterr(AA,3));
            % check dimensions
            triggered.(fn).conf95(1,:,s)      = single(prctile(AA,97.5,3));
            triggered.(fn).conf95(2,:,s)      = single(prctile(AA,2.5,3));
        else
            triggered.(fn).mean(1,:,s)        = single(mean(BB,2));
            triggered.(fn).std(1,:,s)         = single(std(BB,0,2));
            % check dimensions
            triggered.(fn).conf95(1,:,s)      = single(prctile(BB,97.5,2));
            triggered.(fn).conf95(2,:,s)      = single(prctile(BB,2.5,2));
        end
        
        
%         if n_iterations==1
%             triggered.(fn).complete(:,:,s)    = squeeze(AA)';
%         else
%             triggered.(fn).complete(:,:,s)    = BB';
%         end
        
        triggered.(fn).complete(:,:,s)    = BB';
        complete.(fn)(:,:,s)    = BB';
    end
end

if nargin>4 %% surrogate data,add signficance test here !
    for f=1:numel(FN)
        fn=FN{f};
        if ismember(fn,{'powbp','itpcbp'}) % treat each frequency independantly here!
            sig=false(size(observedD.(fn).mean));
            clus_z=[];
            for fb=1:size(complete.(fn),2)
                [signtmp, clus_t_tmp] = ecg_bna_test_vs_surrogate_means(complete.(fn)(:,fb,:),observedD.(fn).mean(:,fb,:));
                sig(:,fb,:)=signtmp;
                clus_z=[clus_z {clus_t_tmp}];
            end
        else
            [sig,clus_z] = ecg_bna_test_vs_surrogate_means(complete.(fn),observedD.(fn).mean);
        end
        significance.(fn)=sig;
        cluster_z_values.(fn)=clus_z;
    end
end


end


