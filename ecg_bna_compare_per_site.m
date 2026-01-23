function difference = ecg_bna_compare_per_site( C1, C2, cfg)

significance=struct;

switch cfg.to_trigger
    case 'LFP';
        FN={'pow','powbp','evoked','pha','itpc','itpcbp'};
    case 'MUA';
        FN={'evoked'};
    case 'ECG';
        FN={'evoked'};
end
datatype=lower(cfg.to_trigger); % mua or lfp


%tfs=[site.(datatype)];
%% check pre-allocation time saving (?)

for f=1:numel(FN)
    fn=FN{f};
%     
%     S2=size(tfs.(fni),1);
%     
%     triggered.(fn).mean        = zeros(1,S2,S3);
%     triggered.(fn).std         = zeros(1,S2,S3);
%     triggered.(fn).conf95      = zeros(2,S2,S3);
 observed=C1.observed.(fn).mean   - C2.observed.(fn).mean;
 surrogate=C1.surrogate.(fn).complete   - C2.surrogate.(fn).complete;
    
 difference.time=C1.time;
 difference.tfr_time=C1.tfr_time;
 difference.observed.ntriggers=(C1.observed.ntriggers+C2.observed.ntriggers)/2;
 difference.surrogate.ntriggers=(C1.surrogate.ntriggers+C2.surrogate.ntriggers)/2;
 
 difference.observed.(fn).mean=observed;
 difference.surrogate.(fn).mean=mean(surrogate,1); 
 difference.surrogate.(fn).conf95(1,:,:)      = single(prctile(surrogate,97.5,1));
 difference.surrogate.(fn).conf95(2,:,:)      = single(prctile(surrogate,2.5,1));       
 if ismember(fn,{'evoked'})
 difference.observed.(fn).sterr=sqrt(C1.observed.(fn).sterr.^2 + C2.observed.(fn).sterr.^2);
 end
 
 %% this is not good yet (?)
 difference.observed.(fn).std=sqrt(C1.observed.(fn).std.^2 + C2.observed.(fn).std.^2);
 difference.surrogate.(fn).std=sqrt(C1.surrogate.(fn).std.^2 + C2.surrogate.(fn).std.^2);
 
 
if ismember(fn,{'powbp','itpcbp'}) % treat each frequency independently here!
    sig=false(size(observed));
    clus_z=[];
    for fb=1:size(observed,2)
        [signtmp, clus_t_tmp] = ecg_bna_test_vs_surrogate_means(surrogate(:,fb,:),observed(:,fb,:));
        sig(:,fb,:)=signtmp;
        clus_z=[clus_z {clus_t_tmp}];
    end
else
    [sig,clus_z] = ecg_bna_test_vs_surrogate_means(surrogate,observed);
end
significance.(fn)=sig;
cluster_z_values.(fn)=clus_z;
    

    
end
difference.significance=significance;
difference.cluster_z_values=cluster_z_values;





end
