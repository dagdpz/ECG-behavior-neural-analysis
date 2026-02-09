function [significance, cluster_T] = ecg_bna_test_vs_surrogate_means(sur,obsmean)
pthreshold=0.05;
clusterthreshold=0.05;
n_sur=size(sur,1);
s_ix=1:n_sur;
minsize=3;

zThreshold = 1.96;

%tThreshold = (1-pthreshold)/2;


% find (for each surrogate) the max cluster significance

T_pos = cell(size(s_ix));
T_neg = cell(size(s_ix));
for s=s_ix
    dat_s=squeeze(sur(s,:,:));
    %T_map=(dat_s-est_mean)./est_std; %% we actually care about pos/neg here
    est_mean=squeeze(mean(sur(s_ix~=s,:,:),1));
    est_var=squeeze(var(sur(s_ix~=s,:,:),1));
    T_map=(dat_s-est_mean)./sqrt(est_var); %% we actually care about pos/neg here%
    %     conf_map=squeeze(sum(sur(s_ix~=s,:,:)<repmat(sur(s,:,:),[size(sur,1)-1,1,1]),1));
    %     conf_map=conf_map/(size(sur,1)-1)-0.5;
    
    sigpos=T_map>zThreshold;
    signeg=T_map<-zThreshold;
    
    CCpos = bwconncomp(sigpos);
    CCneg = bwconncomp(signeg);
    sumTpos = cellfun(@(x) sum(T_map(x)),CCpos.PixelIdxList);
    sumTneg = cellfun(@(x) sum(T_map(x)),CCneg.PixelIdxList);
    
    
    
    T_pos{s} = sort(sumTpos,'descend');
    T_neg{s} = sort(sumTneg,'ascend');
end


est_mean=squeeze(mean(sur,1));
est_std=squeeze(std(sur,0,1));

T_pos_n_dims=cellfun(@numel,T_pos);
T_neg_n_dims=cellfun(@numel,T_neg);

T_pos_max_dim=max(T_pos_n_dims);
T_neg_max_dim=max(T_neg_n_dims);

T_pos_full=zeros(n_sur,T_pos_max_dim);
for dim=1:T_pos_max_dim
    clus_idx=T_pos_n_dims>=dim;
    n=sum(clus_idx);
    T_pos_full(1:n,dim)=cellfun(@(x) x(dim),T_pos(clus_idx));
end
T_neg_full=zeros(n_sur,T_neg_max_dim);
for dim=1:T_neg_max_dim
    clus_idx=T_neg_n_dims>=dim;
    n=sum(clus_idx);
    T_neg_full(1:n,dim)=cellfun(@(x) x(dim),T_neg(clus_idx));
end

%T_map=(squeeze(obsmean)-est_mean)./sqrt(est_var);
T_map=(squeeze(obsmean)-est_mean)./est_std;

% conf_map=squeeze(sum(sur<repmat(obsmean,[size(sur,1),1,1]),1));
% conf_map=conf_map/size(sur,1)-0.5;

sigpos=T_map>zThreshold;
signeg=T_map<-zThreshold;

pos_clusters = bwconncomp(sigpos);
neg_clusters = bwconncomp(signeg);

numTpos = cellfun(@numel,pos_clusters.PixelIdxList);
numTneg = cellfun(@numel,neg_clusters.PixelIdxList);

pos_clusters.PixelIdxList(numTpos<minsize) = [];
neg_clusters.PixelIdxList(numTneg<minsize) = [];


sumTpos = cellfun(@(x) sum(T_map(x)),pos_clusters.PixelIdxList);
sumTneg = cellfun(@(x) sum(T_map(x)),neg_clusters.PixelIdxList);

% re-order descending
[sumTpos, oi]=sort(sumTpos,'descend');
pos_clusters.PixelIdxList=pos_clusters.PixelIdxList(oi);
[sumTneg, oi]=sort(sumTneg,'ascend');
neg_clusters.PixelIdxList=neg_clusters.PixelIdxList(oi);

p_values_pos = zeros(1, numel(sumTpos));
for i = 1:min(T_pos_max_dim,numel(sumTpos))
    p_values_pos(i) = mean(T_pos_full(:,i) >= sumTpos(i));
    if p_values_pos(i)>clusterthreshold
        % make all smaller clusters insignificant
        p_values_pos(1:numel(sumTpos)>i)=deal(1);
        break;
    end
end
% For negative clusters
p_values_neg = zeros(1, numel(sumTneg));
for i = 1:min(T_neg_max_dim,numel(sumTneg))
    p_values_neg(i) = mean(T_neg_full(:,i) <= sumTneg(i));
    if p_values_neg(i)>clusterthreshold
        % make all smaller clusters insignificant
        p_values_neg(1:numel(sumTneg)>i)=deal(1);
        break;
    end
end

sigpixpos=vertcat(pos_clusters.PixelIdxList{p_values_pos<clusterthreshold});
sigpixneg=vertcat(neg_clusters.PixelIdxList{p_values_neg<clusterthreshold});

significance=zeros(size(obsmean));
significance(sigpixpos)=1;
significance(sigpixneg)=-1;


cluster_T=[sumTpos(p_values_pos<clusterthreshold) sumTneg(p_values_neg<clusterthreshold)];
end
