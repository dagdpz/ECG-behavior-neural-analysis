
function [significance, cluster_T]=do_test2(data1,data2,fn)

numPermutation = 1000;
pthreshold=0.05;
clusterthreshold=0.05;



n1=size(data1,1); % trials x freq x time
n2=size(data2,1);
dataall=[data1;data2];

t_thresh_pos = (1-pthreshold/2);
t_thresh_neg = (pthreshold/2);
    
T_pos=zeros(1,numPermutation);
T_neg=zeros(1,numPermutation);

[~,RX]=sort(rand(size(dataall,1),numPermutation),1);
RX=[(1:(n1+n2))' ,RX];

sur=zeros(numPermutation+1,size(dataall,2),size(dataall,3));

for s = 1:numPermutation+1 % loop through shuffles - first
    d1=dataall(RX(1:n1,s),:,:);
    d2=dataall(RX(n1+1:end,s),:,:);
    
    if ismember(fn,{'itpc','itpcbp'})
        BB=single(abs(mean(d1,1)))-single(abs(mean(d2,1)));
    elseif ismember(fn,{'pha'})
        BB=single(angle(mean(d1,1)))-single(angle(mean(d2,1)));
    else
        BB=single(mean(d1,1))-single(mean(d2,1));
    end
    
%     T_map = simpleTTest2(d1,d2,1);
%     
%     pos_thresh = T_map>t_thresh_pos ;
%     neg_thresh = T_map<t_thresh_neg ;
%     
%     pos_clusters = bwconncomp(pos_thresh);
%     neg_clusters = bwconncomp(neg_thresh);
%     % max sum of the tmap
%     sumTpos = cellfun(@(x) sum(T_map(x)),pos_clusters.PixelIdxList);
%     sumTneg = cellfun(@(x) sum(T_map(x)),neg_clusters.PixelIdxList);
%     T_pos(s) = max([sumTpos 0]);
%     T_neg(s) = min([sumTneg 0]);
    
    sur(s,:,:)=squeeze(BB);
end
obs=sur(1,:,:);
sur(1,:,:)=[];

s_ix=1:numPermutation;
for s = s_ix % loop through shuffles - observed is in the end here so we can use the stored sig and tvals

    dat_s=squeeze(sur(s,:,:));
    est_mean=squeeze(mean(sur(s_ix~=s,:,:),1));
    est_var=squeeze(var(sur(s_ix~=s,:,:),1));
    T_map=(dat_s-est_mean)./sqrt(est_var); %% we actually care about pos/neg here
    
    conf_map=squeeze(sum(sur(s_ix~=s,:,:)<repmat(sur(s,:,:),[size(sur,1)-1,1,1]),1));
    conf_map=conf_map/(size(sur,1)-1);
    
    sigpos=conf_map>t_thresh_pos;
    signeg=conf_map<-t_thresh_neg;
    
    pos_clusters = bwconncomp(sigpos);
    neg_clusters = bwconncomp(signeg);
    sumTpos = cellfun(@(x) sum(T_map(x)),pos_clusters.PixelIdxList);
    sumTneg = cellfun(@(x) sum(T_map(x)),neg_clusters.PixelIdxList);
    T_pos(s) = max([sumTpos 0]);
    T_neg(s) = min([sumTneg 0]);
end
   

est_mean=squeeze(mean(sur,1));
est_var=squeeze(var(sur,1));
T_map=(squeeze(obs)-est_mean)./sqrt(est_var);

conf_map=squeeze(sum(sur<repmat(obs,[size(sur,1),1,1]),1));
conf_map=conf_map/(size(sur,1));

sigpos=conf_map>t_thresh_pos;
signeg=conf_map<-t_thresh_neg;

pos_clusters = bwconncomp(sigpos);
neg_clusters = bwconncomp(signeg);

% max sum of the tmap
sumTpos = cellfun(@(x) sum(T_map(x)),pos_clusters.PixelIdxList);
sumTneg = cellfun(@(x) sum(T_map(x)),neg_clusters.PixelIdxList);

% Compute p-values using vectorized operations
% For positive clusters
p_values_pos = zeros(1, numel(sumTpos));
for i = 1:numel(sumTpos)
    p_values_pos(i) = mean(T_pos >= sumTpos(i));
end
% For negative clusters
p_values_neg = zeros(1, numel(sumTneg));
for i = 1:numel(sumTneg)
    p_values_neg(i) = mean(T_neg <= sumTneg(i));
end

sigpixpos=vertcat(pos_clusters.PixelIdxList{p_values_pos<clusterthreshold});
sigpixneg=vertcat(neg_clusters.PixelIdxList{p_values_neg<clusterthreshold});

significance=zeros(size(T_map));
significance(sigpixpos)=1;
significance(sigpixneg)=-1;

cluster_T=[sumTpos(p_values_pos<clusterthreshold) sumTneg(p_values_neg<clusterthreshold)*-1];
end
