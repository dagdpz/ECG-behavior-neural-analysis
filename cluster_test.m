function [mean_diff,significance,t_map]=cluster_test(c1,c2,pairedorunpaired)



%%
% Permutation Based on Lukas function
numPermutation = 1000;
pthreshold = 0.01;
clusterthreshold = 0.01;

n1=size(c1,1);
n2=size(c2,1);

switch pairedorunpaired
    case 'paired'
        c_all=c1-c2;
    case 'unpaired'
        c_all=[c1;c2];
end

T_pos=zeros(1,numPermutation);
T_neg=zeros(1,numPermutation);


stderr_diff = std(c_all, 0, 1) / sqrt(size(c_all,1));

for s = 1:numPermutation % loop through shuffles
    switch pairedorunpaired
        case 'paired'
            inverseindexes = repmat((randi(2,1,size(c_all,1))-1.5)'*2,[1,size(c_all,2),size(c_all,3)]);
            realRtmp=c_all.*inverseindexes;
            [~,P]=ttest(realRtmp);
            mean_diff = mean(realRtmp, 1);
            
        case 'unpaired'
            N_total=size(c_all,1);
            [~,six]=sort(randn(1,N_total));
            ix1=six(1:n1);
            ix2=six(n1+1:end);
            [~,P]=ttest2(c_all(ix1,:,:),c_all(ix2,:,:));
            mean_diff = mean(c_all(ix1,:,:), 1) - mean(c_all(ix2,:,:), 1);
    end
    
    
    
    P= squeeze(P);
    
    t_map = squeeze(mean_diff ./ (stderr_diff + eps));
    
    pos_thresh = (P < pthreshold) & t_map>0 ;
    neg_thresh = (P < pthreshold) & t_map<0 ;
    
    % Find clusters using connected components
    pos_clusters = bwconncomp(pos_thresh); % i would prefer 4 here
    neg_clusters = bwconncomp(neg_thresh); % i would prefer 4 here
    
    % max sum of the tmap
    sumTpos = cellfun(@(x) sum(t_map(x)),pos_clusters.PixelIdxList);
    sumTneg = cellfun(@(x) sum(t_map(x)),neg_clusters.PixelIdxList);
    T_pos(s) = max([sumTpos 0]);
    T_neg(s) = min([sumTneg 0]);
end


%% now find significant clusters in the real data


switch pairedorunpaired
    case 'paired'
        [~,P]=ttest(c_all);
        mean_diff=mean(c1-c2,1);
        
    case 'unpaired'
        [~,P]=ttest2(c1,c2);
        mean_diff=mean(c1,1)-mean(c2,1);
end

P= squeeze(P);
t_map = squeeze(mean_diff ./ (stderr_diff + eps));

pos_thresh = (P < pthreshold) & t_map>0 ;
neg_thresh = (P < pthreshold) & t_map<0 ;


% Find clusters using connected components
pos_clusters = bwconncomp(pos_thresh); % i would prefer 4 here
neg_clusters = bwconncomp(neg_thresh); % i would prefer 4 here

% max sum of the tmap
sumTpos = cellfun(@(x) sum(t_map(x)),pos_clusters.PixelIdxList);
sumTneg = cellfun(@(x) sum(t_map(x)),neg_clusters.PixelIdxList);


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

significance=zeros(size(t_map));
significance(sigpixpos)=1;
significance(sigpixneg)=-1;

%% store T (p-value)
significance=significance';
end
