function [significance, pmap]=ecg_bna_population_significance_normalized_parametric(norm_data)

%%
% Permutation Based on Lukas function
numPermutation = 1000;
%ntoflip=round(size(C1,3)*0.3);
pthreshold = 0.01;
clusterthreshold = 0.01;
% convert p threshold into a t threshold
%tThreshold = (1-pthreshold)/2;






% surr_P = permute(surr, [3, 2, 1]);  % [trials x  time x frequencies]
% diff_data=real_P-surr_P;
%
% df = size(diff_data,1) - 1;
% t_thresh_pos = tinv(1-pthreshold/2, df);
% t_thresh_neg = tinv(pthreshold/2, df);


T_pos=zeros(1,numPermutation);
T_neg=zeros(1,numPermutation);

%generating the shuffled files :
% real_shf=zeros([numPermutation,size(real_P,2),size(real_P,3)]);
% surr_shf=zeros([numPermutation,size(surr_P,2),size(surr_P,3)]);


real_P = permute(norm_data, [3, 2, 1]);  % [trials x  time x frequencies]

for s = 1:numPermutation % loop through shuffles
    inverseindexes = repmat((randi(2,1,size(real_P,1))-1.5)'*2,[1,size(real_P,2),size(real_P,3)]);
    realRtmp=real_P.*inverseindexes;
    [~,P]=ttest(realRtmp);
    P= squeeze(P);
    
    mean_diff = mean(realRtmp, 1);
    stderr_diff = std(realRtmp, 0, 1) / sqrt(size(realRtmp,1));
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

    [~,P]=ttest(real_P);
    P= squeeze(P);
    mean_diff = mean(real_P, 1);
    stderr_diff = std(real_P, 0, 1) / sqrt(size(real_P,1));
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
pmap=t_map;
significance=significance';


end
