function [significance, pmap]=ecg_bna_population_significance(real,surr)

%%
% Permutation Based on Lukas function
numPermutation = 1000;
%ntoflip=round(size(C1,3)*0.3);
pthreshold = 0.05;
clusterthreshold = 0.05;
% convert p threshold into a t threshold
%tThreshold = (1-pthreshold)/2;

surr_P = permute(surr, [3, 2, 1]);  % [trials x  time x frequencies]
real_P = permute(real, [3, 2, 1]);  % [trials x  time x frequencies]
diff_data=real_P-surr_P;
df = size(diff_data,1) - 1;
t_thresh_pos = tinv(1-pthreshold/2, df);
t_thresh_neg = tinv(pthreshold/2, df);


T_pos=zeros(1,numPermutation);
T_neg=zeros(1,numPermutation);

%generating the shuffled files :
% real_shf=zeros([numPermutation,size(real_P,2),size(real_P,3)]);
% surr_shf=zeros([numPermutation,size(surr_P,2),size(surr_P,3)]);
for s = 1:numPermutation % loop through shuffles
    inverseindexes = repmat((randi(2,1,size(real_P,1))-1.5)'*2,[1,size(real_P,2),size(real_P,3)]);
    realRtmp=real_P;
    surrRtmp=surr_P;
    realRtmp(inverseindexes==-1)=surr_P(inverseindexes==-1);
    surrRtmp(inverseindexes==-1)=real_P(inverseindexes==-1);
    
    % Compute this permutations t-statistics    
    diff_data=realRtmp-surrRtmp;
    mean_diff = mean(diff_data, 1);
    stderr_diff = std(diff_data, 0, 1) / sqrt(size(diff_data,1));
    t_map = squeeze(mean_diff ./ (stderr_diff + eps));
    
    % Find significant pixels 
    pos_thresh = t_map > t_thresh_pos;
    neg_thresh = t_map < t_thresh_neg;
    
    % Find clusters using connected components
    pos_clusters = bwconncomp(pos_thresh, 8); % i would prefer 4 here
    neg_clusters = bwconncomp(neg_thresh, 8); % i would prefer 4 here
    
    % max sum of the tmap
    sumTpos = cellfun(@(x) sum(t_map(x)),pos_clusters.PixelIdxList);
    sumTneg = cellfun(@(x) sum(t_map(x)),neg_clusters.PixelIdxList);
    T_pos(s) = max([sumTpos 0]);
    T_neg(s) = min([sumTneg 0]);
end

%% now find significant clusters in the real data
diff_data=real_P-surr_P;

% Compute original t-statistics
mean_diff = mean(diff_data, 1);
stderr_diff = std(diff_data, 0, 1) / sqrt(size(diff_data,1));
t_map = squeeze(mean_diff ./ (stderr_diff + eps));
df = size(diff_data,1) - 1;

% Find clusters in original data
t_thresh_pos = tinv(1-clusterthreshold/2, df);
t_thresh_neg = tinv(clusterthreshold/2, df);

pos_thresh = t_map > t_thresh_pos;
neg_thresh = t_map < t_thresh_neg;

pos_clusters = bwconncomp(pos_thresh); % i would prefer 4 here
neg_clusters = bwconncomp(neg_thresh); % i would prefer 4 here

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
% 
% % Create significance map
% significance_map = false(n_time, n_freq);
% for i = 1:n_clusters
%     if p_values(i) < clusterthreshold
%         significance_map(clusters{i}) = true;
%     end
% end



significance=zeros(size(t_map));
significance(sigpixpos)=1;
significance(sigpixneg)=-1;

%% store T (p-value)
pmap=t_map;
significance=significance';    
  
% 
% 
% % weird nonparametric part here
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This file is updated on 23.07.25
% % based on the configurations of the data that was generated for:
% % "grand_average_wo combine_hemispheres_LFP real_New GrandAvg version "
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%
% % Permutation Based on Lukas function
% numPermutation = 1000;
% %ntoflip=round(size(C1,3)*0.3);
% pthreshold = 0.05;
% clusterthreshold = 0.05;
% % convert p threshold into a t threshold
% tThreshold = (1-pthreshold)/2;
% 
% surr_P = permute(surr, [3, 2, 1]);  % [trials x  time x frequencies]
% real_P = permute(real, [3, 2, 1]);  % [trials x  time x frequencies]
% 
% % generating the shuffled files :
% real_shf=zeros([numPermutation,size(real_P,2),size(real_P,3)]);
% surr_shf=zeros([numPermutation,size(surr_P,2),size(surr_P,3)]);
% for s = 1:numPermutation % loop through shuffles
% 
%     realtmp=real_P;
%     surrtmp=surr_P;
%     ntoflip=randsample(0:size(real,3),1);
%     if ntoflip>=1
%     toflip=randsample(1:size(real,3),ntoflip);
%     realtmp(toflip,:,:)=surr_P(toflip,:,:);
%     surrtmp(toflip,:,:)=real_P(toflip,:,:);
%     end
% 
% %     inverseindexes = repmat((randi(2,1,size(task_R,1))-1.5)'*2,[1,size(task_R,2),size(task_R,3)]);
% %     taskRtmp=task_R;
% %     restRtmp=rest_R;
% %     taskRtmp(inverseindexes==-1)=rest_R(inverseindexes==-1);
% %     restRtmp(inverseindexes==-1)=task_R(inverseindexes==-1);
% 
% %     task_shf(s,:,:) = squeeze(mean(task_R(:,:,:).*inverseindexes,1));
% %     rest_shf(s,:,:) = squeeze(mean(rest_R(:,:,:).*inverseindexes,1));
%     real_shf(s,:,:) = squeeze(mean(realtmp,1));
%     surr_shf(s,:,:) = squeeze(mean(surrtmp,1));
% end
% shf=real_shf-surr_shf;
% 
% % %% entirely different approach: flip all together!
% % diff_data=real_P-surr_P;
% % sizediff=size(diff_data);
% % %rand_signs = 2 * (rand(sizediff(1),sizediff(2),sizediff(3), numPermutation) > 0.5) - 1;
% % rand_signs = 2 * (rand(sizediff(1), numPermutation) > 0.5) - 1;
% % rand_signs = repmat(rand_signs,[1,1,sizediff(2),sizediff(3)]);
% % rand_signs =permute(rand_signs,[1,3,4,2]);
% % % Initialize arrays for null distribution
% % T_pos = zeros(1, numPermutation);
% % T_neg = zeros(1, numPermutation);
% %
% % % Reshape diff_data for broadcasting
% % diff_data_expanded = repmat(diff_data, [1, 1, 1, numPermutation]);
% %
% % % Perform sign flipping for all permutations at once
% % shf = permute(squeeze(mean(diff_data_expanded .* rand_signs,1)), [3, 1, 2]);
% % clear diff_data_expanded
% 
% s_ix=1:numPermutation;
% for s = 1:numPermutation % loop through shuffles
%     T=squeeze(sum(shf(s_ix~=s,:,:)<repmat(shf(s,:,:),[size(shf,1)-1,1,1]),1));
%     T=T/(size(shf,1)-1)-0.5;
% 
%     sigpos=T>tThreshold;
%     signeg=T<-tThreshold;
%     CCpos = bwconncomp(sigpos);
%     CCneg = bwconncomp(signeg);
%     % this max sum
%     sumTpos = cellfun(@(x) sum(T(x)),CCpos.PixelIdxList);
%     sumTneg = cellfun(@(x) sum(T(x)),CCneg.PixelIdxList);
%     T_pos(s) = max([sumTpos 0]);
%     T_neg(s) = min([sumTneg 0]);
% end
% 
% 
% %% now find significant clusters in the real data
% T_max_thr=prctile(T_pos,100*(1-clusterthreshold));   % sum of
% T_min_thr=prctile(T_neg,100*clusterthreshold);    % sum of
% 
% T=squeeze(sum(shf<repmat(mean((real_P-surr_P),1),[size(shf,1),1,1]),1));
% T=T/size(shf,1)-0.5;
% 
% sigpos=T>tThreshold;
% signeg=T<-tThreshold;
% 
% CCpos = bwconncomp(sigpos);
% CCneg = bwconncomp(signeg);
% sigCpos = cellfun(@(x) sum(T(x))>=T_max_thr,CCpos.PixelIdxList);
% sigCneg = cellfun(@(x) sum(T(x))<=T_min_thr,CCneg.PixelIdxList);
% sigpixpos=vertcat(CCpos.PixelIdxList{sigCpos});
% sigpixneg=vertcat(CCneg.PixelIdxList{sigCneg});
% % significance_pos=zeros(size(T));
% % significance_neg=zeros(size(T));
% 
% 
% significance=zeros(size(T));
% significance(sigpixpos)=1;
% significance(sigpixneg)=-1;
% 
% %% store T (p-value)
% pmap=1-abs(T*2);
% significance=significance';

end
