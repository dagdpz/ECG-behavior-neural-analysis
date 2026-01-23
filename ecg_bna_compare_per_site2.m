function [ out] = ecg_bna_compare_per_site( trig, data_label, cfg, ncomp )
% Assuming the data is loaded in the following variables:
% data_task: data for the task condition (freq x time x trials)
% data_rest: data for the rest condition (freq x time x trials)
% tfr_time: A vector representing the time points (1 x time)
% tfr_freq: A vector representing the frequency points (1 x freq)
%path_to_trig = fullfile(cfg.fldr.LFP_sites);

c1=cfg.lfp.compare_conditions{ncomp}(1);
c2=cfg.lfp.compare_conditions{ncomp}(2);
% events = {trig.condition(1).event.event_name};
for e = 1: length(cfg.events)
    
    out.event(e).time=trig.condition(c1).event(e).time;
    out.event(e).tfr_time=trig.condition(c1).event(e).tfr_time;
    out.event(e).event_name=trig.condition(c1).event(e).event_name;
    out.event(e).observed.ntriggers=trig.condition(c1).event(e).observed.ntriggers;
    out.event(e).surrogate.ntriggers=trig.condition(c2).event(e).observed.ntriggers;
    for fi = 1:length(data_label)
        Fin = data_label{fi};
        data1 = trig.condition(c1).event(e).observed.(Fin).complete;
        data2 = trig.condition(c2).event(e).observed.(Fin).complete;
        if ismember(Fin,{'powbp','itpcbp'}) % treat each frequwncy independantly here!
            cluster_t_values=[];
            significance=false(1,size(data1,2),size(data1,3));
            surrogate.mean=zeros(size(significance));
            surrogate.std=zeros(size(significance));
            surrogate.sterr=zeros(size(significance));
            surrogate.conf95=zeros(2,size(data1,2),size(data1,3));
            observed.mean=zeros(size(significance));
            observed.std=zeros(size(significance));
            observed.sterr=zeros(size(significance));
            for f=1:size(data1,2)
                [signtmp, cluster_t_tmp, obstmp, surrtmp]=do_test2(data1(:,f,:),data2(:,f,:),Fin);
                significance(:,f,:)=signtmp;
                observed.mean(:,f,:)=obstmp.mean;
                observed.std(:,f,:)=obstmp.std;
                observed.sterr(:,f,:)=obstmp.sterr;
                surrogate.mean(:,f,:)=surrtmp.mean;
                surrogate.std(:,f,:)=surrtmp.std;
                surrogate.sterr(:,f,:)=surrtmp.sterr;
                surrogate.conf95(:,f,:)=surrtmp.conf95;
                cluster_t_values=[cluster_t_values {cluster_t_tmp}];
            end
        else
            [significance, cluster_t_values, observed, surrogate]=do_test2(data1,data2,Fin);
        end
        
        %% store this shit
        out.event(e).observed.(Fin)=observed;
        out.event(e).surrogate.(Fin)=surrogate;
        out.event(e).significance.(Fin)=significance;
        out.event(e).cluster_t_values.(Fin)=cluster_t_values;
        
    end
end
end

function [significance, cluster_T]=do_test(data1,data2)

numPermutation = 1000;
pthreshold=0.05;
clusterthreshold=0.05;

n1=size(data1,1); % trials x freq x time
n2=size(data2,1);
dataall=[data1;data2];
df = n1+n2-2;

    t_thresh_pos = tinv(1-pthreshold/2, df);
    t_thresh_neg = tinv(pthreshold/2, df);
    
T_pos=zeros(1,numPermutation);
T_neg=zeros(1,numPermutation);

[~,RX]=sort(rand(size(dataall,1),numPermutation),1);

for s = 1:numPermutation % loop through shuffles
    d1=dataall(RX(1:n1,s),:,:);
    d2=dataall(RX(n1+1:end,s),:,:);
    
    T_map = simpleTTest2(d1,d2,1);
    %[~,P,~,stats] = ttest2(d1,d2);
    
%     pos_thresh = (P < pthreshold) & T_map>0 ;
%     neg_thresh = (P < pthreshold) & T_map<0 ;
    pos_thresh = T_map>t_thresh_pos ;
    neg_thresh = T_map<t_thresh_neg ;
    
    pos_clusters = bwconncomp(pos_thresh);
    neg_clusters = bwconncomp(neg_thresh);
    % max sum of the tmap
    sumTpos = cellfun(@(x) sum(T_map(x)),pos_clusters.PixelIdxList);
    sumTneg = cellfun(@(x) sum(T_map(x)),neg_clusters.PixelIdxList);
    T_pos(s) = max([sumTpos 0]);
    T_neg(s) = min([sumTneg 0]);
end

[~,P,~,stats] = ttest2(data1,data2);
T_map=stats.tstat;
pos_thresh = (P < pthreshold) & T_map>0 ;
neg_thresh = (P < pthreshold) & T_map<0 ;


% Find clusters using connected components
pos_clusters = bwconncomp(pos_thresh); % i would prefer 4 here
neg_clusters = bwconncomp(neg_thresh); % i would prefer 4 here

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

function [significance, cluster_T, out,est]=do_test2(data1,data2,fn)
%tic
numPermutation = 100;
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
RX=[(1:(n1+n2))', RX];

sur=zeros(numPermutation+1,size(dataall,2),size(dataall,3));
% toc
% a=1;
% tic
for s = 1:numPermutation+1 % loop through shuffles - first is observed    
    d1=dataall(RX(1:n1,s),:,:);
    d2=dataall(RX(n1+1:end,s),:,:);    
    if ismember(fn,{'itpc','itpcbp'})
        BB=single(abs(mean(d1,1)))-single(abs(mean(d2,1)));
    elseif ismember(fn,{'pha'})
        BB=single(angle(mean(d1,1)))-single(angle(mean(d2,1)));
    else
       BB=single(mean(d1,1))-single(mean(d2,1));
    end
    sur(s,:,:)=squeeze(BB);
end
obs=sur(1,:,:);
sur(1,:,:)=[];
% toc
% a=1;
% tic
s_ix=1:numPermutation;
for s = s_ix % loop through shuffles - observed is in the end here so we can use the stored sig and tvals?

    dat_s=squeeze(sur(s,:,:));
    est.mean=squeeze(mean(sur(s_ix~=s,:,:),1));
    est_var=squeeze(var(sur(s_ix~=s,:,:),1));
    T_map=(dat_s-est.mean)./sqrt(est_var); %% we actually care about pos/neg here
    
    conf_map=squeeze(sum(sur(s_ix~=s,:,:)<repmat(sur(s,:,:),[size(sur,1)-1,1,1]),1));
    conf_map=conf_map/(size(sur,1)-1);
    
    sigpos=conf_map>t_thresh_pos;
    signeg=conf_map<t_thresh_neg;
    
    pos_clusters = bwconncomp(sigpos);
    neg_clusters = bwconncomp(signeg);
    sumTpos = cellfun(@(x) sum(T_map(x)),pos_clusters.PixelIdxList);
    sumTneg = cellfun(@(x) sum(T_map(x)),neg_clusters.PixelIdxList);
    T_pos(s) = max([sumTpos 0]);
    T_neg(s) = min([sumTneg 0]);
end
% toc
% 
% a=1;
% 
% tic

out.mean=obs;
out.std=mean([std(data1,0,1);std(data2,0,1)],1);
out.sterr=mean([sterr(data1,1);sterr(data2,1)],1);


est.mean=mean(sur,1);
est.std=std(sur,0,1);
est.sterr=sterr(sur,1);
est.conf95=single(prctile(sur,97.5,1));
est.conf95(2,:,:)=single(prctile(sur,2.5,1));

est_var=squeeze(var(sur,1));
T_map=(squeeze(obs)-squeeze(est.mean))./sqrt(est_var);

conf_map=squeeze(sum(sur<repmat(obs,[size(sur,1),1,1]),1));
conf_map=conf_map/(size(sur,1));

sigpos=conf_map>t_thresh_pos;
signeg=conf_map<t_thresh_neg;

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

significance=zeros(size(sur(1,:,:)));
significance(sigpixpos)=1;
significance(sigpixneg)=-1;

cluster_T=[sumTpos(p_values_pos<clusterthreshold) sumTneg(p_values_neg<clusterthreshold)*-1];
% toc
% a=1;
end

