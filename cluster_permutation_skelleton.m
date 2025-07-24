 
% convert p threshold into a t threshold
tThreshold = (1-pthreshold)/2;   

%% potentially you have to loop this through all time frequency bins
% assumeing dimesnions are (site X time X frequency)

for time=timebins
    for freq=freqbins
        [H,P,CI,STATS]=ttest(task_R(:,time,freq),rest_R(:,time,freq));
        %[H,P,CI,STATS]=ttest2(task_R(:,time,freq),rest_R(:,time,freq));
        T_R(time,freq)=STATS.tstat;
    end
end

%T=T/(size(dat,1)-1)-0.5;
sigpos_R=T>tThreshold;
signeg_R=T<-tThreshold;

% find (for each surrogate) the max cluster significance (somehow i
% feel this should be done PER frequency, hmmmm....)

taskandrest=cat(1,task_R,rest_R);
for s=s_ix % loop through shuffles
    
%     nsitest=size(task,1); %% check if really first dimension    
%     nsitesr=size(rest,1); %% check if really first dimension   
%     idxt=randsample(nsitest+nsitesr,nsitest);
%     idxr=setdiff(1:(nsitest+nsitesr),idxt);  
%     task=taskandrest(idxt,:,:);
%     rest=taskandrest(idxr,:,:);

    inverseindexes=(randi(2,1,size(task_R,1))-1.5)*2;    
    task=task_R*inverseindexes;
    rest=rest_R*inverseindexes;  
    
    for time=timebins
        for freq=freqbins
            [H,P,CI,STATS]=ttest(task(:,time,freq),rest(:,time,freq));
            %[H,P,CI,STATS]=ttest2(task(:,time,freq),rest(:,time,freq));
            T(time,freq)=STATS.tstat;
        end
    end
        
    sigpos=T>tThreshold;
    signeg=T<-tThreshold;
    
    CCpos = bwconncomp(sigpos);
    CCneg = bwconncomp(signeg);
    
    % this max sum
    sumTpos = cellfun(@(x) sum(T(x)),CCpos.PixelIdxList);
    sumTneg = cellfun(@(x) sum(T(x)),CCneg.PixelIdxList);
    
    T_pos(s) = max([sumTpos 0]);
    T_neg(s) = min([sumTneg 0]);
end




%% now find significant clusters in the real data
T_max_thr=prctile(T_pos,100*(1-clusterthreshold));   % sum of
T_min_thr=prctile(T_neg,100*clusterthreshold);    % sum of


CCpos = bwconncomp(sigpos_R);
CCneg = bwconncomp(signeg_R);
sigCpos = cellfun(@(x) sum(T_R(x))>=T_max_thr,CCpos.PixelIdxList);
sigCneg = cellfun(@(x) sum(T_R(x))<=T_min_thr,CCneg.PixelIdxList);
sigpixpos=vertcat(CCpos.PixelIdxList{sigCpos});
sigpixneg=vertcat(CCneg.PixelIdxList{sigCneg});

        
        
        
        
%         significance.(fn)=zeros(size(realD.(fn).mean));
%         significance.(fn)(sigpixpos)=true;
%         significance.(fn)(sigpixneg)=true;


% 
% 
%         
%         %% now find significant clusters in the real data        
%         T_max_thr=prctile(T_pos,100*(1-clusterthreshold));   % sum of     
%         T_min_thr=prctile(T_neg,100*clusterthreshold);    % sum of    
%         est_mean=squeeze(mean(dat,1));
%         est_var=squeeze(var(dat,1));
%         
%         %T=(squeeze(realD.(fn).mean)-est_mean)./sqrt(est_var); 
%         
%         T=squeeze(sum(dat<repmat(realD.(fn).mean,[size(dat,1),1,1]),1));
%         T=T/(size(dat,1)-1)-0.5;
%         
%         sigpos=T>tThreshold;        
%         signeg=T<-tThreshold;
      