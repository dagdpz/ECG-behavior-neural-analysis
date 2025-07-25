
close all,
clear all
clc

%%
monkey ='Magnus';%'Bacchus';

if strcmp(monkey,'Bacchus')
    path_to_save = 'Y:\Projects\Pulv_bodysignal\LFP\ECG_Bacchus_TaskRest_generalTrig_finalRun\not filtered\grand_average_wo combine_hemispheres_LFP real_New GrandAvg version\';
    load('Y:\Projects\Pulv_bodysignal\LFP\ECG_Bacchus_TaskRest_generalTrig_finalRun\not filtered\grand_average_wo combine_hemispheres_LFP real_New GrandAvg version\Bacchus_Rpeak_Triggered_target_wise_Grand_grand_avg_sessions_sitesall.mat')
elseif strcmp(monkey,'Magnus')
    path_to_save = 'Y:\Projects\Pulv_bodysignal\LFP\ECG_Magnus_TaskRest_generalTrig_finalRun\not filtered\grand_average_wo combine_hemispheres_LFP real_New GrandAvg version\';
    load('Y:\Projects\Pulv_bodysignal\LFP\ECG_Magnus_TaskRest_generalTrig_finalRun\not filtered\grand_average_wo combine_hemispheres_LFP real_New GrandAvg version\Magnus_Rpeak_Triggered_target_wise_Grand_grand_avg_sessions_sitesall.mat')
end

%%
targets = {'VPL', 'dPul', 'MD'};% unique({sites.target});
cond = {'Rest','Task'};
%%
cfg.monkey = monkey;
cfg.lfp.foi = logspace(log10(2),log10(120),60); %
cfg.lfp.smoothWin = 5;
cfg.lfp.frequency_bands = [2 4; 4 8; 8 14; 14 30; 30 50; 70 120];

% Smoothing Kernel here:
win = 1:cfg.lfp.smoothWin; win=win-(numel(win)+1)/2;
half_win = ceil(size(win,2)/2)-1;
gaussian_kernel=normpdf(win,0,numel(win)/6);
gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);

% cfg.analyse_states = {'R',    'Rpeak',1,-0.25, 0.25;...
%                       'R_low', 'Rpeak_lowIBI',1,-0.25, 0.25;...
%                       'R_high', 'Rpeak_highIBI',1,-0.25, 0.25;...
%                       'Cue',  'state',4,-0.10, 0.4};
cfg.analyse_states = {'R',    'Rpeak',1,-0.25, 0.25};
%% Computing the Target-wise averaging of total available sites

% potential cutoff -> probably startfreq should go into settings at somepoint
startfreq=3.5;
colsbp=jet(length(cfg.lfp.frequency_bands));
fbx=find(cfg.lfp.frequency_bands(:,1)>=startfreq,1,'first');
frequency_bands=cfg.lfp.frequency_bands(fbx:end,:);
%freqName = cfg.lfp.freqName(fbx:end);
colsbp=colsbp(fbx:end,:);
freqb=num2cell(strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz'),2);
fx=find(cfg.lfp.foi>=startfreq,1,'first');
freq = cfg.lfp.foi(fx:end);

%%
for t = 1: length(targets)
    sites_for_this_target=arrayfun(@(x) any(strfind(x.target,targets{t})),sites);
    target_sites = sites(sites_for_this_target);
    condition=vertcat(target_sites.condition);
    
    for c = 1:size(condition,2)
        events=cat(1,condition(:,c).event);
        
        
        for e=1:size(cfg.analyse_states,1)
            E              = [];
            E.pow          = [];
            E.itpc         = [];
            E.lfp          = [];
            E.itpcbp       = [];
            E.powbp        = [];
            E.tfr_time     = [];
            E.lfp_time     = [];
            E.lfp_std      = [];
            
            E.nTriggers = mean([events(:,e).nTriggers]);
            E.pow       = (cat(3,events(:,e).pow));
            E.itpc      = (cat(3,events(:,e).itpc));
            E.lfp       = (cat(3,events(:,e).lfp));
            
            E.itpcbp    = (cat(3,events(:,e).itpcbp));
            E.powbp     = (cat(3,events(:,e).powbp));
            E.tfr_time  = events(1,e).tfr_time;
            E.lfp_time  = events(1,e).time;
            E.lfp_std   = std(cat(3,events(:,e).lfp),0,3);
            
            E.pow       = E.pow(fx:end,:,:);
            E.itpc      = E.itpc(fx:end,:,:);
            E.itpcbp    = E.itpcbp(fbx:end,:,:);
            E.powbp     = E.powbp(fbx:end,:,:);
            
            
            eventname = cfg.analyse_states{e,1};
            tar(t).con(c).(eventname)=E;
        end
        tar(t).con(c).cond_name = cond{c};
    end
    tar(t).target = targets{t};
    tar(t).nSites = length(target_sites);
end
%%
numPermutation = 1000;
pthreshold = 0.05;
clusterthreshold = 0.05;
% convert p threshold into a t threshold
tThreshold = (1-pthreshold)/2;
Fin = {'itpc','pow'};
tin = {[1,2],[1,3],[2,3]}; %1: length(targets);
c = 1; % comparing nuclei only in Rest
% %% Method 1
% for Tin = 1:length(tin)
%     t = tin{Tin};
%     for fin = 1:length(Fin)
%         for e = 1:size(cfg.analyse_states,1)
%             
%             data_targ1 = tar(t(1)).con(c).(eventname).(Fin{fin}); % grand adverage of MD target, Rest Condition , across grand_avg(2).nSites;
%             data_targ2 = tar(t(2)).con(c).(eventname).(Fin{fin}); % grand adverage of MD target, Task Condition ,across grand_avg(2).nSites;
%             tfr_freq = freq;
%             tfr_time = tar(1).con(c).(eventname).tfr_time;
%             
%             
%             targ1_R = permute(data_targ1, [3, 2, 1]);  % [sites x  time x frequencies]
%             targ2_R = permute(data_targ2, [3, 2, 1]);  % [sites x  time x frequencies]
%             
%             tic
%             % potentially you have to loop this through all time frequency bins
%             % assumeing dimesnions are (site X time X frequency)
%             timebins = length(tfr_time);
%             freqbins = length(tfr_freq);
%             for time = 1:timebins
%                 for freqi = 1:freqbins
%                     [H,P,CI,STATS]=ttest2(targ1_R(:,time,freqi),targ2_R(:,time,freqi));
%                     T_R(time,freqi)=STATS.tstat;
%                 end
%             end
%             toc
%             
%             %T=T/(size(dat,1)-1)-0.5;
%             sigpos_R = T_R>tThreshold;
%             signeg_R = T_R<-tThreshold;
%             
%             concatTarg = cat(1,targ2_R,targ1_R);
%             nsite1 = size(targ1_R,1);
%             nsite2 = size(targ2_R,1);
%             targ1_shf = zeros([numPermutation,nsite1,size(targ1_R,2),size(targ1_R,3)]);
%             targ2_shf = zeros([numPermutation,nsite2,size(targ2_R,2),size(targ2_R,3)]);
%             tic
%             for s = 1:numPermutation % loop through shuffles
%                 idx1 = randsample(nsite1+nsite2,nsite1);
%                 idx2 = setdiff(1:(nsite1+nsite2),idx1);
%                 targ1_shf(s,:,:,:) = squeeze(concatTarg(idx1,:,:));
%                 targ2_shf(s,:,:,:) = squeeze(concatTarg(idx2,:,:));
%             end
%             toc
%             
%             tic
%             for s = 1:numPermutation % loop through shuffles
%                 targ1 = squeeze(targ1_shf(s,:,:,:));
%                 targ2 = squeeze(targ2_shf(s,:,:,:));
%                 
%                 for time = 1:timebins
%                     for freqi = 1:freqbins
%                         [H,P,CI,STATS]=ttest2(targ1(:,time,freqi),targ2(:,time,freqi));
%                         T(time,freqi)=STATS.tstat;
%                     end
%                 end
%                 
%                 sigpos=T>tThreshold;
%                 signeg=T<-tThreshold;
%                 
%                 CCpos = bwconncomp(sigpos);
%                 CCneg = bwconncomp(signeg);
%                 
%                 % this max sum
%                 sumTpos = cellfun(@(x) sum(T(x)),CCpos.PixelIdxList);
%                 sumTneg = cellfun(@(x) sum(T(x)),CCneg.PixelIdxList);
%                 
%                 T_pos(s) = max([sumTpos 0]);
%                 T_neg(s) = min([sumTneg 0]);
%                 
%             end
%             toc
%             
%             %% now find significant clusters in the real data
%             T_max_thr=prctile(T_pos,100*(1-clusterthreshold));   % sum of
%             T_min_thr=prctile(T_neg,100*clusterthreshold);    % sum of
%             
%             
%             CCpos = bwconncomp(sigpos_R);
%             CCneg = bwconncomp(signeg_R);
%             sigCpos = cellfun(@(x) sum(T_R(x))>=T_max_thr,CCpos.PixelIdxList);
%             sigCneg = cellfun(@(x) sum(T_R(x))<=T_min_thr,CCneg.PixelIdxList);
%             sigpixpos=vertcat(CCpos.PixelIdxList{sigCpos});
%             sigpixneg=vertcat(CCneg.PixelIdxList{sigCneg});
%             significance_pos=zeros(size(T));
%             significance_neg=zeros(size(T));
%             significance_pos(sigpixpos)=true;
%             significance_neg(sigpixneg)=true;
%             
%             % ==============================================================
%             % plotting the significance
%             h(e) = figure;
%             h(e).WindowState = 'maximized';
%             
%             % Plotting the Task data
%             subplot(131)
%             imagesc(tfr_time, 1:numel(tfr_freq),(mean(data_targ2,3)));
%             set(gca,'YDir','normal');
%             %     set(gca, 'YScale', 'log');  % If you want to show the y-axis in logarithmic scale
%             fbandstart = unique(frequency_bands(:))';
%             set(gca,'Xlim',[-.25 .25]);
%             line([0 0], ylim, 'color', 'k');
%             % horizontal lines to separate frequency bands
%             fbandstart_idx = zeros(size(fbandstart));
%             for f = fbandstart
%                 f_idx = find(abs(tfr_freq - f) == min(abs(tfr_freq - f)), 1, 'first');
%                 line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
%                 fbandstart_idx(fbandstart == f) = f_idx;
%             end
%             set(gca,'TickDir','out')
%             set(gca, 'ytick', fbandstart_idx);
%             set(gca, 'yticklabel', fbandstart);
%             set(gca, 'ylim', [0.5,numel(tfr_freq) + 0.5]);
%             axis square;
%             colormap(jet);  % Change the colormap if desired
%             colorbar;  % Show colorbar
%             title([targets{t(1)},' - ',Fin{fin},' - ',cond{c},' - nSite: ', num2str(size(data_targ2,3))],'fontsize', 8,'Interpreter', 'none');
%             % ==============================================================
%             % Plotting the Rest data
%             subplot(132)
%             imagesc(tfr_time, 1:numel(tfr_freq),(mean(data_targ1,3)));
%             set(gca,'YDir','normal');
%             %     set(gca, 'YScale', 'log');  % If you want to show the y-axis in logarithmic scale
%             fbandstart = unique(frequency_bands(:))';
%             set(gca,'Xlim',[-.25 .25]);
%             line([0 0], ylim, 'color', 'k');
%             % horizontal lines to separate frequency bands
%             fbandstart_idx = zeros(size(fbandstart));
%             for f = fbandstart
%                 f_idx = find(abs(tfr_freq - f) == min(abs(tfr_freq - f)), 1, 'first');
%                 line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
%                 fbandstart_idx(fbandstart == f) = f_idx;
%             end
%             set(gca,'TickDir','out')
%             set(gca, 'ytick', fbandstart_idx);
%             set(gca, 'yticklabel', fbandstart);
%             set(gca, 'ylim', [0.5,numel(tfr_freq) + 0.5]);
%             axis square;
%             colormap(jet);  % Change the colormap if desired
%             colorbar;  % Show colorbar
%             title([targets{t(2)},' - ',Fin{fin},' - ',cond{c},' - nSite: ', num2str(size(data_targ1,3))],'fontsize', 8,'Interpreter', 'none');
%             
%             % ==============================================================
%             % Plotting the Stats data
%             subplot(133)
%             imagesc(tfr_time, 1:numel(tfr_freq),(mean(data_targ2,3)- mean(data_targ1,3)));
%             %         imagesc(stat_tfs.time, 1:numel(stat_tfs.freq),stat_tfs.stat);
%             set(gca,'YDir','normal');
%             %     set(gca, 'YScale', 'log');  % If you want to show the y-axis in logarithmic scale
%             fbandstart = unique(frequency_bands(:))';
%             set(gca,'Xlim',[-.25 .25]);
%             line([0 0], ylim, 'color', 'k');
%             % horizontal lines to separate frequency bands
%             fbandstart_idx = zeros(size(fbandstart));
%             for f = fbandstart
%                 f_idx = find(abs(tfr_freq - f) == min(abs(tfr_freq - f)), 1, 'first');
%                 line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
%                 fbandstart_idx(fbandstart == f) = f_idx;
%             end
%             set(gca,'TickDir','out')
%             set(gca, 'ytick', fbandstart_idx);
%             set(gca, 'yticklabel', fbandstart);
%             set(gca, 'ylim', [0.5,numel(tfr_freq) + 0.5]);
%             axis square;
%             colormap(jet);  % Change the colormap if desired
%             colorbar;  % Show colorbar
%             title([Fin{fin},' - ',targets{t(1)},' vs. ',targets{t(2)},' difference - with Sig. Clusters' ],'fontsize', 8,'Interpreter', 'none');
%             
%             hold on
%             %         contour(stat_tfs.time, 1:numel(stat_tfs.freq), stat_tfs.mask, [0.5, 0.5], 'LineColor', 'k', 'LineWidth', 2); % Overlay significant clusters
%             contour(tfr_time, 1:numel(tfr_freq), significance_pos',  1.5 , 'k'); % Overlay significant clusters
%             contour(tfr_time, 1:numel(tfr_freq), significance_neg',  1.5 , 'k'); % Overlay significant clusters
%             data_diff = (mean(data_targ2,3)- mean(data_targ1,3));
%             set(gca,'clim', [min(min(data_diff)) max(max(data_diff))])
%             
%             sgtitle([monkey,'-',Fin{fin},'-',targets{t(1)},' vs. ',targets{t(2)},'- Significant difference in: ',cond{c},'- numPerm =',...
%                 num2str(numPermutation)],'fontsize', 12,'Interpreter', 'none');
%             results_file = fullfile([path_to_save,filesep,monkey,'-',targets{t(1)},' vs. ',targets{t(2)},'-',Fin{fin},'_Significant difference']);
%             export_fig(h(e),[results_file,'.pdf']);
%         end
%     end
% end

%% Method 2
for Tin = 1:length(tin)
    t = tin{Tin};
    for fin = 1:length(Fin)
        for e = 1:size(cfg.analyse_states,1)
            
            data_targ1 = tar(t(1)).con(c).(eventname).(Fin{fin}); % grand adverage of MD target, Rest Condition , across grand_avg(2).nSites;
            data_targ2 = tar(t(2)).con(c).(eventname).(Fin{fin}); % grand adverage of MD target, Task Condition ,across grand_avg(2).nSites;
            
            
            
            tfr_freq = freq;
            tfr_time = tar(1).con(c).(eventname).tfr_time;
            
            
            targ1_R = permute(data_targ1, [3, 2, 1]);  % [sites x  time x frequencies]
            targ2_R = permute(data_targ2, [3, 2, 1]);  % [sites x  time x frequencies]
            
                nonnan1=mean(targ1_R,1);nonnan1(isnan(nonnan1))=[];
                nonnan2=mean(targ2_R,1);nonnan2(isnan(nonnan2))=[];
                collim=max(abs([min(nonnan1(:)) max(nonnan1(:)) min(nonnan2(:)) max(nonnan2(:))]));
           
            
            concatTarg = cat(1,targ2_R,targ1_R);
            targ1_shf = zeros([numPermutation,size(targ1_R,2),size(targ1_R,3)]);
            targ2_shf = zeros([numPermutation,size(targ2_R,2),size(targ2_R,3)]);
            
            % generating the shuffled files :
            tic
            for s = 1:numPermutation % loop through shuffles
                nsite1 = size(targ1_R,1);
                nsite2 = size(targ2_R,1);
                idx1 = randsample(nsite1+nsite2,nsite1);
                idx2 = setdiff(1:(nsite1+nsite2),idx1);
                targ1_shf(s,:,:) = squeeze(mean(concatTarg(idx1,:,:),1));
                targ2_shf(s,:,:) = squeeze(mean(concatTarg(idx2,:,:),1));
            end
            toc
         
            tic
            shf = targ1_shf - targ2_shf;
            s_ix = 1:numPermutation;
            for s = 1:numPermutation % loop through shuffles
                
                T=squeeze(sum(shf(s_ix~=s,:,:)<repmat(shf(s,:,:),[size(shf,1)-1,1,1]),1));
                T=T/(size(shf,1)-1)-0.5;
                
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
            toc
            
            
            %% now find significant clusters in the real data
            T_max_thr=prctile(T_pos,100*(1-clusterthreshold));   % sum of
            T_min_thr=prctile(T_neg,100*clusterthreshold);    % sum of
            
            T=squeeze(sum(shf<repmat((mean(targ1_R,1)-mean(targ2_R,1)),[size(shf,1),1,1]),1));
            T=T/(size(shf,1)-1)-0.5;
            
            sigpos=T>tThreshold;
            signeg=T<-tThreshold;
            
            CCpos = bwconncomp(sigpos);
            CCneg = bwconncomp(signeg);
            sigCpos = cellfun(@(x) sum(T(x))>=T_max_thr,CCpos.PixelIdxList);
            sigCneg = cellfun(@(x) sum(T(x))<=T_min_thr,CCneg.PixelIdxList);
            sigpixpos=vertcat(CCpos.PixelIdxList{sigCpos});
            sigpixneg=vertcat(CCneg.PixelIdxList{sigCneg});
            significance_pos=zeros(size(T));
            significance_neg=zeros(size(T));
            significance_pos(sigpixpos)=true;
            significance_neg(sigpixneg)=true;
            
            
            % ==============================================================
            % plotting the significance
            h(e) = figure;
            h(e).WindowState = 'maximized';
            
            % Plotting the Task data
            subplot(131)
            imagesc(tfr_time, 1:numel(tfr_freq),(mean(data_targ2,3)));
            set(gca,'YDir','normal');
            %     set(gca, 'YScale', 'log');  % If you want to show the y-axis in logarithmic scale
            fbandstart = unique(frequency_bands(:))';
            set(gca,'Xlim',[-.25 .25]);
            line([0 0], ylim, 'color', 'k');
            % horizontal lines to separate frequency bands
            fbandstart_idx = zeros(size(fbandstart));
            for f = fbandstart
                f_idx = find(abs(tfr_freq - f) == min(abs(tfr_freq - f)), 1, 'first');
                line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                fbandstart_idx(fbandstart == f) = f_idx;
            end
            set(gca,'TickDir','out')
            set(gca, 'ytick', fbandstart_idx);
            set(gca, 'yticklabel', fbandstart);
            set(gca, 'ylim', [0.5,numel(tfr_freq) + 0.5]);            
            set(gca,'CLim',[-collim collim]);
            axis square;
            colormap(jet);  % Change the colormap if desired
            colorbar;  % Show colorbar
            title([targets{t(1)},' - ',Fin{fin},' - ',cond{c},' - nSite: ', num2str(size(data_targ2,3))],'fontsize', 8,'Interpreter', 'none');
            % ==============================================================
            % Plotting the Rest data
            subplot(132)
            imagesc(tfr_time, 1:numel(tfr_freq),(mean(data_targ1,3)));
            set(gca,'YDir','normal');
            %     set(gca, 'YScale', 'log');  % If you want to show the y-axis in logarithmic scale
            fbandstart = unique(frequency_bands(:))';
            set(gca,'Xlim',[-.25 .25]);
            line([0 0], ylim, 'color', 'k');
            % horizontal lines to separate frequency bands
            fbandstart_idx = zeros(size(fbandstart));
            for f = fbandstart
                f_idx = find(abs(tfr_freq - f) == min(abs(tfr_freq - f)), 1, 'first');
                line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                fbandstart_idx(fbandstart == f) = f_idx;
            end
            set(gca,'TickDir','out')
            set(gca, 'ytick', fbandstart_idx);
            set(gca, 'yticklabel', fbandstart);
            set(gca, 'ylim', [0.5,numel(tfr_freq) + 0.5]);            
            set(gca,'CLim',[-collim collim]);
            axis square;
            colormap(jet);  % Change the colormap if desired
            colorbar;  % Show colorbar
            title([targets{t(2)},' - ',Fin{fin},' - ',cond{c},' - nSite: ', num2str(size(data_targ1,3))],'fontsize', 8,'Interpreter', 'none');
            
            % ==============================================================
            % Plotting the Stats data
            subplot(133)
            imagesc(tfr_time, 1:numel(tfr_freq),(mean(data_targ2,3)- mean(data_targ1,3)));
            %         imagesc(stat_tfs.time, 1:numel(stat_tfs.freq),stat_tfs.stat);
            set(gca,'YDir','normal');
            %     set(gca, 'YScale', 'log');  % If you want to show the y-axis in logarithmic scale
            fbandstart = unique(frequency_bands(:))';
            set(gca,'Xlim',[-.25 .25]);
            line([0 0], ylim, 'color', 'k');
            % horizontal lines to separate frequency bands
            fbandstart_idx = zeros(size(fbandstart));
            for f = fbandstart
                f_idx = find(abs(tfr_freq - f) == min(abs(tfr_freq - f)), 1, 'first');
                line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                fbandstart_idx(fbandstart == f) = f_idx;
            end
            set(gca,'TickDir','out')
            set(gca, 'ytick', fbandstart_idx);
            set(gca, 'yticklabel', fbandstart);
            set(gca, 'ylim', [0.5,numel(tfr_freq) + 0.5]); 
            
            nonnan2=mean(data_targ2,3)- mean(data_targ1,3);nonnan2(isnan(nonnan2))=[];
            collim=max(abs([min(nonnan2(:)) max(nonnan2(:))]));
            set(gca,'CLim',[-collim collim]);
            axis square;
            colormap(jet);  % Change the colormap if desired
            colorbar;  % Show colorbar
            title([Fin{fin},' - ',targets{t(1)},' vs. ',targets{t(2)},' difference - with Sig. Clusters' ],'fontsize', 8,'Interpreter', 'none');
            
            hold on
            %         contour(stat_tfs.time, 1:numel(stat_tfs.freq), stat_tfs.mask, [0.5, 0.5], 'LineColor', 'k', 'LineWidth', 2); % Overlay significant clusters
            contour(tfr_time, 1:numel(tfr_freq), significance_pos',  1.5 , 'k'); % Overlay significant clusters
            contour(tfr_time, 1:numel(tfr_freq), significance_neg',  1.5 , 'k'); % Overlay significant clusters
            data_diff = (mean(data_targ2,3)- mean(data_targ1,3));
            set(gca,'clim', [min(min(data_diff)) max(max(data_diff))])
            
            sgtitle([monkey,'-',Fin{fin},'-',targets{t(1)},' vs. ',targets{t(2)},'- Significant difference in: ',cond{c},'- numPerm =',...
                num2str(numPermutation)],'fontsize', 12,'Interpreter', 'none');
            results_file = fullfile([path_to_save,filesep,monkey,'-',targets{t(1)},' vs. ',targets{t(2)},'-',Fin{fin},'_Significant difference_mthd2']);
            export_fig(h(e),[results_file,'.pdf']);
        end
    end
end