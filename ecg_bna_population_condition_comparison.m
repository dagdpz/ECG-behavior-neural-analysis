function ecg_bna_population_condition_comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is updated on 23.07.25
% based on the configurations of the data that was generated for:
% "grand_average_wo combine_hemispheres_LFP real_New GrandAvg version "
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all,
clear all
clc

%%
filtered = true;
monkey ='Magnus';%'Bacchus';

if strcmp(monkey,'Bacchus')
    path_to_save = 'Y:\Projects\Pulv_bodysignal\LFP\ECG_Bacchus_TaskRest_generalTrig_finalRun\not filtered\grand_average_wo combine_hemispheres_LFP real_New GrandAvg version\';
    if filtered
        load('Y:\Projects\Pulv_bodysignal\LFP\ECG_Bacchus_TaskRest_generalTrig_finalRun\4hz filtered\grand_average_wo combine_hemispheres_LFP real_New GrandAvg version\Bacchus_Rpeak_Triggered_target_wise_Grand_grand_avg_sessions_sitesall.mat')
    else
        load('Y:\Projects\Pulv_bodysignal\LFP\ECG_Bacchus_TaskRest_generalTrig_finalRun\not filtered\grand_average_wo combine_hemispheres_LFP real_New GrandAvg version\Bacchus_Rpeak_Triggered_target_wise_Grand_grand_avg_sessions_sitesall.mat')
    end
elseif strcmp(monkey,'Magnus')
    path_to_save = 'Y:\Projects\Pulv_bodysignal\LFP\ECG_Magnus_TaskRest_generalTrig_finalRun\not filtered\grand_average_wo combine_hemispheres_LFP real_New GrandAvg version\';
    if filtered
        load('Y:\Projects\Pulv_bodysignal\LFP\ECG_Magnus_TaskRest_generalTrig_finalRun\4hz filtered\grand_average_wo combine_hemispheres_LFP real_New GrandAvg version\Magnus_Rpeak_Triggered_target_wise_Grand_grand_avg_sessions_sitesall.mat')
    else
        load('Y:\Projects\Pulv_bodysignal\LFP\ECG_Magnus_TaskRest_generalTrig_finalRun\not filtered\grand_average_wo combine_hemispheres_LFP real_New GrandAvg version\Magnus_Rpeak_Triggered_target_wise_Grand_grand_avg_sessions_sitesall.mat')
    end
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
% UPdates on 24.07.25 >>>>>>>>>>>>>>>>>>>>>>>>>>>
% Permutation Based on Lukas function
numPermutation = 1000;
pthreshold = 0.05;
clusterthreshold = 0.05;
% convert p threshold into a t threshold
tThreshold = (1-pthreshold)/2;
%%
Fin = {'itpc','pow'};
for t = 1: length(targets)
    for fin = 1:length(Fin)
        for e=1:size(cfg.analyse_states,1)
            % Assuming the data is loaded in the following variables:
            % itpc_data_task: ITPC data for the task condition (freq x time x trials)
            % itpc_data_rest: ITPC data for the rest condition (freq x time x trials)
            data_rest = tar(t).con(1).(eventname).(Fin{fin}); % grand adverage of MD target, Rest Condition , across grand_avg(2).nSites;
            data_task = tar(t).con(2).(eventname).(Fin{fin}); % grand adverage of MD target, Task Condition ,across grand_avg(2).nSites;
            tfr_freq = freq;
            tfr_time = tar(t).con(1).(eventname).tfr_time;
            
            
            rest_R = permute(data_rest, [3, 2, 1]);  % [trials x  time x frequencies]
            task_R = permute(data_task, [3, 2, 1]);  % [trials x  time x frequencies]
            
            nonnan1=mean(rest_R,1);nonnan1(isnan(nonnan1))=[];
            nonnan2=mean(task_R,1);nonnan2(isnan(nonnan2))=[];
            collim=max(abs([min(nonnan1(:)) max(nonnan1(:)) min(nonnan2(:)) max(nonnan2(:))]));
            
            % generating the shuffled files :
            tic
            task_shf=zeros([numPermutation,size(task_R,2),size(task_R,3)]);
            rest_shf=zeros([numPermutation,size(rest_R,2),size(rest_R,3)]);
            for s = 1:numPermutation % loop through shuffles
                inverseindexes = repmat((randi(2,1,size(task_R,1))-1.5)'*2,[1,size(task_R,2),size(task_R,3)]);
                task_shf(s,:,:) = squeeze(mean((task_R(:,:,:)).*inverseindexes,1));
                rest_shf(s,:,:) = squeeze(mean((rest_R(:,:,:)).*inverseindexes,1));
            end
            toc
            
            tic
            shf=task_shf-rest_shf;
            s_ix=1:numPermutation;
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
            
            T=squeeze(sum(shf<repmat(mean((task_R-rest_R),1),[size(shf,1),1,1]),1));
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
            imagesc(tfr_time, 1:numel(tfr_freq),(mean(data_task,3)));
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
            title([Fin{fin},' - ',cond{2},' - nSite: ', num2str(size(data_task,3))],'fontsize', 8,'Interpreter', 'none');
            % ==============================================================
            % Plotting the Rest data
            subplot(132)
            imagesc(tfr_time, 1:numel(tfr_freq),(mean(data_rest,3)));
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
            title([Fin{fin},' - ',cond{1},' - nSite: ', num2str(size(data_rest,3))],'fontsize', 8,'Interpreter', 'none');
            
            % ==============================================================
            % Plotting the Stats data
            subplot(133)
            imagesc(tfr_time, 1:numel(tfr_freq),(mean(data_task,3)- mean(data_rest,3)));
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
            axis square;
            colormap(jet);  % Change the colormap if desired
            colorbar;  % Show colorbar
            title([Fin{fin},' - ',cond{2},' vs. ',cond{1},' difference - with Sig. Clusters' ],'fontsize', 8,'Interpreter', 'none');
            
            hold on
            %         contour(stat_tfs.time, 1:numel(stat_tfs.freq), stat_tfs.mask, [0.5, 0.5], 'LineColor', 'k', 'LineWidth', 2); % Overlay significant clusters
            contour(tfr_time, 1:numel(tfr_freq), significance_pos',  1.5 , 'k'); % Overlay significant clusters
            contour(tfr_time, 1:numel(tfr_freq), significance_neg',  1.5 , 'k'); % Overlay significant clusters
           
            nonnan2=mean(data_task,3)- mean(data_rest,3);nonnan2(isnan(nonnan2))=[];
            collim=max(abs([min(nonnan2(:)) max(nonnan2(:))]));
            set(gca,'CLim',[-collim collim]);
            
            sgtitle([monkey,'-',targets{t},'-',Fin{fin},'-',cond{2},'-',cond{1},'- Significant difference - numPerm =',...
                num2str(numPermutation)],'fontsize', 12,'Interpreter', 'none');
            results_file = fullfile([path_to_save,filesep,monkey,'-',targets{t},'-',Fin{fin},'_Significant difference']);
            export_fig(h(e),[results_file,'.pdf']);
        end
    end
end
close all,
%%
% adjusted part for LFP evoked or MUA

Fin = {'lfp','mua'};
for t = 1: length(targets)
    for fin = 1%:length(Fin)
        for e=1:size(cfg.analyse_states,1)
            % Assuming the data is loaded in the following variables:
            % itpc_data_task: ITPC data for the task condition (freq x time x trials)
            % itpc_data_rest: ITPC data for the rest condition (freq x time x trials)
            data_rest = squeeze(tar(t).con(1).(eventname).(Fin{fin})); % grand adverage of MD target, Rest Condition , across grand_avg(2).nSites;
            data_task = squeeze(tar(t).con(2).(eventname).(Fin{fin})); % grand adverage of MD target, Task Condition ,across grand_avg(2).nSites;
            tfr_freq = freq;
            tfr_time = tar(t).con(1).(eventname).tfr_time;
            
            
            rest_R = permute(data_rest, [2, 1]);  % [trials x  time ]
            task_R = permute(data_task, [2, 1]);  % [trials x  time ]
            
            
            % generating the shuffled files :
            tic
            task_shf = zeros(numPermutation, size(task_R,2));
            rest_shf = zeros(numPermutation, size(rest_R,2));
            for s = 1:numPermutation
                inverseindexes = (randi(2, size(task_R,1), 1) - 1.5) * 2;  % [trials x 1]
                task_shf(s,:) = mean(task_R .* inverseindexes, 1);
                rest_shf(s,:) = mean(rest_R .* inverseindexes, 1);
            end
            toc
            
            tic
            shf=task_shf-rest_shf;
            s_ix=1:numPermutation;
            for s = 1:numPermutation
                T = sum(shf(s_ix~=s,:) < repmat(shf(s,:), [numPermutation-1,1]), 1) / (numPermutation-1) - 0.5;
                
                sigpos = T > tThreshold;
                signeg = T < -tThreshold;
                
                CCpos = bwconncomp(sigpos);  % 1D cluster
                CCneg = bwconncomp(signeg);
                
                sumTpos = cellfun(@(x) sum(T(x)), CCpos.PixelIdxList);
                sumTneg = cellfun(@(x) sum(T(x)), CCneg.PixelIdxList);
                
                T_pos(s) = max([sumTpos 0]);
                T_neg(s) = min([sumTneg 0]);
            end
            toc
            
            
            %% now find significant clusters in the real data
            T_max_thr=prctile(T_pos,100*(1-clusterthreshold));   % sum of
            T_min_thr=prctile(T_neg,100*clusterthreshold);    % sum of
            
            T=squeeze(sum(shf<repmat(mean((task_R-rest_R),1),[size(shf,1),1,1]),1));
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
            
            % Plotting the Task and Rest data
            linePropsT ={'color',[1 0 0]};
            linePropsR ={'color',[0 0 1]};
            
            smoothed_meanT = smoothit(mean(task_R,1),half_win,gaussian_kernel);
            smoothed_stdT = smoothit(std(task_R,1),half_win,gaussian_kernel);
            
            smoothed_meanR = smoothit(mean(rest_R,1),half_win,gaussian_kernel);
            smoothed_stdR = smoothit(std(rest_R,1),half_win,gaussian_kernel);
            
            shadedErrorBar(tfr_time,smoothed_meanT,smoothed_stdT,linePropsT,1); hold on;
            shadedErrorBar(tfr_time,smoothed_meanR,smoothed_stdR,linePropsR,1);
            
            if strcmp(monkey,'Magnus')
                y1 = [-3 3];
            else
                y1 = [-15 10];
                %             yl = ylim;  % get current y-limits of the plot
            end
            
            for i = 1:CCpos.NumObjects
                idx = CCpos.PixelIdxList{i};
                t_start = tfr_time(idx(1));
                t_end = tfr_time(idx(end));
                
                % Draw a semi-transparent rectangle over the significant cluster
                patch([t_start t_end t_end t_start], ...
                    [yl(1) yl(1) yl(2) yl(2)], ...
                    [0.1 0 0.1], ...       % [0.8 0.8 0.8], ...       % gray color
                    'FaceAlpha', 0.3, ...    % transparency
                    'EdgeColor', 'none');
            end
            
            for i = 1:CCneg.NumObjects
                idx = CCneg.PixelIdxList{i};
                t_start = tfr_time(idx(1));
                t_end = tfr_time(idx(end));
                
                % Draw a semi-transparent rectangle over the significant cluster
                patch([t_start t_end t_end t_start], ...
                    [yl(1) yl(1) yl(2) yl(2)], ...
                    [0 0.1 0.1], ...       % [0.8 0.8 0.8], ...       % gray color
                    'FaceAlpha', 0.3, ...    % transparency
                    'EdgeColor', 'none');
            end
            
            set(gca,'ylim',y1);
            xlabel('Time(s)');
            line([0 0], yl, 'color', 'k');
            title([monkey,'-',targets{t},'-',Fin{fin},' - nSite: ', num2str(size(data_task,2)),'-',cond{2},'-',cond{1},'- Significant difference - numPerm =',...
                num2str(numPermutation)],'fontsize', 12,'Interpreter', 'none');
            if filtered
                ylabel('LFP evoked Potential (4Hz filtered)');
                results_file = fullfile([path_to_save,filesep,monkey,'-',targets{t},'-',Fin{fin},'-',cond{2},' vs. ',cond{1},'_Significant difference_4Hzfiltered']);
            else
                ylabel('LFP evoked Potential');
                results_file = fullfile([path_to_save,filesep,monkey,'-',targets{t},'-',Fin{fin},'-',cond{2},' vs. ',cond{1},'_Significant difference']);
            end
            
            export_fig(h(e),[results_file,'.pdf']);
        end
    end
end
close all
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
end

function out=smoothit(in,half_win,gaussian_kernel)
concat_input = cat(2,(in(:,half_win:-1:1)),in);
concat_input = cat(2,concat_input, (in(:,end:-1:end-half_win+1)));
for k=1:size(concat_input,1)
    jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');
end
out = jnk(:,half_win+1:end-half_win);
end
