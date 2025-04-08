function ecg_bna_plot_PEPH(session_info,cfg)

matlab_year = version('-release');
matlab_year = str2double(matlab_year(1:end-1));

dataFolder = [cfg.SPK_root_results_fldr filesep 'cardioballistic' filesep];
condition_colors={cfg.condition.color};
fileList = dir([dataFolder session_info.Monkey(1:3) '_' session_info.Date '*spikes_ECGphase.mat']);
N_conditions = length({cfg.condition.name});

phase_bins = cfg.phase.phase_bin_centers;
phase_bin_edges = [0 phase_bins+(phase_bins(2)-phase_bins(1))/2];
Nbins      = numel(phase_bins);
for u = 1:length(fileList)
    
    load([dataFolder filesep fileList(u).name], 'data')
    sgtitleText = [data.unitId '_' data.target '; ',...
        'FR, Hz: ' num2str(data.FR) '; SNR: ' num2str(data.quantSNR) ...
        '; Fano Factor: ' num2str(data.stability_rating) '; % of ISIs < 3 ms: '  num2str(100 * data.Single_rating) '%'];
    
    %'; ch ' num2str(data.channel) ' AMP_MI = ' num2str(data.(cfg.condition(1).name).AMP_MI(1)) '; ' cfg.condition(1).name ' p = ' num2str(data.(cfg.condition(1).name).AMP_MI(2))
    
    %% R-peak triggered spike response
    figure;
    set(gcf, 'Position', [50    40   450*N_conditions   300])
    
    to_export=0;
    for c=1:N_conditions
        L=cfg.condition(c).name;
        subplot(1,N_conditions,c)
        if data.(L).NrTrials>0
            to_export=1;
            plot_phase_response(data.(L),cfg,c)
        end
    end
    
    if matlab_year < 2016
        mtit(sgtitleText,'interpreter','none'); 
    else
        sgtitle(sgtitleText,'interpreter','none')
    end
    if to_export
    filename= [data.unitId, '__' data.target '_PEPH' ];
    export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    end
    close(gcf);
    
    
    %% Overall picture of criteria
    f1 = figure;
    title('')
    %set(f1, 'Position', [274 148 1452 797])
    
        
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        ttl = ['\color[rgb]{' num2str(cfg.condition(c).color) '}' L];
        
        if data.(L).NrTrials>0
            thr=data.thresholds_microV;
            
            
            % Amp <-> FR correlation things
            % to get a useful polt, we rescale amps
            ampmean=mean(data.(L).AMP_byBin);
            SDsubmean=mean(data.(L).SDsubtractedSDP);
            meandiff=mean(data.(L).AMP_byBin)-SDsubmean;
            scaling= (max(data.(L).SDsubtractedSDP)-min(data.(L).SDsubtractedSDP))/...
                     (max(data.(L).AMP_byBin)-min(data.(L).AMP_byBin(data.(L).AMP_byBin~=0)));
            amp_rescaled    =(data.(L).AMP_byBin    -ampmean)*scaling + SDsubmean;
            ampmax_rescaled =(data.(L).AMP_byBin_max-ampmean)*scaling + SDsubmean;
            ampmin_rescaled =(data.(L).AMP_byBin_min-ampmean)*scaling + SDsubmean;
            thr_rescaled1   =(thr(1)-ampmean)*scaling  + SDsubmean;
            thr_rescaled2   =(thr(2)-ampmean)*scaling  + SDsubmean;
            subplot(2,2,c)
            hold on
            line(phase_bins, data.(L).SDsubtractedSDP, 'Color', condition_colors{c}, 'linewidth',3)
            line(phase_bins, amp_rescaled, 'Color', 'k', 'linewidth',3)
            line(phase_bins, ampmax_rescaled, 'Color', 'k', 'linewidth',1)
            line(phase_bins, ampmin_rescaled, 'Color', 'k', 'linewidth',1)
            line(phase_bins, repmat(thr_rescaled1,size(phase_bins)), 'Color', 'k', 'linestyle',':')
            line(phase_bins, repmat(thr_rescaled2,size(phase_bins)), 'Color', 'k', 'linestyle',':')
            xlabel('Spike Time, radians')
            xlim([0,2*pi]);
            ylabel('FR/Voltage [uV]')
            legend('Normalized FR',['Amp' '-' num2str(meandiff) '*' num2str(scaling)],'threshold','Location','Best')
            title([ttl ', CC: ' num2str(data.(L).FR_amp_cc) ', p: ' num2str(data.(L).FR_amp_p)]);
            
            % threshold things
            subplot(2,2,c+2)
            hold on
            plot(1:250,data.(L).AMP_voltageBinned, 'Color', condition_colors{c});
            lims=ylim(gca);
            plot([thr(1) thr(1)],lims, 'Color', 'k');
            plot([thr(2) thr(2)],lims, 'Color', 'k');
            xlabel('Amplitude [uV]')
            ylabel('Count (N spikes)')
            title([ttl 'Distance : ' num2str(data.(L).distance2thr) 'all high' num2str(data.(L).AMP_all_high)]);
        else
            subplot(2,2,c)
            subplot(2,2,c+2)
        end
    end
    figure(f1)
    if matlab_year < 2016
        mtit(sgtitleText,'xoff',0,'yoff',0.05,'interpreter','none'); 
    else
        sgtitle(sgtitleText,'interpreter','none')
    end
    
    if to_export
    filename= [data.unitId, '_' data.target '_Criteria' ];
    export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    end
    close(gcf);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%     
%     for c=1:numel(cfg.condition)
%         L=cfg.condition(c).name;
%         ttl{c} = ['\color[rgb]{' num2str(cfg.condition(c).color) '}' L];
%         
%         if data.(L).NrTrials>0
%             %ttladd{c}=['_ ' num2str(size(data.(L).waveforms_upsampled_microvolts(1:25:end,:),1)) ' out of ' num2str(size(data.(L).waveforms_upsampled_microvolts,1))];
%         
%             
%             
%             %ax_sp = subplot(2,3,1);
%             subplot(2,3,1);
%             box on
%             hold on
% %             line(cfg.phase.wf_times_interp_ms, data.(L).waveforms_upsampled_microvolts(1:25:end, :), 'Color', [condition_colors{c} 0.1])
% %             line(repmat([cfg.phase.wf_times_interp_ms(1) cfg.phase.wf_times_interp_ms(end)],4,1)', repmat(data.thresholds_microV,1,2)', 'Color', 'r')
%             xlim([cfg.phase.wf_times_interp_ms(1) cfg.phase.wf_times_interp_ms(end)])
%             xlabel('Time, ms')
%             ylabel('Voltage, μV')
%             
%             subplot(2,3,2)
%             box on
%             hold on
%             
%             h = data.(L).SD;
%             p=polar([phase_bins phase_bins(1)],[h h(1)]);%, 'Color', condition_colors{c}(1:3))
%             set(p, 'Color', condition_colors{c}(1:3));
%             
%             
%             subplot(2,3,4)
%             box on
%             hold on
%             %if sum(isnan(data.(L).waveforms_byBin_microvolts))==0 %%???
%             line(cfg.phase.wf_times_interp_ms, data.(L).waveforms_byBin_microvolts', 'Color', [condition_colors{c} 0.1])
%             %end
%             xlabel('Spike Time, ms')
%             ylabel('Voltage, μV')
%             %set(gca, 'XLim', get(ax_sp,'XLim'), 'YLim', get(ax_sp,'YLim'))
%             
%             
% %             %%% 5/6
% %             subplot(2,3,c+4)
% %             title({ttl{c}, 'WF Voltage Change over Heart Cycle'});
% %             imagesc(bsxfun(@minus,data.(L).waveforms_byBin_microvolts,mean(data.(L).waveforms_byBin_microvolts,1)))
% %             cbar = colorbar;
% %             cbar.Title.String = '\Delta Spike Amplitude, μV';
% %             caxis([-3 3])
% %             xlabel('Spike Time')
% %             ylabel('Heart-Cycle Phase')
% %             set(gca, 'XTickLabel', [], 'YTickLabel', [])
%         end
%     end
%     
%     subplot(2,3,1);
%     %title([{'Example Waveforms: '}, ttl, ttladd]);
%     
%     subplot(2,3,2)
%     title([{'ECG-triggered polar PSTH (spike counts / mean spike counts per bin)'}, ttl])
%     
%     subplot(2,3,4)
%     title(['Average WFs by bin:' , ttl])
%     
%     filename= ['Spikes_PolarPSTH__' data.unitId, '_' data.target];
%     export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
%     close(gcf);
    
    
    
    
    
%     %% Distributions of features
%     f2 = figure;
%     set(f2, 'Position', [784   148   942   797])
%     if matlab_year >= 2020
%         sgtitle(sgtitleText, 'interpreter', 'none')
%     end
%     
%     subplot(2,2,1)
%     box on
%     hold on
%     for c=1:numel(cfg.condition)
%         L=cfg.condition(c).name;
%         stairs(data.(L).AMP_voltageBins, data.(L).AMP_voltageBinned, 'Color', condition_colors{c})
%     end
%     hold off
%     line([data.thresholds_microV(1) data.thresholds_microV(1)],ylim, 'Color', 'k')
%     line([data.thresholds_microV(2) data.thresholds_microV(2)],ylim, 'Color', 'k')
%     xlabel('AMP, μV');
%     ylabel('Spike Counts')
%     xlim([0 350])
%     hold off
%     title({'AMP: Same X-axis for All Units', [num2str(data.(L).distance2thr) ' \muV to the closest threshold']})
%     legend({cfg.condition.name}, 'Location', 'Best')
%     
%     subplot(2,2,2)
%     box on
%     hold on
%     for c=1:numel(cfg.condition)
%         L=cfg.condition(c).name;
%         [h,bins]=hist(abs(data.(L).HW_ms), [0:0.005:0.5]);
%         stairs(bins,h,'Color', condition_colors{c})
%     end
%     hold off
%     xlim([0 0.5])
%     xlabel('HW, ms');
%     ylabel('Spike Counts')
%     title('HW: Same X-axis for All Units')
%     
%     subplot(2,2,3)
%     box on
%     hold on
%     for c=1:numel(cfg.condition)
%         L=cfg.condition(c).name;
%         [h,bins]=hist(abs(data.(L).TPW_ms), [0.2:0.005:1]);
%         stairs(bins,h,'Color', condition_colors{c})
%     end
%     hold off
%     xlabel('TPW, ms');
%     ylabel('Spike Counts')
%     title('TPW: X-axis Adjusts for Each Unit')
%     
%     subplot(2,2,4)
%     box on
%     hold on
%     for c=1:numel(cfg.condition)
%         L=cfg.condition(c).name;
%         [h,bins]=hist(abs(data.(L).REP_ms),[0:0.005:0.5]);
%         stairs(bins,h,'Color', condition_colors{c})
%     end
%     hold off
%     xlabel('REP, ms');
%     ylabel('Spike Counts')
%     title('REP: X-axis Adjusts for Each Unit')
%     
%     filename= ['Distributions_AMP_HW_TPW_REP__' data.unitId, '_' data.target];
%     export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
%     close(gcf);
%     
%     %% Real and reshuffled data
%     f3 = figure;
%     set(f3, 'Position', [1 41 1920 963])
%     if matlab_year >= 2020
%         sgtitle(sgtitleText, 'interpreter', 'none')
%     end
%     feature_list = {'AMP', 'HW', 'TPW', 'REP'};
%     feature_bin_list = {'AMP_microV_byBin', 'HW_ms_byBin', 'TPW_ms_byBin', 'REP_ms_byBin'};
%     smoothed_feature_list = {'AMP_microV_byBin_smoothed', 'HW_ms_byBin_smoothed', 'TPW_ms_byBin_smoothed', 'REP_ms_byBin_smoothed'};
%     subplot_numbers = {[1 2], [4 5], [11 12], [14 15]};
%     
%     for featureNum = 1:4
%         for c=1:numel(cfg.condition)
%             L=cfg.condition(c).name;
%             subplot(3,5,subplot_numbers{featureNum}(c))
%             yyaxis left
%             bar(phase_bins, data.(L).spike_phases_histogram, 'FaceColor', condition_colors{c}(1:3)/2)
%             hold on
%             plot(phase_bins, data.(L).spike_phases_histogram_smoothed, '-', 'Color', condition_colors{c}(1:3), 'LineWidth',2)
%             ylabel('Spike Counts')
%             currYLims = get(gca, 'YLim');
%             currYLims(2) = 2 * currYLims(2);
%             set(gca, 'YLim', currYLims)
%             yyaxis right
%             filledArea = fill([phase_bins fliplr(phase_bins) phase_bins(1)], ...
%                 ...
%                 [data.(L).([feature_list{featureNum} '_upperPrctile_97_5']) ...
%                 fliplr(data.(L).([feature_list{featureNum} '_lowerPrctile_2_5'])) ...
%                 data.(L).([feature_list{featureNum} '_upperPrctile_97_5'])(1)], ...
%                 [0 0 0], 'FaceALpha', 0.15, 'EdgeColor', 'none');
%             hold on;
%             p1 = plot(phase_bins, data.(L).(feature_bin_list{featureNum}),'-k','LineWidth',2);
%             p2 = plot(phase_bins, data.(L).(smoothed_feature_list{featureNum}), '-', 'Color', condition_colors{c}(1:3), 'LineWidth',2);
%             currYLims = get(gca, 'YLim');
%             currYLims(1) = currYLims(1) - diff(currYLims);
%             set(gca, 'YLim', currYLims)
%             title({[feature_list{featureNum} ': MI = ' num2str(data.(L).([feature_list{featureNum} '_MI'])(1)) '; p = ' num2str(data.(L).([feature_list{featureNum} '_MI'])(2))], ...
%                 ['Max consec. bins: ' num2str(data.(L).([feature_list{featureNum} '_max_consec_bins']))], ...
%                 ['Corr. w/PSTH: cc = ' num2str(data.(L).cc_PSTH_feature(featureNum)) '; p = ' num2str(data.(L).pp_PSTH_feature(featureNum))]})
%             xlim([0 2*pi])
%             xlabel('Heart-cycle Phase (0-2pi)')
%             ylabel([feature_list{featureNum} ', microvolts'])
%             if c == 1
%                 %                 legend([filledArea p1 p2], {'95% Confidence Interval', 'Average by Bin', 'Rlowess-Smoothed'}, 'Location', 'southoutside')
%             end
%         end
%     end
%     
%     filename= ['PSTH_overlaid_Feature_Dynamics__' data.unitId, '_' data.target];
%     export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
%     close(gcf);
%     
%     % cosine fitted plots
%     f4 = figure;
%     set(f4, 'Position', [1 41 1920 963])
%     if matlab_year >= 2020
%         sgtitle(sgtitleText, 'interpreter', 'none')
%     end
%     
%     for c=1:numel(cfg.condition)
%         col=condition_colors{c};
%         L=cfg.condition(c).name;
%         subplot(3,5,c)
%         
%         yyaxis left
%         set(gca,'YColor',col/2);
%         plot(data.(L).spike_phases_radians, data.(L).AMP_microV, '.', 'Color', col/2)
%         % adjust yaxis
%         ylims = get(gca, 'YLim');
%         ylims(2) = ylims(2) + (ylims(2) - ylims(1));
%         set(gca, 'YLim', ylims)
%         ylabel('Abs. Spike AMP, \muV')
%         
%         yyaxis right
%         set(gca,'YColor','k');
%         plot(phase_bins, 100*data.(L).AMP_microV_byBin'/nanmean(data.(L).AMP_microV_byBin_smoothed),':k','LineWidth',1);
%         hold on
%         plot(phase_bins, ...
%             100*(1+data.(L).AMP_MI(1)*cos(phase_bins-data.(L).AMP_MI(3))), 'Color', col,'LineWidth',2);
%         plot(phase_bins, ...
%             100*data.(L).AMP_microV_byBin_smoothed / nanmean(data.(L).AMP_microV_byBin_smoothed), ...
%             '--', 'Color', col,'LineWidth',2)
%         title(['Motion Index, %: ' num2str(100 * data.(L).AMP_MI(1)) ...
%             ', p = ' num2str(data.(L).AMP_MI(2))])
%         xlabel('Heart-Cycle Phase (0-2pi)')
%         ylabel({'% AMP Signal Change,', 'signal / average)'})
%         % adjust yaxis
%         if ~isnan(data.(L).spike_phases_radians)
%             ylims = get(gca, 'YLim');
%             ylims(1) = ylims(1) - (ylims(2) - ylims(1));
%             set(gca, 'YLim', ylims)
%         end
%         
%         subplot(3,5,3+c)
%         yyaxis left
%         set(gca,'YColor',col/2);
%         plot(data.(L).spike_phases_radians, data.(L).HW_ms, '.', 'Color', col/2)
%         % adjust yaxis
%         ylims = get(gca, 'YLim');
%         ylims(2) = ylims(2) + (ylims(2) - ylims(1));
%         set(gca, 'YLim', ylims)
%         ylabel('Spike HW, \muV')
%         
%         yyaxis right
%         set(gca,'YColor','k');
%         plot(phase_bins, 100*data.(L).HW_ms_byBin'/nanmean(data.(L).HW_ms_byBin_smoothed),':k','LineWidth',1);
%         hold on
%         plot(phase_bins, ...
%             100*(1+data.(L).HW_MI(1)*cos(phase_bins-data.(L).HW_MI(3))), ...
%             'Color', col,'LineWidth',2);
%         plot(phase_bins, ...
%             100*data.(L).HW_ms_byBin_smoothed / nanmean(data.(L).HW_ms_byBin_smoothed), ...
%             '--', 'Color', col,'LineWidth',2)
%         title(['Motion Index, %: ' num2str(100 * data.(L).HW_MI(1)) ...
%             ', p = ' num2str(data.(L).HW_MI(2))])
%         xlabel('Heart-Cycle Phase (0-2pi)')
%         ylabel({'% HW Signal Change,', 'signal / average'})
%         % adjust yaxis
%         if ~isnan(data.(L).spike_phases_radians)
%             ylims = get(gca, 'YLim');
%             ylims(1) = ylims(1) - (ylims(2) - ylims(1));
%             set(gca, 'YLim', ylims)
%         end
%         
%         subplot(3,5,10+c)
%         yyaxis left
%         set(gca,'YColor',col/2);
%         plot(data.(L).spike_phases_radians, data.(L).HW_ms, '.', 'Color', col/2)
%         % adjust yaxis
%         ylims = get(gca, 'YLim');
%         ylims(2) = ylims(2) + (ylims(2) - ylims(1));
%         set(gca, 'YLim', ylims)
%         ylabel('Spike TPW, \muV')
%         
%         yyaxis right
%         set(gca,'YColor','k');
%         plot(phase_bins, 100*data.(L).TPW_ms_byBin'/nanmean(data.(L).TPW_ms_byBin_smoothed),':k','LineWidth',1);
%         hold on
%         plot(phase_bins, ...
%             100*(1+data.(L).TPW_MI(1)*cos(phase_bins-data.(L).TPW_MI(3))), ...
%             'Color', col,'LineWidth',2);
%         plot(phase_bins, ...
%             100*data.(L).TPW_ms_byBin_smoothed / nanmean(data.(L).TPW_ms_byBin_smoothed), ...
%             '--', 'Color', col,'LineWidth',2)
%         title(['Motion Index, %: ' num2str(100 * data.(L).TPW_MI(1)) ...
%             ', p = ' num2str(data.(L).TPW_MI(2))])
%         xlabel('Heart-Cycle Phase (0-2pi)')
%         ylabel({'% TPW Signal Change,', 'signal / average'})
%         % adjust yaxis
%         if ~isnan(data.(L).spike_phases_radians)
%             ylims = get(gca, 'YLim');
%             ylims(1) = ylims(1) - (ylims(2) - ylims(1));
%             set(gca, 'YLim', ylims)
%         end
%         
%         subplot(3,5,13+c)
%         yyaxis left
%         set(gca,'YColor',col/2);
%         plot(data.(L).spike_phases_radians, data.(L).REP_ms, '.', 'Color', col/2)
%         % adjust yaxis
%         ylims = get(gca, 'YLim');
%         ylims(2) = ylims(2) + (ylims(2) - ylims(1));
%         set(gca, 'YLim', ylims)
%         ylabel('Spike REP, \muV')
%         
%         yyaxis right
%         set(gca,'YColor','k');
%         plot(phase_bins, 100*data.(L).REP_ms_byBin'/nanmean(data.(L).REP_ms_byBin_smoothed),':k','LineWidth',1);
%         hold on
%         plot(phase_bins, ...
%             100*(1+data.(L).REP_MI(1)*cos(phase_bins-data.(L).REP_MI(3))), ...
%             'Color', col,'LineWidth',2);
%         plot(phase_bins, ...
%             100*data.(L).REP_ms_byBin_smoothed / nanmean(data.(L).REP_ms_byBin_smoothed), ...
%             '--', 'Color', col,'LineWidth',2)
%         title(['Motion Index, %: ' num2str(100 * data.(L).REP_MI(1)) ...
%             ', p = ' num2str(data.(L).REP_MI(2))])
%         xlabel('Heart-Cycle Phase (0-2pi)')
%         ylabel({'% REP Signal Change,', 'signal / average'})
%         % adjust yaxis
%         if ~isnan(data.(L).spike_phases_radians)
%             ylims = get(gca, 'YLim');
%             ylims(1) = ylims(1) - (ylims(2) - ylims(1));
%             set(gca, 'YLim', ylims)
%         end
%     end
%     
%     filename= ['Cosine_Fitted__' data.unitId, '_' data.target];
%     export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
%     close(gcf);
    

%% THIS HERE IS REMOVED


    %% histogram with AIC and BIC
    %     figure,
    %     set(gcf, 'Position', [381 607 821 389])
    %     if matlab_year >= 2020
    %         sgtitle(sgtitleText, 'interpreter', 'none')
    %     end
    %
    %     for c=1:numel(cfg.condition)
    %         L=cfg.condition(c).name;
    %
    %         [maxAIC, maxAICid] = max([data.(L).linear.aic data.(L).cosine.aic data.(L).vonMisesPos.aic data.(L).vonMisesNeg.aic]);
    %         [maxBIC, maxBICid] = max([data.(L).linear.bic data.(L).cosine.bic data.(L).vonMisesPos.bic data.(L).vonMisesNeg.bic]);
    %
    %         subplot(1,2,c)
    %         bar([data.(L).linear.aic data.(L).cosine.aic data.(L).vonMisesPos.aic data.(L).vonMisesNeg.aic; ...
    %             data.(L).linear.bic data.(L).cosine.bic data.(L).vonMisesPos.bic data.(L).vonMisesNeg.bic]')
    %         hold on
    %         if ~isnan(maxAIC)
    %             plot(maxAICid, 0, 'or')
    %         end
    %         if ~isnan(maxBIC)
    %             plot(maxBICid, 0, 'xb')
    %         end
    %         set(gca, 'XTickLabel', {'linear', 'cosine', 'pos. von Mises', 'neg. von Mises'})
    %         xlabel('Type of Fit')
    %         ylabel('AIC or BIC value')
    %     end
    %
    %     filename= ['AIC_BIC_' data.unitId, '_' data.target];
    %     export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    %     close(gcf);
    
end


end

function plot_phase_response(data,cfg,c)

phase_bins=cfg.phase.phase_bin_centers;
hold on
box on
xlim([0 2*pi])
col=cfg.condition(c).color;
lineProps={'color',col,'linewidth',1};
shadedErrorBar(phase_bins, data.SD, data.SD_SEM,lineProps,1);
lineProps={'color',col,'linewidth',1,'linestyle',':'};
shadedErrorBar(phase_bins, data.SDP, [data.SDPCu-data.SDP; data.SDP-data.SDPCL], lineProps,1);
%             line([0 0],ylim,'color','k');
if data.sig_sign==-1
    ypos=min(data.SD)*-1;
elseif data.sig_sign==1
    ypos=max(data.SD);
else
    ypos=NaN;
end
to_plot=data.sig;
to_plot(to_plot~=data.sig_sign)=NaN;
plot(phase_bins,to_plot*ypos,'color',col,'linewidth',5);

%             y_lims=get(gca,'ylim');
%             text(BINS(10),y_lims(2)-diff(y_lims)*c*1/20, [L ': trials = ' ,num2str(data.(L).NrTrials), ' events = ' ,num2str(data.(L).NrEvents) ],'Color',condition_colors{c}, 'Interpreter', 'none');

%             set(gca, 'XTick', [-200 -100 0 100 200], 'XTickLabel', [-200 -100 0 100 200])
ylabel('Firing Rate, Hz');
xlabel('Phase of the Heart Cycle [0 2\pi]');

end


    
    %% fits from Mosher, and cosine and von Mises by Luba
    %     f0 = figure;
    %     set(f0, 'Position', [2 38 500*numel(cfg.condition) 958])
    %     if matlab_year >= 2020
    %         sgtitle(sgtitleText, 'interpreter', 'none')
    %     end
    %     for c=1:numel(cfg.condition)
    %         L=cfg.condition(c).name;
    %
    %         subplot(4,length(cfg.condition),c)
    %         % cosine fit of the phase histogram (Mosher's)
    %         bar(phase_bins, data.(L).spike_phases_histogram, 'FaceColor', cfg.condition(c).color/2)
    %         hold on
    %         line(phase_bins, data.(L).spike_phases_histogram_smoothed, 'Color', cfg.condition(c).color, 'LineWidth', 2)
    %         plot(phase_bins, ...
    %             mean(data.(L).spike_phases_histogram_smoothed)* (1+data.(L).histogram_MI*cos(phase_bins-data.(L).histogram_phase)), ...
    %             '--', 'Color', cfg.condition(c).color, 'LineWidth', 2);
    %         if ~isnan(data.(L).histogram_phase)
    %             xline(data.(L).histogram_phase, ':', 'Color', cfg.condition(c).color)
    %         end
    %         xlim([0 2*pi])
    %         ymin = min(data.(L).spike_phases_histogram);
    %         if ~isnan(ymin)
    %             ylims = get(gca, 'YLim');
    %             ylims(1) = ymin;% update the lower ylim
    %             ylim(ylims)
    %         end
    %         title({[cfg.condition(c).name ': MI = ' num2str(data.(L).histogram_MI) '; p = ' num2str(data.(L).histogram_p)], ...
    %             ['Mean Phase = ' num2str(data.(L).histogram_phase) '; R-squared = ' num2str(data.(L).rsquared)]})
    %         xlabel('Phase, radians')
    %         ylabel('Spike Counts')
    %         legend({'Phase PSTH', 'Smoothed', 'Cosine Fit (Mosher)', 'Circ.Mean Phase'}, 'Location', 'best')
    %
    %         distList = {'cosine', 'vonMisesPos', 'vonMisesNeg'};
    %         for distNum = 1:3
    %             currFit = distList{distNum};
    %             subplot(4,length(cfg.condition),distNum*length(cfg.condition)+c)
    %             % cosine fits by Luba
    %             yyaxis left
    %             imagesc(phase_bins, 1:size(data.(L).spike_phases_histogram2,2), data.(L).spike_phases_histogram2')
    %             colormap(flipud(gray))
    %             ylabel('# Heart Cycle')
    %
    %             yyaxis right
    %             line(phase_bins, data.(L).(currFit).average, 'Color', cfg.condition(c).color, 'LineWidth', 2)
    %             hold on
    %             plot(phase_bins, data.(L).(currFit).yfit, ...
    %                 '--', 'Color', cfg.condition(c).color, 'LineWidth', 2)
    %
    %             if strcmp(currFit, 'cosine')
    %                 if ~isnan(data.(L).(currFit).coefs(2))
    %                     xline(data.(L).(currFit).coefs(2), ':', 'Color', cfg.condition(c).color)
    %                 end
    %                 title({['Mean Phase = ' num2str(data.(L).(currFit).coefs(2)) '; R^2 = ' num2str(data.(L).(currFit).rsquared) '; p = ' num2str(data.(L).(currFit).pvalue)], ...
    %                     ['Linear: R^2 = ' num2str(data.(L).linear.rsquared) '; p = ' num2str(data.(L).linear.pvalue(2))]})
    %
    %                 % plot results of the linear fit
    %                 plot(phase_bins, data.(L).linear.yfit, ':k', 'LineWidth', 1.5)
    %             else
    %                 if ~isnan(data.(L).(currFit).coefs(4))
    %                     xline(data.(L).(currFit).coefs(4), ':', 'Color', cfg.condition(c).color)
    %                 end
    %                 title({['Mean Phase = ' num2str(data.(L).(currFit).coefs(4)) '; R^2 = ' num2str(data.(L).(currFit).rsquared) '; p = ' num2str(data.(L).(currFit).pvalue)], ...
    %                     ['a1 = ' num2str(data.(L).(currFit).coefs(1)) '; d1 = ' num2str(data.(L).(currFit).coefs(2)) '; \kappa = ' num2str(data.(L).(currFit).coefs(3))]})
    %             end
    %
    %             xlim([0 2*pi])
    %             ylabel('Spikes per bin across heart-cycles')
    %             legend({'Phase PSTH', 'Circular Fit', 'Circ.Mean Phase'}, 'Location', 'best')
    %             box on
    %         end
    %
    %     end
    %
    %     filename = ['PhasePSTH__' data.unitId, '_' data.target];
    %     export_fig(gcf, [dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    %     close(gcf);