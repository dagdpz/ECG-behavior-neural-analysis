function ecg_bna_plot_session_ECG_related_spikePhase(session_info,cfg)

matlab_year = version('-release');
matlab_year = str2double(matlab_year(1:end-1));

dataFolder = [cfg.SPK_root_results_fldr filesep 'cardioballistic' filesep];

condition_colors={cfg.condition.color};

fileList = dir([dataFolder session_info.Monkey(1:3) '_' session_info.Date '*spikes_ECGphase.mat']);

for untNum = 1:length(fileList)
    
    load([dataFolder filesep fileList(untNum).name], 'data')
    
    % with the current Fs and number of samples per spike I have 1.4 ms per
    % spike
    
    %% title with all the spike analysis parameters
    sgtitleText = {[data.unitId '_' data.target '; ch ' num2str(data.channel) ';'  ' unit ' data.unit], ... %  
        ['FR, Hz: ' num2str(data.FR) '; SNR: ' num2str(data.quantSNR) '; Fano Factor: ' num2str(data.stability_rating) '; % of ISIs < 3 ms: '  num2str(100 * data.Single_rating) '%'], ...
        [cfg.condition(1).name ' AMP_MI = ' num2str(data.(cfg.condition(1).name).AMP_MI(1)) '; ' cfg.condition(1).name ' p = ' num2str(data.(cfg.condition(1).name).AMP_MI(2))]};
    
    %% fits from Mosher, and cosine and von Mises by Luba
    f0 = figure;
    set(f0, 'Position', [2 38 500*numel(cfg.condition) 958])
    if matlab_year >= 2020
        sgtitle(sgtitleText, 'interpreter', 'none')
    end
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        
        subplot(4,length(cfg.condition),c)
        % cosine fit of the phase histogram (Mosher's)
        bar(cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram, 'FaceColor', cfg.condition(c).color/2)
        hold on
        line(cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram_smoothed, 'Color', cfg.condition(c).color, 'LineWidth', 2)
        plot(cfg.spk.phase_bin_centers, ...
            mean(data.(L).spike_phases_histogram_smoothed)* (1+data.(L).histogram_MI*cos(cfg.spk.phase_bin_centers-data.(L).histogram_phase)), ...
            '--', 'Color', cfg.condition(c).color, 'LineWidth', 2);
        if ~isnan(data.(L).histogram_phase)
            xline(data.(L).histogram_phase, ':', 'Color', cfg.condition(c).color)
        end
        xlim([0 2*pi])
        ymin = min(data.(L).spike_phases_histogram);
        if ~isnan(ymin)
            ylims = get(gca, 'YLim');
            ylims(1) = ymin;% update the lower ylim
            ylim(ylims)
        end
        title({[cfg.condition(c).name ': MI = ' num2str(data.(L).histogram_MI) '; p = ' num2str(data.(L).histogram_p)], ...
            ['Mean Phase = ' num2str(data.(L).histogram_phase) '; R-squared = ' num2str(data.(L).rsquared)]})
        xlabel('Phase, radians')
        ylabel('Spike Counts')
        legend({'Phase PSTH', 'Smoothed', 'Cosine Fit (Mosher)', 'Circ.Mean Phase'}, 'Location', 'best')
        
        distList = {'cosine', 'vonMisesPos', 'vonMisesNeg'};
        for distNum = 1:3
            currFit = distList{distNum};
            subplot(4,length(cfg.condition),distNum*length(cfg.condition)+c)
            % cosine fits by Luba
            yyaxis left
            imagesc(cfg.spk.phase_bin_centers, 1:size(data.(L).spike_phases_histogram2,2), data.(L).spike_phases_histogram2')
            colormap(flipud(gray))
            ylabel('# Heart Cycle')
%         plot(cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram2, '.', 'Color', [cfg.condition(c).color 0.1])
        
            yyaxis right
            line(cfg.spk.phase_bin_centers, data.(L).(currFit).average, 'Color', cfg.condition(c).color, 'LineWidth', 2)
            hold on
            plot(cfg.spk.phase_bin_centers, data.(L).(currFit).yfit, ...
                '--', 'Color', cfg.condition(c).color, 'LineWidth', 2)
%         line(cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram_smoothed, 'Color', cfg.condition(c).color, 'LineWidth', 2)
            if strcmp(currFit, 'cosine')
                if ~isnan(data.(L).(currFit).coefs(2))
                    xline(data.(L).(currFit).coefs(2), ':', 'Color', cfg.condition(c).color)
                end
                title({['Mean Phase = ' num2str(data.(L).(currFit).coefs(2)) '; R^2 = ' num2str(data.(L).(currFit).rsquared) '; p = ' num2str(data.(L).(currFit).pvalue)], ...
                    ['Linear: R^2 = ' num2str(data.(L).linear.rsquared) '; p = ' num2str(data.(L).linear.pvalue(2))]})
                
                % plot results of the linear fit
                plot(cfg.spk.phase_bin_centers, data.(L).linear.yfit, ':k', 'LineWidth', 1.5)
            else
                if ~isnan(data.(L).(currFit).coefs(4))
                    xline(data.(L).(currFit).coefs(4), ':', 'Color', cfg.condition(c).color)
                end
                title({['Mean Phase = ' num2str(data.(L).(currFit).coefs(4)) '; R^2 = ' num2str(data.(L).(currFit).rsquared) '; p = ' num2str(data.(L).(currFit).pvalue)], ...
                    ['a1 = ' num2str(data.(L).(currFit).coefs(1)) '; d1 = ' num2str(data.(L).(currFit).coefs(2)) '; \kappa = ' num2str(data.(L).(currFit).coefs(3))]})
            end
            
            xlim([0 2*pi])
            ylabel('Spikes per bin across heart-cycles')
            legend({'Phase PSTH', 'Circular Fit', 'Circ.Mean Phase'}, 'Location', 'best')
            box on
        end
        
    end
    
    filename = ['PhasePSTH__' data.unitId, '_' data.target];
    export_fig(gcf, [dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    %% Overall picture of waveforms and PSTHs
    f1 = figure;
    set(f1, 'Position', [274 148 1452 797])
    if matlab_year >= 2020
        sgtitle(sgtitleText, 'interpreter', 'none')
    end
    
    ax_sp = subplot(2,3,1);
    box on
    hold on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        line(cfg.spk.wf_times_interp_ms, data.(L).waveforms_upsampled_microvolts(1:25:end, :), 'Color', [condition_colors{c} 0.1])
        line(repmat([cfg.spk.wf_times_interp_ms(1) cfg.spk.wf_times_interp_ms(end)],4,1)', repmat(data.thresholds_microV,1,2)', 'Color', 'r')
    end
    xlim([cfg.spk.wf_times_interp_ms(1) cfg.spk.wf_times_interp_ms(end)])
    xlabel('Time, ms')
    ylabel('Voltage, μV')
    % create title
    ttl = {};
    for c = 1:length(cfg.condition)
        ttl{c} = ['\color[rgb]{' num2str(cfg.condition(c).color) '}' ...
            cfg.condition(c).name '_ ' num2str(size(data.(cfg.condition(c).name).waveforms_upsampled_microvolts(1:25:end,:),1)) ...
            ' out of ' num2str(size(data.(cfg.condition(c).name).waveforms_upsampled_microvolts,1))];
    end
    title([{'Example Waveforms: '}, ttl])
    
    subplot(2,3,2)
    box on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        h = hist(data.(L).spike_phases_radians, cfg.spk.phase_bin_centers);
        h = h / mean(h);
        p=polar([cfg.spk.phase_bin_centers cfg.spk.phase_bin_centers(1)],[h h(1)]);%, 'Color', condition_colors{c}(1:3))
        set(p, 'Color', condition_colors{c}(1:3));
        hold on
    end
    hold off
    ttl = {};
    for c = 1:length(cfg.condition)
        ttl{c} = ['\color[rgb]{' num2str(cfg.condition(c).color) '}' cfg.condition(c).name];
    end
    title([{'ECG-triggered polar PSTH (spike counts / mean spike counts per bin)'}, ttl])
%     axis tight
    
    subplot(2,3,4)
    title('Average WFs by bin: \color{blue}Rest and \color{red}Task')
    hold on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        if sum(isnan(data.(L).waveforms_byBin_microvolts))==0
            line(cfg.spk.wf_times_interp_ms, data.(L).waveforms_byBin_microvolts', 'Color', [condition_colors{c} 0.1])
        end
    end
    hold off
    box on
    xlabel('Spike Time, ms')
    ylabel('Voltage, μV')
    set(gca, 'XLim', get(ax_sp,'XLim'), 'YLim', get(ax_sp,'YLim'))
    
    subplot(2,3,5)
    imagesc(bsxfun(@minus,data.(cfg.condition(1).name).waveforms_byBin_microvolts,mean(data.(cfg.condition(1).name).waveforms_byBin_microvolts,1)))
    cbar = colorbar;
    cbar.Title.String = '\Delta Spike Amplitude, μV';
    caxis([-3 3])
    title({['{\color[rgb]{' num2str(cfg.condition(1).color) '}[' cfg.condition(1).name ']}'], 'WF Voltage Change over Heart Cycle'})
    xlabel('Spike Time')
    ylabel('Heart-Cycle Phase')
    set(gca, 'XTickLabel', [], 'YTickLabel', [])
    
    subplot(2,3,6)
    imagesc(bsxfun(@minus,data.(cfg.condition(2).name).waveforms_byBin_microvolts,mean(data.(cfg.condition(2).name).waveforms_byBin_microvolts,1)))
    cbar = colorbar;
    cbar.Title.String = '\Delta Spike Amplitude, μV';
    caxis([-3 3])
    title({['{\color[rgb]{' num2str(cfg.condition(2).color) '}[' cfg.condition(2).name ']}'], 'WF Voltage Change over Heart Cycle'})
    xlabel('Spike Time')
    ylabel('Heart-Cycle Phase')
    set(gca, 'XTickLabel', [], 'YTickLabel', [])
    
    filename= ['Spikes_PolarPSTH__' data.unitId, '_' data.target];
    export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    %% Distributions of features
    f2 = figure;
    set(f2, 'Position', [784   148   942   797])
    if matlab_year >= 2020
        sgtitle(sgtitleText, 'interpreter', 'none')
    end
    
    subplot(2,2,1)
    box on
    hold on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        stairs(data.(L).AMP_voltageBins, data.(L).AMP_voltageBinned, 'Color', condition_colors{c})
    end
    hold off
    line([data.thresholds_microV(1) data.thresholds_microV(1)],ylim, 'Color', 'k')
    line([data.thresholds_microV(2) data.thresholds_microV(2)],ylim, 'Color', 'k')
    xlabel('AMP, μV');
    ylabel('Spike Counts')
    xlim([0 350])
    hold off
    title({'AMP: Same X-axis for All Units', [num2str(data.(L).distance2thr) ' \muV to the closest threshold']})
    legend({cfg.condition.name}, 'Location', 'Best')
    
    subplot(2,2,2)
    box on
    hold on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        [h,bins]=hist(abs(data.(L).HW_ms), [0:0.005:0.5]);
        stairs(bins,h,'Color', condition_colors{c})
    end
    hold off
    xlim([0 0.5])
    xlabel('HW, ms');
    ylabel('Spike Counts')
    title('HW: Same X-axis for All Units')
    
    subplot(2,2,3)
    box on
    hold on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        [h,bins]=hist(abs(data.(L).TPW_ms), [0.2:0.005:1]);
        stairs(bins,h,'Color', condition_colors{c})
    end
    hold off
    xlabel('TPW, ms');
    ylabel('Spike Counts')
    title('TPW: X-axis Adjusts for Each Unit')
    
    subplot(2,2,4)
    box on
    hold on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        [h,bins]=hist(abs(data.(L).REP_ms),[0:0.005:0.5]);
        stairs(bins,h,'Color', condition_colors{c})
    end
    hold off
    xlabel('REP, ms');
    ylabel('Spike Counts')
    title('REP: X-axis Adjusts for Each Unit')
    
    filename= ['Distributions_AMP_HW_TPW_REP__' data.unitId, '_' data.target];
    export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    %% Real and reshuffled data
    f3 = figure;
    set(f3, 'Position', [1 41 1920 963])
    if matlab_year >= 2020
        sgtitle(sgtitleText, 'interpreter', 'none')
    end
    feature_list = {'AMP', 'HW', 'TPW', 'REP'};
    feature_bin_list = {'AMP_microV_byBin', 'HW_ms_byBin', 'TPW_ms_byBin', 'REP_ms_byBin'};
    smoothed_feature_list = {'AMP_microV_byBin_smoothed', 'HW_ms_byBin_smoothed', 'TPW_ms_byBin_smoothed', 'REP_ms_byBin_smoothed'};
    subplot_numbers = {[1 2], [4 5], [11 12], [14 15]};
    
    for featureNum = 1:4
        for c=1:numel(cfg.condition)
            L=cfg.condition(c).name;
            subplot(3,5,subplot_numbers{featureNum}(c))
            yyaxis left
            histogram(data.(L).spike_phases_radians, cfg.spk.phase_bin_centers, 'FaceColor', condition_colors{c}(1:3))
            ylabel('Spike Counts')
            currYLims = get(gca, 'YLim');
            currYLims(2) = 2 * currYLims(2);
            set(gca, 'YLim', currYLims)
            yyaxis right
            filledArea = fill([cfg.spk.phase_bin_centers fliplr(cfg.spk.phase_bin_centers) cfg.spk.phase_bin_centers(1)], ...
                ...
                [data.(L).([feature_list{featureNum} '_upperPrctile_97_5']) ...
                fliplr(data.(L).([feature_list{featureNum} '_lowerPrctile_2_5'])) ...
                data.(L).([feature_list{featureNum} '_upperPrctile_97_5'])(1)], ...
                [0 0 0], 'FaceALpha', 0.15, 'EdgeColor', 'none');
            hold on;
            p1 = plot(cfg.spk.phase_bin_centers, data.(L).(feature_bin_list{featureNum}),'-k','LineWidth',2);
            p2 = plot(cfg.spk.phase_bin_centers, data.(L).(smoothed_feature_list{featureNum}), '-', 'Color', condition_colors{c}(1:3), 'LineWidth',2);
            currYLims = get(gca, 'YLim');
            currYLims(1) = currYLims(1) - diff(currYLims);
            set(gca, 'YLim', currYLims)
            title({[feature_list{featureNum} ': MI = ' num2str(data.(L).([feature_list{featureNum} '_MI'])(1)) '; p = ' num2str(data.(L).([feature_list{featureNum} '_MI'])(2))], ...
                ['Max consec. bins: ' num2str(data.(L).([feature_list{featureNum} '_max_consec_bins']))], ...
                ['Corr. w/PSTH: cc = ' num2str(data.(L).cc_PSTH_feature(featureNum)) '; p = ' num2str(data.(L).pp_PSTH_feature(featureNum))]})
            xlim([0 2*pi])
            xlabel('Heart-cycle Phase (0-2pi)')
            ylabel([feature_list{featureNum} ', microvolts'])
            if c == 1
%                 legend([filledArea p1 p2], {'95% Confidence Interval', 'Average by Bin', 'Rlowess-Smoothed'}, 'Location', 'southoutside')
            end
        end
    end
    
    filename= ['PSTH_overlaid_Feature_Dynamics__' data.unitId, '_' data.target];
    export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    % cosine fitted plots
    f4 = figure;
    set(f4, 'Position', [1 41 1920 963])
    if matlab_year >= 2020
        sgtitle(sgtitleText, 'interpreter', 'none')
    end
    
    for c=1:numel(cfg.condition)
        col=condition_colors{c};
        L=cfg.condition(c).name;
        subplot(3,5,c)
        
        yyaxis left
        set(gca,'YColor',col/2);
        plot(data.(L).spike_phases_radians, data.(L).AMP_microV, '.', 'Color', col/2)
        % adjust yaxis
        ylims = get(gca, 'YLim');
        ylims(2) = ylims(2) + (ylims(2) - ylims(1));
        set(gca, 'YLim', ylims)
        ylabel('Abs. Spike AMP, \muV')
        
        yyaxis right
        set(gca,'YColor','k');
        plot(cfg.spk.phase_bin_centers, 100*data.(L).AMP_microV_byBin'/nanmean(data.(L).AMP_microV_byBin_smoothed),':k','LineWidth',1);
        hold on
        plot(cfg.spk.phase_bin_centers, ...
            100*(1+data.(L).AMP_MI(1)*cos(cfg.spk.phase_bin_centers-data.(L).AMP_MI(3))), 'Color', col,'LineWidth',2);
        plot(cfg.spk.phase_bin_centers, ...
            100*data.(L).AMP_microV_byBin_smoothed / nanmean(data.(L).AMP_microV_byBin_smoothed), ...
            '--', 'Color', col,'LineWidth',2)
        title(['Motion Index, %: ' num2str(100 * data.(L).AMP_MI(1)) ...
            ', p = ' num2str(data.(L).AMP_MI(2))])
        xlabel('Heart-Cycle Phase (0-2pi)')
        ylabel({'% AMP Signal Change,', 'signal / average)'})
        % adjust yaxis
        if ~isnan(data.(L).spike_phases_radians)  
            ylims = get(gca, 'YLim');
            ylims(1) = ylims(1) - (ylims(2) - ylims(1));
            set(gca, 'YLim', ylims)
        end
        
        subplot(3,5,3+c)
        yyaxis left
        set(gca,'YColor',col/2);
        plot(data.(L).spike_phases_radians, data.(L).HW_ms, '.', 'Color', col/2)
        % adjust yaxis
        ylims = get(gca, 'YLim');
        ylims(2) = ylims(2) + (ylims(2) - ylims(1));
        set(gca, 'YLim', ylims)
        ylabel('Spike HW, \muV')
        
        yyaxis right
        set(gca,'YColor','k');
        plot(cfg.spk.phase_bin_centers, 100*data.(L).HW_ms_byBin'/nanmean(data.(L).HW_ms_byBin_smoothed),':k','LineWidth',1);
        hold on
        plot(cfg.spk.phase_bin_centers, ...
            100*(1+data.(L).HW_MI(1)*cos(cfg.spk.phase_bin_centers-data.(L).HW_MI(3))), ...
            'Color', col,'LineWidth',2);
        plot(cfg.spk.phase_bin_centers, ...
            100*data.(L).HW_ms_byBin_smoothed / nanmean(data.(L).HW_ms_byBin_smoothed), ...
            '--', 'Color', col,'LineWidth',2)
        title(['Motion Index, %: ' num2str(100 * data.(L).HW_MI(1)) ...
            ', p = ' num2str(data.(L).HW_MI(2))])
        xlabel('Heart-Cycle Phase (0-2pi)')
        ylabel({'% HW Signal Change,', 'signal / average'})
        % adjust yaxis
        if ~isnan(data.(L).spike_phases_radians)  
            ylims = get(gca, 'YLim');
            ylims(1) = ylims(1) - (ylims(2) - ylims(1));
            set(gca, 'YLim', ylims)
        end
        
        subplot(3,5,10+c)
        yyaxis left
        set(gca,'YColor',col/2);
        plot(data.(L).spike_phases_radians, data.(L).HW_ms, '.', 'Color', col/2)
        % adjust yaxis
        ylims = get(gca, 'YLim');
        ylims(2) = ylims(2) + (ylims(2) - ylims(1));
        set(gca, 'YLim', ylims)
        ylabel('Spike TPW, \muV')
        
        yyaxis right
        set(gca,'YColor','k');
        plot(cfg.spk.phase_bin_centers, 100*data.(L).TPW_ms_byBin'/nanmean(data.(L).TPW_ms_byBin_smoothed),':k','LineWidth',1);
        hold on
        plot(cfg.spk.phase_bin_centers, ...
            100*(1+data.(L).TPW_MI(1)*cos(cfg.spk.phase_bin_centers-data.(L).TPW_MI(3))), ...
            'Color', col,'LineWidth',2);
        plot(cfg.spk.phase_bin_centers, ...
            100*data.(L).TPW_ms_byBin_smoothed / nanmean(data.(L).TPW_ms_byBin_smoothed), ...
            '--', 'Color', col,'LineWidth',2)
        title(['Motion Index, %: ' num2str(100 * data.(L).TPW_MI(1)) ...
            ', p = ' num2str(data.(L).TPW_MI(2))])
        xlabel('Heart-Cycle Phase (0-2pi)')
        ylabel({'% TPW Signal Change,', 'signal / average'})
        % adjust yaxis
        if ~isnan(data.(L).spike_phases_radians)  
            ylims = get(gca, 'YLim');
            ylims(1) = ylims(1) - (ylims(2) - ylims(1));
            set(gca, 'YLim', ylims)
        end
        
        subplot(3,5,13+c)
        yyaxis left
        set(gca,'YColor',col/2);
        plot(data.(L).spike_phases_radians, data.(L).REP_ms, '.', 'Color', col/2)
        % adjust yaxis
        ylims = get(gca, 'YLim');
        ylims(2) = ylims(2) + (ylims(2) - ylims(1));
        set(gca, 'YLim', ylims)
        ylabel('Spike REP, \muV')
        
        yyaxis right
        set(gca,'YColor','k');
        plot(cfg.spk.phase_bin_centers, 100*data.(L).REP_ms_byBin'/nanmean(data.(L).REP_ms_byBin_smoothed),':k','LineWidth',1);
        hold on
        plot(cfg.spk.phase_bin_centers, ...
            100*(1+data.(L).REP_MI(1)*cos(cfg.spk.phase_bin_centers-data.(L).REP_MI(3))), ...
            'Color', col,'LineWidth',2);
        plot(cfg.spk.phase_bin_centers, ...
            100*data.(L).REP_ms_byBin_smoothed / nanmean(data.(L).REP_ms_byBin_smoothed), ...
            '--', 'Color', col,'LineWidth',2)
        title(['Motion Index, %: ' num2str(100 * data.(L).REP_MI(1)) ...
            ', p = ' num2str(data.(L).REP_MI(2))])
        xlabel('Heart-Cycle Phase (0-2pi)')
        ylabel({'% REP Signal Change,', 'signal / average'})
        % adjust yaxis
        if ~isnan(data.(L).spike_phases_radians)  
            ylims = get(gca, 'YLim');
            ylims(1) = ylims(1) - (ylims(2) - ylims(1));
            set(gca, 'YLim', ylims)
        end
    end
    
    filename= ['Cosine_Fitted__' data.unitId, '_' data.target];
    export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    %% results of correlation analysis
    f5 = figure;
    set(f5, 'Position', [381 424 1450 398])
    if matlab_year >= 2020
        sgtitle(sgtitleText, 'interpreter', 'none')
    end
    
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        subplot(1,numel(cfg.condition)+2,c)
        scatter(data.(L).FRbyRR_Hz, data.(L).cycleDurations_s, [], condition_colors{c}(1:3), 'Marker', '.')
        hold on
        linTrend = lsline(gca);
        linTrend.Color = condition_colors{c}(1:3);
        xlabel('Firing Rate per Heart Cycle, Hz')
        ylabel('Heart Cycle Duration, s')
        box on
        title([L ': cc = ' num2str(data.(L).pearson_r(4)) '; p = ' num2str(data.(L).permuted_p(4))])
        legend({'Real Data', 'Least-Square Fit'}, 'Location', 'Best')
    end
    
    subplot(1,numel(cfg.condition)+2,numel(cfg.condition)+2)
    box on
    hold on
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        plot(data.cc_lag_list, data.(L).pearson_r, '-o', 'Color', condition_colors{c}(1:3))
    end
    ylim([-0.4 0.4])
    xlabel('Lag: Number of Heart Cycles Shifted')
    ylabel('Correlation Coefficient')
    legend({cfg.condition.name})
    
    filename= ['FR_RR_Correlations__' data.unitId, '_' data.target];
    export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    %% unit FR per heart cycle as a function of time
	figure;
    set(gcf, 'Position', [1 41 1920 482])

    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        yyaxis left
        plot(data.(L).timeRRstart, data.(L).FRbyRR_Hz, 'o-', 'Color', cfg.condition(c).color)
        ylabel('FR per Heart Cycle, Hz')
        yyaxis right
        plot(data.(L).timeRRstart, data.(L).cycleDurations_s, 'kx--')
        ylabel('Heart-Cycle Duration, s')
    end
    xlabel('Time from Session Start, s')
    if matlab_year >= 2020
        sgtitle(sgtitleText, 'interpreter', 'none')
    end
    
    filename= ['Time_vs_FRperCycle_' data.unitId, '_' data.target];
    export_fig(gcf,[dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    close all
    
end
