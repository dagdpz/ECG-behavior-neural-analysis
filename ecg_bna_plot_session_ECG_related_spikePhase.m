function ecg_bna_plot_session_ECG_related_spikePhase(session_info,ecg_bna_cfg)

dataFolder = [session_info.SPK_fldr filesep 'cardioballistic' filesep];

Fs = 2.441406250000000e+04; % sampling frequency of BB signal, Hz
wf_times_ms = 1000 * (1/Fs:1/Fs:32/Fs); % in ms
wf_times_interp_ms = 1000 * (1/4/Fs:1/4/Fs:32/Fs); % in ms
peak_id = 10; % sample number of the trough in the spike waveform
phase_bins = linspace(0, 2*pi, 64);
% phase_bins = phase_bins(1:end-1);

condition_labels={'Rest','Task'};
condition_colors={[0 0 1 0.1],[1 0 0 0.1]};

% for ii = 1:116
%     
%     unitNum = ['00' num2str(ii)];
%     if ii < 100
%         unitNum = unitNum(end-1:end);
%     else
%         unitNum = unitNum(end-2:end);
%     end
%     
%     dataFile = [saveFolder 'Mag_20230511_' unitNum '_spikes_ECGphase.mat'];
%     if exist(dataFile,'file')
%         load(dataFile, 'data')
%     else
%         recordingArea{ii}  = '_';
%         MI_AMP(ii)         = NaN;
%         MI_HW(ii)          = NaN;
%         MI_TPW(ii)         = NaN;
%         MI_REP(ii)         = NaN;
%         pearson_r(ii)      = NaN;
%         continue
%     end
%     
%     recordingArea{ii}  = data.target;
%     MI_AMP(ii)         = data.AMP_MI(1);
%     MI_HW(ii)          = data.HW_MI(1);
%     MI_TPW(ii)         = data.TPW_MI(1);
%     MI_REP(ii)         = data.REP_MI(1);
%     pearson_r(ii)      = data.pearson_r;
%     
% end
% 
% unqAreas = unique(recordingArea);
% for ii = 1:length(unqAreas)
%     currArea = unqAreas{ii};
%     idx(ii,:) = cellfun(@ (x) strcmp(x, currArea), recordingArea);
%     
%     MI_AMP_bins(ii,:) = histcounts(100 * MI_AMP(idx(ii,:)), 0:10);
%     MI_HW_bins(ii,:)  = histcounts(100 * MI_HW(idx(ii,:)), 0:10);
%     MI_TPW_bins(ii,:) = histcounts(100 * MI_TPW(idx(ii,:)), 0:10);
%     MI_REP_bins(ii,:) = histcounts(100 * MI_REP(idx(ii,:)), 0:10);
% end
% 
% f0 = figure;
% set(f0, 'Position', [784   148   942   797])
% subplot(2,2,1)
% bar(MI_AMP_bins','stacked')
% title('AMP')
% xlabel('MI, % change')
% ylabel('Unit Count')
% legend(unqAreas)
% 
% subplot(2,2,2)
% bar(MI_HW_bins','stacked')
% title('HW')
% 
% subplot(2,2,3)
% bar(MI_TPW_bins','stacked')
% title('TPW')
% 
% subplot(2,2,4)
% bar(MI_REP_bins', 'stacked')
% title('REP')
% 
% f1 = figure;
% set(f1, 'Position', [784   148   942   797])
% subplot(2,2,1)
% hist(MI_AMP * 100, [0:10])
% xlim([0 10])
% xlabel('AMP Motion Index, % change')
% ylabel('Unit Count')
% 
% subplot(2,2,2)
% hist(MI_HW * 100, [0:10])
% xlim([0 10])
% xlabel('HW Motion Index, % change')
% ylabel('Unit Count')
% 
% subplot(2,2,3)
% hist(MI_TPW * 100, [0:10])
% xlim([0 10])
% xlabel('TPW Motion Index, % change')
% ylabel('Unit Count')
% 
% subplot(2,2,4)
% hist(MI_REP * 100, [0:10])
% xlim([0 10])
% xlabel('REP Motion Index, % change')
% ylabel('Unit Count')
% 
% f2 = figure;
% set(f2, 'Position', [784   148   942   797])
% hist(pearson_r)
% xlim([-0.5 0.5])
% xlabel('Pearson Correlation Coefficient')
% ylabel('Unit Count')

for ii = 1:115
    
    unitNum = ['00' num2str(ii)];
    if ii < 100
        unitNum = unitNum(end-1:end);
    else
        unitNum = unitNum(end-2:end);
    end
    
    dataFile = [dataFolder 'Mag_20230511_' unitNum '_spikes_ECGphase.mat'];
    if exist(dataFile,'file')
        load(dataFile, 'data')
    else
        continue
    end
    % with the current Fs and number of samples per spike I have 1.4 ms per
    % spike
    

    
    %% AMP - absolute magnitude at the trough
    %     %% HW - in ms corresponding to the total time that the EAP trough is below half AMP
    %     HAMP = AMP_abs/2; % half of max amplitude
    %     beyond_halfWidth_idx = bsxfun(@(x,y) sign_trough*x>y, Y(15:40,:), HAMP); % multiply by signum: if spike is negative it will flip, if positiv - not
    % %     figure, imagesc(beyond_halfWidth_idx')
    %     HW = sum(beyond_halfWidth_idx,1)*(1/3/Fs); % in seconds
    %
    %     HWbyPhase      = arrayfun(@(x) mean(HW(bin == x)), unique(bin));
    %     HWbyPhase_std  = arrayfun(@(x) std(HW(bin == x)), unique(bin));
    
    %%
    sgtitleText = {[data.unitId '_' data.target ' ch ' num2str(data.channel) ';'  ' unit ...'], ...
        ['FR, Hz: ' num2str(data.FR) '; SNR: ' num2str(data.SNR) '; Fano Factor: ' num2str(data.stability) '; % of ISIs < 3 ms: '  num2str(100 * data.single_rating) '%']};
    % ['AMP_MI = ' num2str(data.AMP_MI(1)) '; p = ' num2str(data.Task.AMP_MI(2))]
    
    %% Overall picture of waveforms and PSTHs
    f1 = figure;
    set(f1, 'Position', [784   148   942   797])
    sgtitle(sgtitleText, 'interpreter', 'none')
    
    ax_sp = subplot(2,2,1);
    box on
    hold on
    for tasktype = 1:2
        line(wf_times_interp_ms, data.(condition_labels{tasktype}).waveforms_upsampled_microvolts(1:300:end, :), 'Color', condition_colors{tasktype})
        line(repmat([wf_times_interp_ms(1) wf_times_interp_ms(end)],4,1)', repmat(data.(condition_labels{tasktype}).thresholds_microV,1,2)', 'Color', 'r')
    end
    xlim([wf_times_interp_ms(1) wf_times_interp_ms(end)])
    xlabel('Time, ms')
    ylabel('Voltage, μV')
    title({'Example Waveforms: ', ...
        ['\color{blue}Rest: ' num2str(size(data.Rest.waveforms_upsampled_microvolts(1:300:end,:),1)) ' out of ' num2str(size(data.Rest.waveforms_upsampled_microvolts,1))], ...
        ['\color{red}Task: ' num2str(size(data.Task.waveforms_upsampled_microvolts(1:300:end,:),1)) ' out of ' num2str(size(data.Task.waveforms_upsampled_microvolts,1))]})
    
    subplot(2,2,2)
    box on
    for tasktype = 1:2
        h = histcounts(data.(condition_labels{tasktype}).spike_phases_radians, 64);
        h = h / mean(h);
        polarplot(phase_bins,h, 'Color', condition_colors{tasktype}(1:3))
        hold on
    end
    hold off
    title({'ECG-triggered polar PSTH (spike counts / mean spike counts per bin)', ...
        '\color{blue}Rest', '\color{red}Task'})
    axis tight
    
    subplot(2,2,3)
    title('Average WFs by bin: \color{blue}Rest and \color{red}Task')
    hold on
    for tasktype = 1:2
        line(wf_times_interp_ms, data.(condition_labels{tasktype}).waveforms_byBin_microvolts, 'Color', condition_colors{tasktype})
    end
    hold off
    box on
    xlabel('Spike Time, ms')
    ylabel('Voltage, μV')
    set(gca, 'XLim', get(ax_sp,'XLim'), 'YLim', get(ax_sp,'YLim'))
    
    subplot(2,2,4)
    imagesc(data.Rest.waveforms_byBin_microvolts - mean(data.Rest.waveforms_byBin_microvolts,1))
    cbar = colorbar;
    cbar.Title.String = '\Delta Spike Amplitude, μV';
    caxis([-3 3])
    title('{\color{blue}[Only Rest]} WF Voltage Change over Heart Cycle')
    xlabel('Spike Time')
    ylabel('Heart-Cycle Phase')
    set(gca, 'XTickLabel', [], 'YTickLabel', [])
    
    filename= ['Spikes_PolarPSTH__' data.unitId, '_' data.target];
    export_fig([dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    %% Distributions of features
    f2 = figure;
    set(f2, 'Position', [784   148   942   797])
    sgtitle(sgtitleText, 'interpreter', 'none')
    
    subplot(2,2,1)
    box on
%     b = [histcounts(abs(data.Rest.AMP_microV), [0:15:350]); ...
%         histcounts(abs(data.Task.AMP_microV), [0:15:350])];
%     b = ( b ./ max(b,[],2) )';
%     barplot = bar([7.5:15:350], b);
%     for tasktype = 1:2
%         barplot(tasktype).FaceColor = condition_colors{tasktype}(1:3);
%     end
    hold on
    for tasktype = 1:2
        histogram(abs(data.(condition_labels{tasktype}).AMP_microV), [0:5:350], 'FaceColor', condition_colors{tasktype}(1:3))
    end
    hold off
    xline(data.(condition_labels{tasktype}).thresholds_microV(1), 'Color', 'r')
    xline(data.(condition_labels{tasktype}).thresholds_microV(2), 'Color', 'r')
    xlabel('AMP, μV');
    ylabel('Spike Counts')
    xlim([0 350])
    hold off
    title('AMP: Same X-axis for All Units')
    legend(condition_labels, 'Location', 'Best')
    
    subplot(2,2,2)
    box on
    hold on
    for tasktype = 1:2
        histogram(abs(data.(condition_labels{tasktype}).HW_ms), [0:0.005:0.5], 'FaceColor', condition_colors{tasktype}(1:3))
    end
    hold off
    xlim([0 0.5])
    xlabel('HW, ms');
    ylabel('Spike Counts')
    title('HW: Same X-axis for All Units')
    
    subplot(2,2,3)
    box on
    hold on
    for tasktype = 1:2
        histogram(abs(data.(condition_labels{tasktype}).TPW_ms), 100, 'FaceColor', condition_colors{tasktype}(1:3))
    end
    hold off
    xlabel('TPW, ms');
    ylabel('Spike Counts')
    title('TPW: X-axis Adjusts for Each Unit')
    
    subplot(2,2,4)
    box on
    hold on
    for tasktype = 1:2
        histogram(abs(data.(condition_labels{tasktype}).REP_ms), 100, 'FaceColor', condition_colors{tasktype}(1:3))
    end
    hold off
    xlabel('REP, ms');
    ylabel('Spike Counts')
    title('REP: X-axis Adjusts for Each Unit')
    
    filename= ['Distributions_AMP_HW_TPW_REP__' data.unitId, '_' data.target];
    export_fig([dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    %% Real and reshuffled data
    f3 = figure;
    set(f3, 'Position', [1 41 1920 963])
    sgtitle(sgtitleText, 'interpreter', 'none')
    
    for tasktype = 1:2
        subplot(3,5,tasktype)
        yyaxis left
        histogram(abs(data.(condition_labels{tasktype}).spike_phases_radians), phase_bins, 'FaceColor', condition_colors{tasktype}(1:3))
        ylabel('Spike Counts')
        yyaxis right
        filledArea = fill([phase_bins(1:end-1) fliplr(phase_bins(1:end-1)) phase_bins(1)], ...
            [data.(condition_labels{tasktype}).AMP_upperPrctile_97_5 fliplr(data.(condition_labels{tasktype}).AMP_lowerPrctile_2_5) data.(condition_labels{tasktype}).AMP_upperPrctile_97_5(1)], ...
            [0 0 0], 'FaceALpha', 0.5, 'EdgeColor', 'none');
        hold on;
        p1 = plot(phase_bins(1:end-1), data.(condition_labels{tasktype}).AMP_microV_byBin,'-k','LineWidth',2);
        p2 = plot(phase_bins(1:end-1), data.(condition_labels{tasktype}).AMP_microV_byBin_smoothed, '-', 'Color', 'w', 'LineWidth',2);
        title(['AMP: MI = ' num2str(data.(condition_labels{tasktype}).AMP_MI(1)) '; p = ' num2str(data.(condition_labels{tasktype}).AMP_MI(2))])
        xlim([0 2*pi])
        xlabel('Heart-cycle Phase (0-2pi)')
        ylabel('AMP, microvolts')
        if tasktype == 1
            legend([filledArea p1 p2], {'95% Confidence Interval', 'Average by Bin', 'Rlowess-Smoothed'}, 'Location', 'southeast')
        end
    end
    
    for tasktype = 1:2
        subplot(3,5,3+tasktype)
        yyaxis left
        histogram(data.(condition_labels{tasktype}).spike_phases_radians, phase_bins, 'FaceColor', condition_colors{tasktype}(1:3))
        ylabel('Spike Counts')
        yyaxis right
        fill([phase_bins(1:end-1) fliplr(phase_bins(1:end-1)) phase_bins(1)], ...
            [data.(condition_labels{tasktype}).HW_upperPrctile_97_5 fliplr(data.(condition_labels{tasktype}).HW_lowerPrctile_2_5) data.(condition_labels{tasktype}).HW_upperPrctile_97_5(1)], ...
            [0 0 0], 'FaceALpha', 0.5, 'EdgeColor', 'none');
        hold on;
        plot(phase_bins(1:end-1), data.(condition_labels{tasktype}).HW_ms_byBin,'-k','LineWidth',2);
        plot(phase_bins(1:end-1), data.(condition_labels{tasktype}).HW_ms_byBin_smoothed, '-', 'Color', 'w', 'LineWidth',2);
        title(['HW: MI = ' num2str(data.(condition_labels{tasktype}).HW_MI(1)) '; p = ' num2str(data.(condition_labels{tasktype}).HW_MI(2))])
        xlim([0 2*pi])
        xlabel('Heart-cycle Phase (0-2pi)')
        ylabel('HW, ms')
    end
    
    for tasktype = 1:2
        subplot(3,5,10+tasktype)
        yyaxis left
        histogram(data.(condition_labels{tasktype}).spike_phases_radians, phase_bins, 'FaceColor', condition_colors{tasktype}(1:3))
        ylabel('Spike Counts')
        yyaxis right
        fill([phase_bins(1:end-1) fliplr(phase_bins(1:end-1)) phase_bins(1)], ...
            [data.(condition_labels{tasktype}).TPW_upperPrctile_97_5 fliplr(data.(condition_labels{tasktype}).TPW_lowerPrctile_2_5) data.(condition_labels{tasktype}).TPW_upperPrctile_97_5(1)], ...
            [0 0 0], 'FaceALpha', 0.5, 'EdgeColor', 'none');
        hold on;
        plot(phase_bins(1:end-1), data.(condition_labels{tasktype}).TPW_ms_byBin,'-k','LineWidth',2);
        plot(phase_bins(1:end-1), data.(condition_labels{tasktype}).TPW_ms_byBin_smoothed, '-', 'Color', 'w', 'LineWidth',2)
        title(['TPW: MI = ' num2str(data.(condition_labels{tasktype}).TPW_MI(1)) '; p = ' num2str(data.(condition_labels{tasktype}).TPW_MI(2))])
        xlim([0 2*pi])
        xlabel('Heart-cycle Phase (0-2pi)')
        ylabel('TPW, ms')
    end
    
    for tasktype = 1:2
        subplot(3,5,13+tasktype)
        yyaxis left
        histogram(data.(condition_labels{tasktype}).spike_phases_radians, phase_bins, 'FaceColor', condition_colors{tasktype}(1:3))
        ylabel('Spike Counts')
        yyaxis right
        fill([phase_bins(1:end-1) fliplr(phase_bins(1:end-1)) phase_bins(1)], ...
            [data.(condition_labels{tasktype}).REP_upperPrctile_97_5 fliplr(data.(condition_labels{tasktype}).REP_lowerPrctile_2_5) data.(condition_labels{tasktype}).REP_upperPrctile_97_5(1)], ...
            [0 0 0], 'FaceALpha', 0.5, 'EdgeColor', 'none')
        hold on
        plot(phase_bins(1:end-1), data.(condition_labels{tasktype}).REP_ms_byBin,'-k','LineWidth',2);
        plot(phase_bins(1:end-1), data.(condition_labels{tasktype}).REP_ms_byBin_smoothed, '-', 'Color', 'w', 'LineWidth',2)
        title(['AMP: MI = ' num2str(data.(condition_labels{tasktype}).REP_MI(1)) '; p = ' num2str(data.(condition_labels{tasktype}).REP_MI(2))])
        xlim([0 2*pi])
        xlabel('Heart-cycle Phase (0-2pi)')
        ylabel('REP, ms')
    end
    
    filename= ['PSTH_overlaid_Feature_Dynamics__' data.unitId, '_' data.target];
    export_fig([dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    %% cosine fitted plots
    f4 = figure;
    set(f4, 'Position', [1 41 1920 963])
    sgtitle(sgtitleText, 'interpreter', 'none')
    
    for tasktype = 1:2
        subplot(3,5,tasktype)
        fitCardiacModulation(phase_bins(1:end-1), data.(condition_labels{tasktype}).AMP_microV_byBin', {[condition_labels{tasktype} ' AMP']}, 1, 263);
        
        subplot(3,5,3+tasktype)
        fitCardiacModulation(phase_bins(1:end-1), data.(condition_labels{tasktype}).HW_ms_byBin', {[condition_labels{tasktype} ' HW']}, 1, 264);
        
        subplot(3,5,10+tasktype)
        fitCardiacModulation(phase_bins(1:end-1), data.(condition_labels{tasktype}).TPW_ms_byBin', {[condition_labels{tasktype} ' TPW']}, 1, 269);
        
        subplot(3,5,13+tasktype)
        fitCardiacModulation(phase_bins(1:end-1), data.(condition_labels{tasktype}).REP_ms_byBin', {[condition_labels{tasktype} ' REP']}, 1);
    end
    
    filename= ['Cosine_Fitted__' data.unitId, '_' data.target];
    export_fig([dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);

    close all
    
end
