function ecg_bna_plot_session_ECG_related_spikePhase(session_info,ecg_bna_cfg)

dataFolder = [session_info.SPK_fldr filesep 'cardioballistic' filesep];

Fs = 2.441406250000000e+04; % sampling frequency of BB signal, Hz
wf_times_ms = 1000 * (1/Fs:1/Fs:32/Fs); % in ms
wf_times_interp_ms = 1000 * (1/4/Fs:1/4/Fs:32/Fs); % in ms
peak_id = 10; % sample number of the trough in the spike waveform
phase_bins = linspace(0, 2*pi, 63);
% phase_bins = phase_bins(1:end-1);

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
    sgtitleText = {[data.unitId '_' data.target], ...
        ['FR, Hz: ' num2str(data.FR) '; SNR: ' num2str(data.SNR) '; Fano Factor: ' num2str(data.stability) '; % of ISIs < 3 ms: '  num2str(100 * data.single_rating) '%'], ...
        ['AMP_MI = ' num2str(data.AMP_MI(1)) '; p = ' num2str(data.AMP_MI(2))]};
    
    %% Overall picture of waveforms and PSTHs
    f1 = figure;
    set(f1, 'Position', [784   148   942   797])
    sgtitle(sgtitleText, 'interpreter', 'none')
    
    ax_sp = subplot(2,2,1);
    box on
    hold on
    line(wf_times_interp_ms, data.waveforms_upsampled_microvolts(1:300:end, :), 'Color', [0 0 1 0.1])
    line(repmat([wf_times_interp_ms(1) wf_times_interp_ms(end)],4,1)', repmat(data.thresholds_microV,1,2)', 'Color', 'r')
    xlim([wf_times_interp_ms(1) wf_times_interp_ms(end)])
    xlabel('Time, ms')
    ylabel('Voltage, μV')
    title({'Example Waveforms: ', [num2str(size(data.waveforms_upsampled_microvolts(1:300:end,:),1)) ' out of ' num2str(size(data.waveforms_upsampled_microvolts,1))]})
    
    subplot(2,2,2)
    box on
    polarhistogram(data.spike_phases_radians, phase_bins)
    title('ECG-triggered polar PSTH (spike counts)')
    axis tight
    
    subplot(2,2,3)
    
    title('Placeholder of the Average WFs by bin')
    plot(wf_times_interp_ms, data.waveforms_byBin_microvolts)
    xlabel('EAP Time, ms')
    ylabel('Voltage, μV')
    set(gca, 'XLim', get(ax_sp,'XLim'), 'YLim', get(ax_sp,'YLim'))
    
    subplot(2,2,4)
    title('WF Voltage Change over Heart Cycle')
    imagesc(data.waveforms_byBin_microvolts - mean(data.waveforms_byBin_microvolts,1))
    cbar = colorbar;
    cbar.Title.String = '\Delta Spike Amplitude, μV';
    caxis([-3 3])
    xlabel('EAP time')
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
    bar([5:5:350], histcounts(abs(data.AMP_microV), [0:5:350]))
    hold on
    xline(data.thresholds_microV(1), 'Color', 'r')
    xline(data.thresholds_microV(2), 'Color', 'r')
    xlabel('AMP, μV');
    ylabel('Spike Counts')
    xlim([0 350])
    hold off
    
    subplot(2,2,2)
    histogram(abs(data.HW_ms), 100)
    xlim([0 0.5])
    xlabel('HW, ms');
    ylabel('Spike Counts')
    
    subplot(2,2,3)
    histogram(abs(data.TPW_ms), 100)
    hold on
    xlabel('TPW, ms');
    ylabel('Spike Counts')
    
    subplot(2,2,4)
    histogram(abs(data.REP_ms), 100)
    hold on
    xlabel('REP, ms');
    ylabel('Spike Counts')
    
    filename= ['Distributions_AMP_HW_TPW_REP__' data.unitId, '_' data.target];
    export_fig([dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    %% Real and reshuffled data
    f3 = figure;
    set(f3, 'Position', [784   148   942   797])
    sgtitle(sgtitleText, 'interpreter', 'none')
    
    subplot(2,2,1)
    yyaxis left
    histogram(data.spike_phases_radians, phase_bins)
    ylabel('Spike Counts')
    yyaxis right
    filledArea = fill([phase_bins fliplr(phase_bins) phase_bins(1)], ...
        [data.AMP_upperPrctile_97_5 fliplr(data.AMP_lowerPrctile_2_5) data.AMP_upperPrctile_97_5(1)], ...
        'r', 'FaceALpha', 0.3, 'EdgeColor', 'none');
    hold on;
    p1 = plot(phase_bins, data.AMP_microV_byBin,'-k','LineWidth',2);
    p2 = plot(phase_bins, data.AMP_microV_byBin_smoothed, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth',2);
    title(['AMP: MI = ' num2str(data.AMP_MI(1)) '; p = ' num2str(data.AMP_MI(2))])
    xlim([0 2*pi])
    xlabel('Heart-cycle Phase (0-2pi)')
    ylabel('AMP, microvolts')
    legend([filledArea p1 p2], {'95% Confidence Interval', 'Average by Bin', 'Rlowess-Smoothed'}, 'Location', 'southeast')
    
    subplot(2,2,2)
    yyaxis left
    histogram(data.spike_phases_radians, phase_bins)
    ylabel('Spike Counts')
    yyaxis right
    fill([phase_bins fliplr(phase_bins) phase_bins(1)], ...
        [data.HW_upperPrctile_97_5 fliplr(data.HW_lowerPrctile_2_5) data.HW_upperPrctile_97_5(1)], ...
        'r', 'FaceALpha', 0.2, 'EdgeColor', 'none');
    hold on;
    plot(phase_bins, data.HW_ms_byBin,'-k','LineWidth',2);
    plot(phase_bins, data.HW_ms_byBin_smoothed, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth',2);
    title(['HW: MI = ' num2str(data.HW_MI(1)) '; p = ' num2str(data.HW_MI(2))])
    xlim([0 2*pi])
    xlabel('Heart-cycle Phase (0-2pi)')
    ylabel('HW, ms')
    
    subplot(2,2,3)
    yyaxis left
    histogram(data.spike_phases_radians, phase_bins)
    ylabel('Spike Counts')
    yyaxis right
    fill([phase_bins fliplr(phase_bins) phase_bins(1)], ...
        [data.TPW_upperPrctile_97_5 fliplr(data.TPW_lowerPrctile_2_5) data.TPW_upperPrctile_97_5(1)], ...
        'r', 'FaceALpha', 0.2, 'EdgeColor', 'none');
    hold on;
    plot(phase_bins, data.TPW_ms_byBin,'-k','LineWidth',2);
    plot(phase_bins, data.TPW_ms_byBin_smoothed, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth',2)
    title(['TPW: MI = ' num2str(data.TPW_MI(1)) '; p = ' num2str(data.TPW_MI(2))])
    xlim([0 2*pi])
    xlabel('Heart-cycle Phase (0-2pi)')
    ylabel('TPW, ms')
    
    subplot(2,2,4)
    yyaxis left
    histogram(data.spike_phases_radians, phase_bins)
    ylabel('Spike Counts')
    yyaxis right
    fill([phase_bins fliplr(phase_bins) phase_bins(1)], ...
        [data.REP_upperPrctile_97_5 fliplr(data.REP_lowerPrctile_2_5) data.REP_upperPrctile_97_5(1)], ...
        'r', 'FaceALpha', 0.2, 'EdgeColor', 'none')
    hold on
    plot(phase_bins, data.REP_ms_byBin,'-k','LineWidth',2);
    plot(phase_bins, data.REP_ms_byBin_smoothed, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth',2)
    title(['AMP: MI = ' num2str(data.REP_MI(1)) '; p = ' num2str(data.REP_MI(2))])
    xlim([0 2*pi])
    xlabel('Heart-cycle Phase (0-2pi)')
    ylabel('REP, ms')
    
    filename= ['PSTH_overlaid_Feature_Dynamics__' data.unitId, '_' data.target];
    export_fig([dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
    
    %% cosine fitted plots
    f4 = figure;
    set(f4, 'Position', [784   148   942   797])
    sgtitle(sgtitleText, 'interpreter', 'none')
    
    subplot(2,2,1)
    fitCardiacModulation(phase_bins, data.AMP_microV_byBin', {'AMP'}, 1, 263);
           
    subplot(2,2,2)
    fitCardiacModulation(phase_bins, data.HW_ms_byBin', {'HW'}, 1, 264);

    subplot(2,2,3)
    fitCardiacModulation(phase_bins, data.TPW_ms_byBin', {'TPW'}, 1, 269);

    subplot(2,2,4)
    fitCardiacModulation(phase_bins, data.REP_ms_byBin', {'REP'}, 1);
    
    filename= ['Cosine_Fitted__' data.unitId, '_' data.target];
    export_fig([dataFolder, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);

    close all
    
end
