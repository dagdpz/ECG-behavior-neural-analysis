function ecg_bna_plot_PSD_LFP(all_sites, fs)

clear legend_list
h = figure;
hold on
for site = 1: length(all_sites)
    lfp_matrix = all_sites(site,:);
    time = (0:size(lfp_matrix,2)-1) / fs; % Create time vector
    % Initialize FieldTrip data structure
    data = [];
    data.label = arrayfun(@(x) sprintf('LFP%02d', x), 1:size(lfp_matrix, 1), 'UniformOutput', false);
    data.fsample = fs;  % Sampling frequency
    data.time = {time}; % Cell array containing time vector
    data.trial = {lfp_matrix}; % Cell array containing LFP data matrix
    data.sampleinfo = [1 size(lfp_matrix,2)]; % Start and end sample indices
    
    ft_datatype_raw(data) % To Check if the structure is correctly formatted
    ft_data = data;
    
    cfg1 = [];
    cfg1.length  = 2;
    cfg1.overlap = 0.5;
    base_rpt2    = ft_redefinetrial(cfg1, ft_data);
    
    cfg2 = [];
    cfg2.output  = 'pow';
    cfg2.channel =  1; %'all';
    cfg2.method  = 'mtmfft';
    cfg2.taper   = 'boxcar';
    cfg2.foi     = 1:1:120; % 1/cfg1.length  = 1;
    base_freq2   = ft_freqanalysis(cfg2, base_rpt2);
    
    % %     subplot(8,8,site)
    plot(base_freq2.freq, mean(base_freq2.powspctrm(:,:),1))
    legend('2 sec window')
    xlabel('Frequency (Hz)');
    ylabel('absolute power (uV^2)');
    legend_list(site) = {['Sig-', num2str(site)]};
end
hold off

legend(legend_list, 'Location','eastoutside','Interpreter','none');
legend_list= legend_list';
h.WindowState ='maximized';
set(gca,'xlim',[0 20]);
end