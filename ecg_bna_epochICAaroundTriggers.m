function [ica_epoched, avg] = ecg_bna_epochICAaroundTriggers(saving_path,session,comp, trigger_samples, pre, post, doICAPlot,doFFTplots)
% Epoch ICA components around trigger points
% 
% Inputs:
%   comp            - FieldTrip ICA structure (from ft_componentanalysis)
%   trigger_samples - Vector of trigger indices (in same fs as comp)
%   pre             - Time before trigger in seconds (e.g., 0.5)
%   post            - Time after trigger in seconds (e.g., 0.5)
%   doPlot          - Boolean flag to plot averaged components
%
% Outputs:
%   ica_epoched     - FieldTrip structure of epoched ICA data
%   avg             - Timelocked average (FieldTrip structure)

    fs = comp.fsample;
    samples_pre  = round(pre * fs);
    samples_post = round(post * fs);

    trl = zeros(length(trigger_samples), 3);
    for i = 1:length(trigger_samples)
        trl(i,1) = trigger_samples(i) - samples_pre;
        trl(i,2) = trigger_samples(i) + samples_post - 1;
        trl(i,3) = -samples_pre;  % offset (time-zero)
    end

    % Epoch ICA data
    cfg = [];
    cfg.trl = trl;
    ica_epoched = ft_redefinetrial(cfg, comp);

    % Timelock average
    cfg = [];
    cfg.channel = 'all';
    avg = ft_timelockanalysis(cfg, ica_epoched);
    
%     if doPlot
% % %         cfg = [];
% % %         cfg.layout = 'component';
% % %         cfg.showlabels = 'yes';
% % %         cfg.fontsize = 10;
% % %         ft_multiplotER(cfg, avg);
% % %         sgtitle('ICA Components Averaged Around Triggers');
%         cfg = [];
%         cfg.channel = 'all';
%         ft_singleplotER(cfg, avg);
%         title('ICA Components Averaged Around Triggers');
%     end
    if doICAPlot
        fig = figure;
        for cmp = 1: length(comp.label)
            subplot(6,6,cmp)
            plot(avg.time,avg.avg(cmp,:),'LineWidth',1.5)
            xline(0,'k--')
            xlabel('Time(sec)')
            title(['ICA Components: ',num2str(cmp),' Averaged Around Triggers'],'Interpreter','none','FontSize',6);
        end
        sgtitle(session,'Interpreter','none')
        fig.WindowState = 'maximized';
        figName2save = fullfile([saving_path,filesep,session,'_R triggered']);
        saveas(fig,figName2save,'fig');
        saveas(fig,figName2save,'png');
    end
    
    
    Fs = fs;              % Sampling frequency
    T = 1/Fs;             % Sampling period
    L = length(avg.time); % Length of signal
    t = avg.time  ;       % Time vector
    f = Fs*(0:(L/2))/L;   % frequency domain
    if doFFTplots
         fig = figure;
        for cmp = 1: length(comp.label)
            Y = fft(avg.avg(cmp,:));
            P2 = abs(Y/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            
            subplot(6,6,cmp)
            plot(f,P1)
            xlabel('f (Hz)')
            ylabel('|P1(f)|')
            set(gca,'xlim',[0 60])
            title(['ICA Components: ',num2str(cmp)],'Interpreter','none','FontSize',6);
        end
        sgtitle([session,' Single-Sided Amplitude Spectrum'],'Interpreter','none')
        fig.WindowState = 'maximized';
        figName2save = fullfile([saving_path,filesep,session,'_R triggered FFT']);
        saveas(fig,figName2save,'fig');
        saveas(fig,figName2save,'png');
    end
end
