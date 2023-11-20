function ecg_bna_plot_session_spike_histogram(session_info, ecg_bna_cfg)

% check the Matlab version
v = version('-release');
v = str2double(v(1:end-1));

histbins=0.2:0.02:0.8; % bins for RR duration histogram
condition_labels = {'Rest', 'Task'};

% find the current datafile
basepath_to_save=[session_info.SPK_fldr filesep 'per_unit'];
if ~exist(basepath_to_save,'dir')
    mkdir(basepath_to_save);
end

fileList = dir([basepath_to_save, filesep, session_info.Monkey(1:3) '_' session_info.Date '*.mat']);

for flNum = 1:length(fileList)
    
    load([fileList(flNum).folder filesep fileList(flNum).name], 'Output')
    
    % figure out the present units and check that they match in Task and Rest
    if ~isequal(Output.target, Output.target)
        error('Numbers and names of units don''t match for task and rest')
    end
    
%     condition_labels = fieldnames(Output);
    condition_colors={'b','r'};
    BINS=(ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000;
    
    unit_ID = Output.unit_ID;
    target = Output.target;
    
    sgtitleText = {[Output.unit_ID '_' Output.target], ... %  
        ['SNR: ' num2str(Output.quantSNR) '; Fano Factor: ' num2str(Output.stability_rating) '; % of ISIs < 3 ms: '  num2str(100 * Output.Single_rating) '%']};
    
    figure; % raster
    set(gcf, 'Position', [641    40   892   956])
    if v < 2016
        suptitle(sgtitleText,'interpreter','none');
    else
        sgtitle(sgtitleText,'interpreter','none')
    end
    
    for tasktype=1:numel(condition_labels)
        subplot(3,2,tasktype)
        hold on
        set(gca, 'XTick', BINS(1:5:end), 'XTickLabel', BINS(1:5:end))
        xlim([ecg_bna_cfg.analyse_states{3} ecg_bna_cfg.analyse_states{4}]*1000)
        ylim([0 100])
        xlabel('Time from R-peak, s')
        ylabel('Number of Plotted Row')
        box on
        L=condition_labels{tasktype};
        col=condition_colors{tasktype};
        if ~isnan(Output.(L).raster)
            % figure out how many R-peaks we have and if < 100, plot
            % everything; if > 100, choose only 100
            if size(Output.(L).raster,1) <= 100
                stepSize = 1; % just step through all the Rpeaks one by one
            else
                % choose min interval of linearly spaced intervals
                % between 1 and the number of Rpeaks
                stepSize = min(diff(floor(linspace(1, size(Output.(L).raster,1), 100))));
            end
            % plot spikes only for 100 Rpeaks or everything if we have
            % fewer
            a = 1; % introduce row counter
            for RpeakNum = 1:stepSize:size(Output.(L).raster,1)
                x = Output.(L).raster(RpeakNum,:);
                line([BINS; BINS], [x; 2*x]+a-1, 'Color', col, 'LineWidth', 1)
                a = a + 1;
            end
            xline(0)
            title({['Raster for ' L ': 10-ms bins, '], ['100 R-peak intervals out of ' num2str(size(Output.(L).raster,1))], '(linearly spaced)'})
        end
        
    end
    
    for tasktype=1:numel(condition_labels)
        subplot(3,2,tasktype+2)
        hold on
        box on
        xlim([ecg_bna_cfg.analyse_states{3} ecg_bna_cfg.analyse_states{4}]*1000)
        set(gca, 'XTick', BINS(1:5:end), 'XTickLabel', BINS(1:5:end))
        L=condition_labels{tasktype};
        col=condition_colors{tasktype};
        lineProps={'color',col,'linewidth',1};
        shadedErrorBar(BINS, Output.(L).SD, Output.(L).SD_SEM,lineProps,1);
        lineProps={'color',col,'linewidth',1,'linestyle',':'};
        shadedErrorBar(BINS, Output.(L).SDP, [Output.(L).SDPCu; Output.(L).SDPCL], lineProps,1);
        xline(0)
        ypos=NaN;
        if Output.(L).sig_sign==-1;
            ypos=min(Output.(L).SD)*-1;
        elseif Output.(L).sig_sign==1;
            ypos=max(Output.(L).SD);
        end
        to_plot=Output.(L).sig;
        to_plot(to_plot==0)=NaN;
        plot(BINS,to_plot*ypos,col,'linewidth',5);
        
        y_lims=get(gca,'ylim');
        text(BINS(10),y_lims(2)-diff(y_lims)*tasktype*1/20, [L ': trials = ' ,num2str(Output.(L).NrTrials), ' events = ' ,num2str(Output.(L).NrEvents) ],'Color',condition_colors{tasktype});
    end
    
    ylabel('Firing rate');
    xlabel('time to Rpeak');
    
    for tasktype=1:numel(condition_labels)
        L=condition_labels{tasktype};
        col=condition_colors{tasktype};
        subplot(3,2,tasktype+4);
        %         hold on
        box on
        yyaxis left
        histogram(Output.(L).Rds,histbins, 'FaceColor', col)
        ylabel('N');
        yyaxis right
        histogram(Output.(L).Rds_perm,histbins, 'FaceColor', [0.5 0.5 0.5])
        ylabel('N');
        xlabel('RR Duration, s');
        title('grey - jittered RRs');
    end
    filename= ['Raster_PSTH_Rintervals_' unit_ID, '__' target];
    export_fig([basepath_to_save, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
end
end