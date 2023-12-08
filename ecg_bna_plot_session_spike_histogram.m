function ecg_bna_plot_session_spike_histogram(session_info, cfg)

% check the Matlab version
v = version('-release');
v = str2double(v(1:end-1));

histbins=0.2:0.02:0.8; % bins for RR duration histogram
N_conditions=numel(cfg.condition);
condition_colors={cfg.condition.color};
% find the current datafile
basepath_to_save=[cfg.SPK_root_results_fldr filesep 'per_unit'];
if ~exist(basepath_to_save,'dir')
    mkdir(basepath_to_save);
end

fileList = dir([basepath_to_save, filesep, session_info.Monkey(1:3) '_' session_info.Date '*.mat']);

for flNum = 1:length(fileList)
    load([basepath_to_save filesep fileList(flNum).name], 'Output')
    
    % figure out the present units and check that they match in Task and Rest
    if ~isequal(Output.target, Output.target)
        error('Numbers and names of units don''t match for task and rest')
    end
    
%     condition_labels = fieldnames(Output);
    BINS=(cfg.analyse_states{1,3}:cfg.spk.PSTH_binwidth:cfg.analyse_states{1,4})*1000;
    
    unit_ID = Output.unit_ID;
    target = Output.target;
    
    sgtitleText = {[Output.unit_ID '_' Output.target], ... %  
        ['SNR: ' num2str(Output.quantSNR) '; Fano Factor: ' num2str(Output.stability_rating) '; % of ISIs < 3 ms: '  num2str(100 * Output.Single_rating) '%']};
    
    figure; % raster
    set(gcf, 'Position', [50    40   450*N_conditions   956])
   
    
    for c=1:N_conditions
        L=cfg.condition(c).name;
        subplot(3,N_conditions,c)
        hold on
        set(gca, 'XTick', BINS(1:5:end), 'XTickLabel', BINS(1:5:end))
        xlim([cfg.analyse_states{3} cfg.analyse_states{4}]*1000)
        ylim([0 100])
        xlabel('Time from R-peak, s')
        ylabel('Number of Plotted Row')
        box on
        col=condition_colors{c};
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
            line([0 0],ylim,'color','k');
            title({['Raster for ' L ': 10-ms bins, '], ['100 R-peak intervals out of ' num2str(size(Output.(L).raster,1))], '(linearly spaced)'})
        end
        
    end
    
    for c=1:N_conditions
        L=cfg.condition(c).name;
        subplot(3,N_conditions,c+N_conditions)
        hold on
        box on
        xlim([cfg.analyse_states{3} cfg.analyse_states{4}]*1000)
        set(gca, 'XTick', BINS(1:5:end), 'XTickLabel', BINS(1:5:end))
        col=condition_colors{c};
        lineProps={'color',col,'linewidth',1};
        shadedErrorBar(BINS, Output.(L).SD, Output.(L).SD_SEM,lineProps,1);
        lineProps={'color',col,'linewidth',1,'linestyle',':'};
        shadedErrorBar(BINS, Output.(L).SDP, [Output.(L).SDPCu; Output.(L).SDPCL], lineProps,1);
        line([0 0],ylim,'color','k');
        ypos=NaN;
        if Output.(L).sig_sign==-1;
            ypos=min(Output.(L).SD)*-1;
        elseif Output.(L).sig_sign==1;
            ypos=max(Output.(L).SD);
        end
        to_plot=Output.(L).sig;
        to_plot(to_plot==0)=NaN;
        plot(BINS,to_plot*ypos,'color',col,'linewidth',5);
        
        y_lims=get(gca,'ylim');
        text(BINS(10),y_lims(2)-diff(y_lims)*c*1/20, [L ': trials = ' ,num2str(Output.(L).NrTrials), ' events = ' ,num2str(Output.(L).NrEvents) ],'Color',condition_colors{c}, 'Interpreter', 'none');
    end
    
    ylabel('Firing rate');
    xlabel('time to Rpeak');
    
    for c=1:N_conditions
        L=cfg.condition(c).name;
        col=condition_colors{c};
        subplot(3,N_conditions,c+2*N_conditions);
        hold on
        box on
        % yyaxis left % dont even know what this does - but it doesnt work
        % in matlab 2015
        hstR=Output.(L).Rds;
        hstP=Output.(L).Rds_perm;
        hstR=hstR/sum(hstR);
        hstP=hstP/sum(hstP);
        if ~isnan(hstR)
            stairs(cfg.spk.histbins, hstR, 'Color', col, 'LineWidth', 2)
        end
        if ~isnan(hstP)
            stairs(cfg.spk.histbins, hstP, 'Color', [0.5,0.5,0.5], 'LineWidth', 2)
        end
        ylabel('fraction of intervals');
        xlabel('RR Duration, s');
        title('RR durations');
        legend({L, 'jittered'})
    end
    filename= ['Raster_PSTH_Rintervals_' unit_ID, '__' target];
     if v < 2016
        sgtitleText = [Output.unit_ID '_' Output.target,'SNR: ' num2str(Output.quantSNR) '; Fano Factor: ' num2str(Output.stability_rating) '; % of ISIs < 3 ms: '  num2str(100 * Output.Single_rating) '%'];
        mtit(sgtitleText,'xoff', 0, 'yoff', 0.2,'interpreter','none'); %% feel free to make this work with a cell, but suptitle is not a thing in matlab 2015
    else
        sgtitle(sgtitleText,'interpreter','none')
    end
    
    export_fig([basepath_to_save, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
end
end