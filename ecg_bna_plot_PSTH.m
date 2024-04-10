function ecg_bna_plot_PSTH(session_info, EPO,cfg)

% check the Matlab version
v = version('-release');
v = str2double(v(1:end-1));

N_conditions=numel(cfg.condition);
condition_colors={cfg.condition.color};
% find the current datafile
basepath_to_save=[cfg.SPK_root_results_fldr filesep 'per_unit'];
if ~exist(basepath_to_save,'dir')
    mkdir(basepath_to_save);
end

fileList = dir([basepath_to_save, filesep, session_info.Monkey(1:3) '_' session_info.Date '*.mat']);

xticks=EPO.xticks;
xtickslabels=EPO.xtickslabels;
histticks=EPO.histticks;
histtickslabels=EPO.histtickslabels;

eventnames=cfg.analyse_states(:,1);
for flNum = 1:length(fileList)
    load([basepath_to_save filesep fileList(flNum).name], 'Output')
    D=Output;
    
    unit_ID = D.unit_ID;
    target = D.target;
    
    sgtitleText = {[D.unit_ID '_' D.target], ['SNR: ' num2str(D.quantSNR) '; Fano Factor: ' num2str(D.stability_rating) '; % of ISIs < 3 ms: '  num2str(100 * D.Single_rating) '%']};
    
    figure;
    set(gcf, 'Position', [50    40   450*N_conditions   956])
    
    %% raster
    for c=1:N_conditions
        L=cfg.condition(c).name;
        col=condition_colors{c};
        subplot(3,N_conditions,c)
        hold on
    nr_events=repmat({0},numel(eventnames));
        for e=1:numel(eventnames)
            E=eventnames{e};
            if ~isfield(D.(L),E) || all(isnan((D.(L).(E).SD)))
                continue;
            end
            if D.(L).(E).NrEvents <= 100
                stepSize = 1:D.(L).(E).NrEvents; % just step through all the Rpeaks one by one
            else
                stepSize = round(linspace(1, D.(L).(E).NrEvents, 100));
            end
            BINS=EPO.tbinpos{e};
            
            x = D.(L).(E).raster(stepSize,:);
            xvals=repmat(BINS,size(x,1),1);
            xvals=xvals(x);
            yvals=repmat((1:size(x,1))',1,size(x,2),1);
            yvals=yvals(x);
            scatter(xvals,yvals,'.');
            
            y_lims=get(gca,'ylim');
            text(EPO.eventalignment(e), y_lims(2)-diff(y_lims)*1/2,E,'color','k','HorizontalAlignment','center','Interpreter','none' );
            nr_events{e}=[num2str(D.(L).(E).NrEvents) '|'];
        end
        title({[L ',' num2str(D.(L).NrTrials) 'Trials: 10-ms bins, '], ['max 100 out of ' nr_events{:} 'events']})
        
        line([EPO.eventalignment' EPO.eventalignment'],ylim,'color','k');
        %ylim([0 100])
        xlabel('Time to Event, ms')
        ylabel('Number of Plotted Row')
        set(gca, 'XTick', xticks, 'XTickLabel', xtickslabels)
        xlim([xticks(1) xticks(end)])
    end
    
    
    
    %% PSTH
    for c=1:N_conditions
        L=cfg.condition(c).name;
        subplot(3,N_conditions,c+N_conditions)
        hold on
        col=condition_colors{c};
        
        for e=1:numel(eventnames)
            E=eventnames{e};
            if ~isfield(D.(L),E) ||all(isnan((D.(L).(E).SD)))
                continue;
            else
                dat=D.(L).(E);
            end
            BINS=EPO.tbinpos{e};
            lineProps={'color',col,'linewidth',1};
            shadedErrorBar(BINS, dat.SD, dat.SD_SEM,lineProps,1);
            lineProps={'color',col,'linewidth',1,'linestyle',':'};
            shadedErrorBar(BINS, dat.SDP, [dat.SDPCU; dat.SDPCL], lineProps,1);
            
            ypos=NaN;
            if dat.sig_sign==-1;
                ypos=min(dat.SD)*-1;
            elseif dat.sig_sign==1;
                ypos=max(dat.SD);
            end
            to_plot=dat.sig;
            to_plot(to_plot==0)=NaN;
            plot(BINS,to_plot*ypos,'color',col,'linewidth',5);
            
            nr_events{e}=[num2str(D.(L).(E).NrEvents) '|'];
            
            y_lims=get(gca,'ylim');
            text(EPO.eventalignment(e), y_lims(2)-diff(y_lims)*1/2,E,'color','k','HorizontalAlignment','center','Interpreter','none' );
        end
        
        %title({[L ',' num2str(D.(L).NrTrials) 'Trials: 10-ms bins, '], ['max 100 out of ' nr_events{:} 'events']})
        line([EPO.eventalignment' EPO.eventalignment'],ylim,'color','k');
        ylabel('Firing rate');
        xlabel('time to Event, ms');
        set(gca, 'XTick', xticks, 'XTickLabel', xtickslabels)
        xlim([xticks(1) xticks(end)])
    end
    
    
    %% interval histograms
    for c=1:N_conditions
        L=cfg.condition(c).name;
        col=condition_colors{c};
        subplot(3,N_conditions,c+2*N_conditions);
        hold on
        for e=1:numel(eventnames)
            E=eventnames{e};
            if ~isfield(D.(L),E) || all(isnan((D.(L).(E).SD)))
                continue;
            else
                dat=D.(L).(E);
            end
            BINS=EPO.histbinpos{e};
            % yyaxis left % dont even know what this does - but it doesnt work
            % in matlab 2015
            hstR=dat.intervals;
            hstP=dat.shuffled_intervals;
            hstR=hstR/sum(hstR);
            hstP=hstP/sum(hstP);
            if ~isnan(hstR)
                stairs(BINS, hstR, 'Color', col, 'LineWidth', 2)
            end
            if ~isnan(hstP)
                stairs(BINS, hstP, 'Color', [0.5,0.5,0.5], 'LineWidth', 2)
            end
        end
        ylabel('fraction of intervals');
        xlabel('RR Duration, s');
        title('RR durations');
        legend({L, 'jittered'})
        set(gca, 'XTick', histticks, 'XTickLabel', histtickslabels)
        xlim([histticks(1) histticks(end)])
    end
    %text(histticks',zeros(size(histticks')), histtickslabels,'HorizontalAlignment',repmat({'right','left'},[1,numel(eventnames)]))
    
    %%
    filename= ['Raster_PSTH_Rintervals_' unit_ID, '__' target];
    if v < 2016
        sgtitleText = [D.unit_ID '_' D.target,'SNR: ' num2str(D.quantSNR) '; Fano Factor: ' num2str(D.stability_rating) '; % of ISIs < 3 ms: '  num2str(100 * D.Single_rating) '%'];
        mtit(sgtitleText,'xoff', 0, 'yoff', 0.2,'interpreter','none'); %% feel free to make this work with a cell, but suptitle is not a thing in matlab 2015
    else
        sgtitle(sgtitleText,'interpreter','none')
    end
    
    export_fig(gcf,[basepath_to_save, filesep, filename], '-pdf','-transparent') % pdf by run
    close(gcf);
end
end