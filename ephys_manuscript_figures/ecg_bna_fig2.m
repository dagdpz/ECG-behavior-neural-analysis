clear all, close all

% This script plots figure 2 for the Vasileva, Kaduk et al., 2024
% manuscript

%% folder
dir2save = 'Y:\Manuscripts\2024_Thalamus_ephys_heart_brain\';

%% set up parameters
N_conditions     = 2;
N_areas          = 3;
condition_names  = {'Rest', 'Task'};
condition_colors = {'b', 'r'};
area_list        = {'VPL', 'dPul', 'MD'};
BINS             = -250:5:250;

nRpeaks          = 400; % number of R-peaks for rasters
max_Rpeaks       = 5922;

%% load files
data.VPL  = load('Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_Bacchus_TaskRest\per_unit_-0.25-0.25s\Bac_20211222_07_VPL_R.mat','Output');
data.dPul = load('Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_Magnus_TaskRest\per_unit_-0.25-0.25s\Mag_20230524_29_dPul_L.mat','Output');
data.MD   = load('Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_Magnus_TaskRest\per_unit_-0.25-0.25s\Mag_20230621_38_MD_L.mat','Output');

figure;
set(gcf,'Position',[400 400 1164 568])

for a = 1:N_areas
    T = area_list{a};
    for c=1:N_conditions
        L=condition_names{c};
        col=condition_colors{c};
        
        %% plot raster with 1 ms resolution
        subplot(3,N_areas,3*(c-1)+a)
        % take 400 equally spaced R-peaks
        Rpeak_ids = round(1 : data.(T).Output.(L).NrEvents/nRpeaks : data.(T).Output.(L).NrEvents);
        [raster_row, raster_col] = find(data.(T).Output.(L).raster(Rpeak_ids,:)); % row - R-peak number, col - bin number
        scatter(raster_col-250,raster_row,3,'MarkerFaceColor',col,'MarkerFaceAlpha',0.4,'MarkerEdgeColor',col,'MarkerEdgeAlpha',0.4);
        hold on
        line([0 0],ylim,'color','k');
        set(gca, 'XTick', [-200 -100 0 100 200], 'XTickLabel', [-200 -100 0 100 200])
        set(gca, 'YTick', [])
        xlim([-250 250])
        
        box on
        if c == 1
            title(['Monkey ' data.(T).Output.unit_ID(1) '_' data.(T).Output.unit_ID(5:end) '_' T],'Interpreter','none')
            if a == 1
                xlabel('Time from R-peak, ms')
                ylabel('Number of R-peaks')
            end
        end
        
        %% plot R-peak triggered PETHs
        subplot(3,N_areas,6+a)
        hold on
        box on
        xlim([-0.25 0.25]*1000)
        lineProps={'color',col,'linewidth',1};
        s1 = shadedErrorBar(BINS, data.(T).Output.(L).SD, data.(T).Output.(L).SD_SEM,lineProps,1);
        lineProps={'color',col,'linewidth',1,'linestyle',':'};
        s2 = shadedErrorBar(BINS, data.(T).Output.(L).SDP, [data.(T).Output.(L).SDPCu; data.(T).Output.(L).SDPCL], lineProps,1);
        ypos = max(abs(data.(T).Output.(L).SD)) * data.(T).Output.(L).sig_sign;
        to_plot=data.(T).Output.(L).sig;
        to_plot(to_plot==0)=NaN;
        plot(BINS,to_plot*ypos,'color',col,'linewidth',5);
        
        if c == 2
            
            y_lims=get(gca,'ylim');
            
            % extend ylims by 1 on both sides
            y_lims(1) = y_lims(1)-1;
            y_lims(2) = y_lims(2)+1;
            set(gca,'YLim',y_lims)
            
            line([0 0],ylim,'color','k');
            
        end
        
        if a == 1
            
            ylabel('Firing Rate, Hz')
            
        elseif a == 3
            
            legend([s1.mainLine s2.mainLine s2.patch], ...
                {'Average PETH','Average Surrogate Data','95% Confidence Interval'}, ...
                'Location','Best')
            
        end
    end
    
    %% plot rasters with all collected R-peaks
    figure,
    currMaxRpeakNum = max([data.(T).Output.Rest.NrEvents data.(T).Output.Task.NrEvents]);
    set(gcf,'Position',[683    41   560   currMaxRpeakNum*955/max_Rpeaks])
    for c=1:N_conditions
        L=condition_names{c};
        col=condition_colors{c};
        
        subplot(1,2,c)
        [raster_row, raster_col] = find(data.(T).Output.(L).raster); % row - R-peak number, col - bin number
        scatter(raster_col-250,raster_row,3,'MarkerFaceColor',col,'MarkerFaceAlpha',0.4,'MarkerEdgeColor',col,'MarkerEdgeAlpha',0.4);
        hold on
        line([0 0],ylim,'color','k');
        set(gca, 'XTick', [-200 -100 0 100 200], 'XTickLabel', [-200 -100 0 100 200])
        set(gca, 'YTick', [0 data.(T).Output.(L).NrEvents])
        xlim([-250 250])
        ylim([0 data.(T).Output.(L).NrEvents])
        box on
        
        if c == 1
            
            xlabel('Time from R-peak, ms')
            ylabel('Number of Collected R-peaks')
        
        end
        
    end
    save_figure_as(['FigS2_' T],dir2save,1)
end
save_figure_as('Fig2_',dir2save,1)

function save_figure_as(filename,basepath_to_save,savePlot)
if savePlot
    export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
    close(gcf);
end
end
