
%% set up parameters
N_monkeys        = 2;
N_conditions     = 2;
N_areas          = 3;
subset_names     = {'stable', 'selected'};
monkey_names     = {'Bacchus','Magnus'};
condition_names  = {'Rest', 'Task'};
condition_colors = {[0 0 1], [1 0 0]};
area_list        = {'VPL', 'dPul', 'MD'};

output_folder='Y:\Manuscripts\2024_Thalamus_ephys_heart_brain\';

for s = 1:2
    
    figure('Name',sprintf(['BarPlot_Pc_' subset_names{s}]),'Position',[728 347 600 466],'PaperPositionMode', 'auto');
    
    for m = 1:2 % loop through monkeys
        
       % load data
        load(['Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_' monkey_names{m} ...
            '_TaskRest_state4\Population_time_domain_per_unit_-0.25-0.25s__after_SNR_exclusion_' subset_names{s} '_noLow_amplitude_ccs_any_VPL_dPul_MD\Output.mat'])
        
%             load(['Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_' monkey_names{m} ...
%                 '_TaskRest\Population_time_domain_after_SNR_exclusion_selected_noLow_amplitude_ccs_any_VPL_dPul_MD\Output.mat'])
        
        for c=1:N_conditions
            L=condition_names{c};
            
            % prepare color matrix
            bar_colors_merged = [condition_colors{c}; 0.7 0.7 0.7];
            
            % prepare data matrix for plotting
            pc_mat = [];
            nb_mat = [];
            for a = 1: N_areas
                T=area_list{a};
                pc_mat = vertcat(pc_mat, [Out.(T).(L).Pc_SignFR(1)+Out.(T).(L).Pc_SignFR(2) Out.(T).(L).Pc_SignFR(3)]);
                nb_mat = vertcat(nb_mat, [Out.(T).(L).Nb_SignFR(1)+Out.(T).(L).Nb_SignFR(2) Out.(T).(L).Nb_SignFR(3)]);
            end
            subplot(2,N_conditions, 2*(m-1) + c);
            b = bar(pc_mat,'stacked');
            for ii = 1:length(b)
                b(ii).FaceColor = bar_colors_merged(ii,:);
            end
            
            
            
            % add percentages to bars
            xbarCnt = [b(1).XData; b(2).XData];
            ybarTop = [b(1).YData; b(1).YData+b(2).YData];
            ybarCnt = ybarTop - pc_mat'/2;
            
            pc_mat=round(pc_mat'*10)/10;
            nb_mat=nb_mat';
            
            % Create text strings
            %txt = compose('%.1f%%',pc_mat');
            for k=1:numel(pc_mat)
            th = text(xbarCnt(k), ybarCnt(k), {num2str(nb_mat(k)); [sprintf('%.1f%', pc_mat(k)) '%']}, ...
                'HorizontalAlignment', 'center', ....
                'VerticalAlignment', 'middle', ...
                'Color', 'w',....
                'FontSize', 8, ...
                'FontWeight', 'bold');
            end
            title(['Monkey ' monkey_names{m}(1) ': ' L])
            % create x tick labels
            %         for ii = 1:3
            %             xtick_labels{ii} = [area_list{ii} ' (' num2str(length(dt.(area_list{ii}).unitId)) ')'];
            %         end
            %         set(gca,'XTickLabel',xtick_labels)
            
            %         title(L,'interpreter','none');
            set(gca,'XTickLabel',area_list,'fontsize',10);
            ylim([0 100])
            if c == 1
                ylabel('Percentage of Units')
            end
            % legend(b, {'increase FR', 'decrease FR', 'non-significant'}, 'Location', 'Best')
            ax = gca;
            ax.FontSize = 12;
            
            
            %         Tab = table(area_list, nb_mat(:,1), nb_mat(:,2), 'VariableNames', {'Brain Area', 'Resp', 'No Resp'});
            %         savename = [output_folder filesep 'UnitCounts_' L '.xlsx'];
            %         writetable(Tab, savename)
            %         clear Tab
            
        end
    end
    
    export_fig([output_folder 'Fig3_' ,subset_names{s}], '-pdf'); %,'-transparent'
    
        %     save_figure_as('Pc_CardiacRelatedUnits_Merged',output_folder,savePlot)

end

