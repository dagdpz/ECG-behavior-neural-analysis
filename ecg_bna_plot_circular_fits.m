function ecg_bna_plot_circular_fits(cfg)

basepath_to_save = [cfg.SPK_root_results_fldr filesep 'Circular_population_results'];
if ~isfolder(basepath_to_save)
    mkdir(basepath_to_save)
end

%% temporary - do plotting of fit parameters
phase_by_area_cond1_vMpos = [];
phase_by_area_cond2_vMpos = [];
area_id_vMpos             = [];

phase_by_area_cond1_vMneg = [];
phase_by_area_cond2_vMneg = [];
area_id_vMneg             = [];

C = colororder;
C(1,:) = [0 0 0];

var_list_ephys = {'NrEvents'};

var_list_cardioballistic = ...
    {'linear.yfit', 'linear.coefs', 'linear.rsquared', 'linear.pvalue', ...
    'cosine.startPoint', 'cosine.yfit', 'cosine.coefs', 'cosine.rsquared', 'cosine.pvalue', ...
    'vonMisesPos.startPoint', 'vonMisesPos.yfit', 'vonMisesPos.coefs', 'vonMisesPos.rsquared', 'vonMisesPos.pvalue', ...
    'vonMisesNeg.startPoint', 'vonMisesNeg.yfit', 'vonMisesNeg.coefs', 'vonMisesNeg.rsquared', 'vonMisesNeg.pvalue', ...
    ...
    'lowIBI_linear.yfit', 'lowIBI_linear.coefs', 'lowIBI_linear.rsquared', 'lowIBI_linear.pvalue', ...
    'highIBI_linear.yfit', 'highIBI_linear.coefs', 'highIBI_linear.rsquared', 'highIBI_linear.pvalue', ...
    ...
    'lowIBI_cosine.startPoint', 'lowIBI_cosine.yfit', 'lowIBI_cosine.coefs', 'lowIBI_cosine.rsquared', 'lowIBI_cosine.pvalue', ...
    'highIBI_cosine.startPoint', 'highIBI_cosine.yfit', 'highIBI_cosine.coefs', 'highIBI_cosine.rsquared', 'highIBI_cosine.pvalue', ...
    ...
    'lowIBI_vonMisesPos.startPoint', 'lowIBI_vonMisesPos.yfit', 'lowIBI_vonMisesPos.coefs', 'lowIBI_vonMisesPos.rsquared', 'lowIBI_vonMisesPos.pvalue', ...
    'highIBI_vonMisesPos.startPoint', 'highIBI_vonMisesPos.yfit', 'highIBI_vonMisesPos.coefs', 'highIBI_vonMisesPos.rsquared', 'highIBI_vonMisesPos.pvalue', ...
    ...
    'lowIBI_vonMisesNeg.startPoint', 'lowIBI_vonMisesNeg.yfit', 'lowIBI_vonMisesNeg.coefs', 'lowIBI_vonMisesNeg.rsquared', 'lowIBI_vonMisesNeg.pvalue', ...
    'highIBI_vonMisesNeg.startPoint', 'highIBI_vonMisesNeg.yfit', 'highIBI_vonMisesNeg.coefs', 'highIBI_vonMisesNeg.rsquared', 'highIBI_vonMisesNeg.pvalue'}; % ,

load([cfg.SPK_root_results_fldr filesep 'unit_lists\unitInfo_after_exclusion_stableTaskAndRest_noCB_corr.mat'], 'unit_ids', 'targets')

%     unqTargets = unique(targets);
unqTargets = {'VPL', 'dPul', 'MD'};

for groupNum = 1:length(cfg.spk.compare_conditions)
    cond1_num = cfg.spk.compare_conditions{groupNum}(1);
    cond2_num = cfg.spk.compare_conditions{groupNum}(2);
    
    for targNum = 1:length(unqTargets)
        
        currTargIds = cellfun(@(x) strcmp(x, unqTargets{targNum}), targets);
        curr_unit_ids = unit_ids(currTargIds);
        dt = ecg_bna_load_variables(cfg, curr_unit_ids, 'cardioballistic', 'data', var_list_cardioballistic);
        
        ephys_dt = ecg_bna_load_variables(cfg, curr_unit_ids, 'per_unit', 'Output', var_list_ephys);
        
        %% PARAMETERS BY CONDITION
        %% stability vs. linear gof
        plot_stability_vs_linear_gof(dt, {'linear'}, cfg, cond1_num, cond2_num, unqTargets{targNum}, 'Rsq_vs_Stability_', basepath_to_save)
        
        plot_stability_vs_linear_gof(dt, {'lowIBI_linear', 'highIBI_linear'}, cfg, cond1_num, cond2_num, unqTargets{targNum}, 'MedSplit_Rsq_vs_Stability_', basepath_to_save)
        
        %% R-interval numbers
        plot_RpeakNum_histograms(ephys_dt, '', cfg, cond1_num, cond2_num, unqTargets{targNum}, 'HeartCycleNumber_', basepath_to_save)
        
        %% Rsquared histograms
        plot_Rsquared_histograms(dt, {''}, cfg, cond1_num, cond2_num, C, unqTargets{targNum}, 'R-squared_', basepath_to_save)
        
        plot_Rsquared_histograms(dt, {'lowIBI_', 'highIBI_'}, cfg, cond1_num, cond2_num, C, unqTargets{targNum}, 'MedSplit_R-squared_', basepath_to_save)
        
        %% Phases of the peak from different fits
        plot_phase_histograms(dt, {''}, cfg, cond1_num, cond2_num, C, unqTargets{targNum}, 'PeakPhase_' , basepath_to_save)
        
        plot_phase_histograms(dt, {'lowIBI_', 'highIBI_'}, cfg, cond1_num, cond2_num, C, unqTargets{targNum}, 'MedSplit_PeakPhase_' , basepath_to_save)
        
        %% scaling factor for von Mises fits
        plot_scalingFactor_histograms(dt, {''}, cfg, cond1_num, cond2_num, C, unqTargets{targNum}, 'ScalingFactorDifference_', basepath_to_save)
        
        plot_scalingFactor_histograms(dt, {'lowIBI_', 'highIBI_'}, cfg, cond1_num, cond2_num, C, unqTargets{targNum}, 'MedSplit_ScalingFactorDifference_', basepath_to_save)
        
        %% kappa - concentration parameter for von Mises fits
        plot_kappa_histograms(dt, {''}, cfg, cond1_num, cond2_num, C, unqTargets{targNum}, 'Kappas_', basepath_to_save)
        
        plot_kappa_histograms(dt, {'lowIBI_', 'highIBI_'}, cfg, cond1_num, cond2_num, C, unqTargets{targNum}, 'MedSplit_Kappas_', basepath_to_save)
        
        %% kappas vs. phase for von Mises fits
        % take kappas
        vMPos_kappa_cond1 = log(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(:,3));
        vMNeg_kappa_cond1 = log(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(:,3));
        
        vMPos_kappa_cond2 = log(dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(:,3));
        vMNeg_kappa_cond2 = log(dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(:,3));
        
        % take phases
        vMPos_Phase_cond1 = dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(:,4);
        vMNeg_Phase_cond1 = dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(:,4);
        
        vMPos_Phase_cond2 = dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(:,4);
        vMNeg_Phase_cond2 = dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(:,4);
        
        figure,
        set(gcf, 'Position', [557   620   754   206])
        subplot(1,2,1)
        scatter(vMPos_kappa_cond1,vMPos_Phase_cond1, '.', 'MarkerEdgeColor', cfg.condition(1).color)
        hold on
        scatter(vMPos_kappa_cond2,vMPos_Phase_cond2, '.', 'MarkerEdgeColor', cfg.condition(2).color)
        box on
        ylim([0 2*pi])
        title('Pos. von Mises')
        xlabel('log(\kappa)')
        ylabel('Peak Phase [0-2\pi]')
        legend({cfg.condition.name}, 'Location', 'Best')
        
        subplot(1,2,2)
        scatter(vMNeg_kappa_cond1,vMNeg_Phase_cond1, '.', 'MarkerEdgeColor', cfg.condition(1).color)
        hold on
        scatter(vMNeg_kappa_cond2,vMNeg_Phase_cond2, '.', 'MarkerEdgeColor', cfg.condition(2).color)
        box on
        ylim([0 2*pi])
        title('Neg. von Mises')
        
        filename = ['Kappa_vs_Phase_' unqTargets{targNum} '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
        export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
        
        %% figure out the better fit (von Mises Pos. vs. von Mises Neg.)
        % 2. sum of pvalues is lower for the better fit
        % 3. at least one of the conditions (Rest or Task) is fitted
        % significantly (p < 0.01)
        
        % for only one condition - significantly fitted with a linear
        % function
        % 1) pos linear fits - cond1 - rest
        lin_sig_pos_cond1 = dt.(cfg.condition(cond1_num).name).linear.pvalue(:,2) < 0.01 & ... % significant fit
            dt.(cfg.condition(cond1_num).name).linear.coefs(:,2) > 0 & ...
            ( (dt.(cfg.condition(cond1_num).name).linear.rsquared > dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared) & ...
            (dt.(cfg.condition(cond1_num).name).linear.rsquared > dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared) );
        
        % 2) pos linear fits - cond2 - task
        lin_sig_pos_cond2 = dt.(cfg.condition(cond2_num).name).linear.pvalue(:,2) < 0.01 & ... % significant fit
            dt.(cfg.condition(cond2_num).name).linear.coefs(:,2) > 0 & ...
            ( (dt.(cfg.condition(cond2_num).name).linear.rsquared > dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared) & ...
            (dt.(cfg.condition(cond2_num).name).linear.rsquared > dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared) );
        
        % 3) neg linear fits - cond1 - rest
        lin_sig_neg_cond1 = dt.(cfg.condition(cond1_num).name).linear.pvalue(:,2) < 0.01 & ... % significant fit
            dt.(cfg.condition(cond1_num).name).linear.coefs(:,2) < 0 & ...
            ( (dt.(cfg.condition(cond1_num).name).linear.rsquared > dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared) & ...
            (dt.(cfg.condition(cond1_num).name).linear.rsquared > dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared) );
        
        % 4) neg linear fits - cond2 - task
        lin_sig_neg_cond2 = dt.(cfg.condition(cond2_num).name).linear.pvalue(:,2) < 0.01 & ... % significant fit
            dt.(cfg.condition(cond1_num).name).linear.coefs(:,2) < 0 & ...
            ( (dt.(cfg.condition(cond2_num).name).linear.rsquared > dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared) & ...
            (dt.(cfg.condition(cond2_num).name).linear.rsquared > dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared) );
        
        % 5) von Mises pos - cond1 - rest
        % check that these units are fitted better with a circular fit in
        % both condition as I'm going to segregate units with at least one 
        % linear fit
        vmpos_sig_cond1 = dt.(cfg.condition(cond1_num).name).vonMisesPos.pvalue < 0.01 & ...
            dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared > dt.(cfg.condition(cond1_num).name).linear.rsquared & ...
            dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared > dt.(cfg.condition(cond2_num).name).linear.rsquared & ...
            dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared > dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared;
        
        % 6) von Mises pos - cond2 - task
        % check that these units are fitted better with a circular fit in
        % both condition as I'm going to segregate units with at least one 
        % linear fit
        vmpos_sig_cond2 = dt.(cfg.condition(cond2_num).name).vonMisesPos.pvalue < 0.01 & ...
            dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared > dt.(cfg.condition(cond1_num).name).linear.rsquared & ...
            dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared > dt.(cfg.condition(cond2_num).name).linear.rsquared & ...
            dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared > dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared;
        
        % 7) von Mises neg - cond1 - rest
        % check that these units are fitted better with a circular fit in
        % both condition as I'm going to segregate units with at least one 
        % linear fit
        vmneg_sig_cond1 = dt.(cfg.condition(cond1_num).name).vonMisesNeg.pvalue < 0.01 & ...
            dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared > dt.(cfg.condition(cond1_num).name).linear.rsquared & ...
            dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared > dt.(cfg.condition(cond2_num).name).linear.rsquared & ...
            dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared > dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared;
        
        % 8) von Mises neg - cond2 - task
        % check that these units are fitted better with a circular fit in
        % both condition as I'm going to segregate units with at least one 
        % linear fit
        vmneg_sig_cond2 = dt.(cfg.condition(cond2_num).name).vonMisesNeg.pvalue < 0.01 & ...
            dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared > dt.(cfg.condition(cond1_num).name).linear.rsquared & ...
            dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared > dt.(cfg.condition(cond2_num).name).linear.rsquared & ...
            dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared > dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared;
        
        vm_pos_ids = (dt.(cfg.condition(cond1_num).name).vonMisesPos.pvalue + dt.(cfg.condition(cond2_num).name).vonMisesPos.pvalue) < ...
            (dt.(cfg.condition(cond1_num).name).vonMisesNeg.pvalue + dt.(cfg.condition(cond2_num).name).vonMisesNeg.pvalue) & ...
            (min([dt.(cfg.condition(cond1_num).name).vonMisesPos.pvalue dt.(cfg.condition(cond2_num).name).vonMisesPos.pvalue],[],2) < 0.01 | ...
            min([dt.(cfg.condition(cond1_num).name).vonMisesNeg.pvalue dt.(cfg.condition(cond2_num).name).vonMisesNeg.pvalue],[],2) < 0.01);
        vm_neg_ids = (dt.(cfg.condition(cond1_num).name).vonMisesPos.pvalue + dt.(cfg.condition(cond2_num).name).vonMisesPos.pvalue) > ...
            (dt.(cfg.condition(cond1_num).name).vonMisesNeg.pvalue + dt.(cfg.condition(cond2_num).name).vonMisesNeg.pvalue) & ...
            (min([dt.(cfg.condition(cond1_num).name).vonMisesPos.pvalue dt.(cfg.condition(cond2_num).name).vonMisesPos.pvalue],[],2) < 0.01 | ...
            min([dt.(cfg.condition(cond1_num).name).vonMisesNeg.pvalue dt.(cfg.condition(cond2_num).name).vonMisesNeg.pvalue],[],2) < 0.01);
        
%         % figure out significance of fits in rest or task
%         vmpos_sig_cond1  = vm_pos_ids & dt.(cfg.condition(cond1_num).name).vonMisesPos.pvalue < 0.01;
%         vmpos_sig_cond2  = vm_pos_ids & dt.(cfg.condition(cond2_num).name).vonMisesPos.pvalue < 0.01;
%         vmneg_sig_cond1  = vm_neg_ids & dt.(cfg.condition(cond1_num).name).vonMisesNeg.pvalue < 0.01;
%         vmneg_sig_cond2  = vm_neg_ids & dt.(cfg.condition(cond2_num).name).vonMisesNeg.pvalue < 0.01;

        % all nonsig units
        selection.all_nonsig_ids      = ...
            ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
            ~vmpos_sig_cond1 & ~vmpos_sig_cond2 & ~vmneg_sig_cond1 & ~vmneg_sig_cond2;
        
        % pos von Mises

        selection.vmpos_sig_cond1_ids = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
            vmpos_sig_cond1 & ~vmpos_sig_cond2;
        selection.vmpos_sig_cond2_ids = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
            ~vmpos_sig_cond1 & vmpos_sig_cond2;
        selection.vmneg_sig_cond1_ids = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
            vmneg_sig_cond1 & ~vmneg_sig_cond2;
        selection.vmneg_sig_cond2_ids = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
            ~vmneg_sig_cond1 & vmneg_sig_cond2;
        
        selection.vmpos_sig_both_ids  = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
            vmpos_sig_cond1 & vmpos_sig_cond2;
        selection.vmneg_sig_both_ids  = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
            vmneg_sig_cond1 & vmneg_sig_cond2;
        
        selection.vmpos_nonsig_ids    = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
            ~vmpos_sig_cond1 & ~vmpos_sig_cond2;
        selection.vmneg_nonsig_ids    = ~(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) & ...
            ~vmneg_sig_cond1 & ~vmneg_sig_cond2;
        
        % R-squared for significant units
        bin_edges       = 10 .^ (-6:0.2:0);
        % R-squared
        lin_Rsq_cond1   = dt.(cfg.condition(cond1_num).name).linear.rsquared(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2);
        lin_Rsq_cond2   = dt.(cfg.condition(cond1_num).name).linear.rsquared(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2);
        
        vMPos_Rsq_cond1 = dt.(cfg.condition(cond1_num).name).vonMisesPos.rsquared(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids);
        vMNeg_Rsq_cond1 = dt.(cfg.condition(cond1_num).name).vonMisesNeg.rsquared(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids);
        
        vMPos_Rsq_cond2 = dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids);
        vMNeg_Rsq_cond2 = dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids);
        
        % bins
        lin_RsqBins_cond1   = histc(lin_Rsq_cond1, bin_edges);
        lin_RsqBins_cond2   = histc(lin_Rsq_cond2, bin_edges);
        
        vMPos_RsqBins_cond1 = histc(vMPos_Rsq_cond1, bin_edges);
        vMNeg_RsqBins_cond1 = histc(vMNeg_Rsq_cond1, bin_edges);
        
        vMPos_RsqBins_cond2 = histc(vMPos_Rsq_cond2, bin_edges);
        vMNeg_RsqBins_cond2 = histc(vMNeg_Rsq_cond2, bin_edges);
        
        % medians
        med_rest = [median(vMPos_Rsq_cond1) median(vMNeg_Rsq_cond1)];
        med_task = [median(dt.(cfg.condition(cond2_num).name).vonMisesPos.rsquared(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids)) ...
            median(dt.(cfg.condition(cond2_num).name).vonMisesNeg.rsquared(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids))];
        
        vMPos_Rsq_cond1     = vMPos_Rsq_cond1(:);
        vMNeg_RsqBins_cond1 = vMNeg_RsqBins_cond1(:);
        vMNeg_RsqBins_cond2 = vMNeg_RsqBins_cond2(:);
        lin_RsqBins_cond1   = lin_RsqBins_cond1(:);
        lin_RsqBins_cond2   = lin_RsqBins_cond2(:);
        
        figure,
        set(gcf, 'Position', [557   620   754   206])
        colororder(C(2:end,:))
        s1 = subplot(1,2,1);
        stairs(bin_edges, [vMPos_RsqBins_cond1 vMNeg_RsqBins_cond1 lin_RsqBins_cond1], 'LineWidth', 2)
        hold on
        plot(med_rest(1), 0, '.', 'Color', C(2,:), 'MarkerSize', 15)
        plot(med_rest(2), 0, '.', 'Color', C(3,:), 'MarkerSize', 15)
        set(gca,'XScale','log')
        xlabel('log10(R^2)')
        ylabel('Number of Units')
        title(cfg.condition(cond1_num).name)
        legend('von Mises Pos.', 'von Mises Neg.', 'Linear Fit', 'Location', 'Best')
        
        s2 = subplot(1,2,2);
        stairs(bin_edges, [vMPos_RsqBins_cond2 vMNeg_RsqBins_cond2 lin_RsqBins_cond2], 'LineWidth', 2)
        hold on
        plot(med_task(1), 0, '.', 'Color', C(2,:), 'MarkerSize', 15)
        plot(med_task(2), 0, '.', 'Color', C(3,:), 'MarkerSize', 15)
        set(gca,'XScale','log')
        xlabel('log10(R^2)')
        title(cfg.condition(cond2_num).name)
        
        linkaxes([s1 s2], 'xy')
        
        sgtitle([cfg.version ' ' unqTargets{targNum}], 'interpreter', 'none')
        
        filename = ['R-squaredSig_' unqTargets{targNum} '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
        export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
        
        p_wilcox(1) = signrank(vMPos_Rsq_cond1, vMPos_Rsq_cond2);
        
        if ~isempty(vMNeg_Rsq_cond1) & ~isempty(vMNeg_Rsq_cond2)
            p_wilcox(2) = signrank(vMNeg_Rsq_cond1, vMNeg_Rsq_cond2);
        else
            p_wilcox(2) = NaN;
        end
        
        T = table({'Pos. von Mises', 'Neg. von Mises'}', ...
            med_rest', med_task', p_wilcox', ...
            'VariableNames', {'Distribution', [cfg.condition(cond1_num).name ': median R^2'], [cfg.condition(cond2_num).name ': median R^2'], 'Paired Wilcoxon''s p'});
        writetable(T, [basepath_to_save filesep filename '.xls'])
        clear p_wilcox
        
        %% scaling factors of significant fits
        plot_scalingFactor_histograms(dt, {''}, cfg, cond1_num, cond2_num, C, unqTargets{targNum}, 'SigScalingFactorDifference_', basepath_to_save, selection)
        
        plot_scalingFactor_histograms(dt, {'lowIBI_', 'highIBI_'}, cfg, cond1_num, cond2_num, C, unqTargets{targNum}, 'MedSplit_SigScalingFactorDifference_', basepath_to_save, selection)
        
        % [Pos. von Mises] Phases of significant fits
        phase_bin_edges = [0:0.2:2*pi];
        
        cond1_pos.counts_sig_cond1       = histc(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(selection.vmpos_sig_cond1_ids,4), phase_bin_edges);
        cond1_pos.counts_sig_cond2       = histc(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(selection.vmpos_sig_cond2_ids,4), phase_bin_edges);
        cond1_pos.counts_sig_cond1_cond2 = histc(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(selection.vmpos_sig_both_ids,4), phase_bin_edges);
        cond1_pos.counts_nonsig          = histc(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(selection.vmpos_nonsig_ids,4), phase_bin_edges);
        
        cond2_pos.counts_sig_cond1       = histc(dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(selection.vmpos_sig_cond1_ids,4), phase_bin_edges);
        cond2_pos.counts_sig_cond2       = histc(dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(selection.vmpos_sig_cond2_ids,4), phase_bin_edges);
        cond2_pos.counts_sig_cond1_cond2 = histc(dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(selection.vmpos_sig_both_ids,4), phase_bin_edges);
        cond2_pos.counts_nonsig          = histc(dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(selection.vmpos_nonsig_ids,4), phase_bin_edges);
        
        cond1_pos.counts_nonsig    = cond1_pos.counts_nonsig(:);
        cond1_pos.counts_sig_cond1 = cond1_pos.counts_sig_cond1(:);
        cond1_pos.counts_sig_cond2 = cond1_pos.counts_sig_cond2(:);
        cond1_pos.counts_sig_cond1_cond2 = cond1_pos.counts_sig_cond1_cond2(:);
        
        cond2_pos.counts_sig_cond2 = cond2_pos.counts_sig_cond2(:);
        cond2_pos.counts_nonsig    = cond2_pos.counts_nonsig(:);
        cond2_pos.counts_sig_cond2 = cond2_pos.counts_sig_cond2(:);
        cond2_pos.counts_sig_cond1 = cond2_pos.counts_sig_cond1(:);
        cond2_pos.counts_sig_cond1_cond2 = cond2_pos.counts_sig_cond1_cond2(:);
        
        figure,
        set(gcf, 'Position', [557   620   754   206])
        s1 = subplot(1,2,1);
        colororder([0 0 1; 1 0 0; 0.5 0 0.5; 1 1 1]) % blue, red, magenta, white
        bar([0:0.2:2*pi], [cond1_pos.counts_sig_cond1 cond1_pos.counts_sig_cond2 cond1_pos.counts_sig_cond1_cond2 cond1_pos.counts_nonsig], 'stacked')
        xlabel('Heart Cycle Phase [0-2\pi]')
        ylabel('# Units')
        legend({['Only ' cfg.condition(cond1_num).name], ['Only ' cfg.condition(cond2_num).name], ['Both ' cfg.condition(cond2_num).name ' and ' cfg.condition(cond1_num).name], 'None'}, 'Location', 'best')
        title(cfg.condition(cond1_num).name)
        
        s2 = subplot(1,2,2);
        colororder([0 0 1; 1 0 0; 0.5 0 0.5; 1 1 1]) % blue, red, magenta, white
        bar([0:0.2:2*pi], [cond2_pos.counts_sig_cond1 cond2_pos.counts_sig_cond2 cond2_pos.counts_sig_cond1_cond2 cond2_pos.counts_nonsig], 'stacked')
        xlabel('Heart Cycle Phase [0-2\pi]')
        title(cfg.condition(cond2_num).name)
        
        linkaxes([s1 s2], 'xy')
        
        sgtitle([unqTargets{targNum} ': Pos. von Mises'])
        
        filename = ['SigPhasesVMPos_' unqTargets{targNum} '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
        export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
        
        % [Neg. von Mises]
        cond1_neg.counts_sig_cond1       = histc(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(selection.vmneg_sig_cond1_ids,4), phase_bin_edges);
        cond1_neg.counts_sig_cond2       = histc(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(selection.vmneg_sig_cond2_ids,4), phase_bin_edges);
        cond1_neg.counts_sig_cond1_cond2 = histc(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(selection.vmneg_sig_both_ids,4), phase_bin_edges);
        cond1_neg.counts_nonsig          = histc(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(selection.vmneg_nonsig_ids,4), phase_bin_edges);
        
        cond2_neg.counts_sig_cond1       = histc(dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(selection.vmneg_sig_cond1_ids,4), phase_bin_edges);
        cond2_neg.counts_sig_cond2       = histc(dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(selection.vmneg_sig_cond2_ids,4), phase_bin_edges);
        cond2_neg.counts_sig_cond1_cond2 = histc(dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(selection.vmneg_sig_both_ids,4), phase_bin_edges);
        cond2_neg.counts_nonsig          = histc(dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(selection.vmneg_nonsig_ids,4), phase_bin_edges);
        
        cond1_neg.counts_nonsig          = cond1_neg.counts_nonsig(:);
        cond1_neg.counts_sig_cond1       = cond1_neg.counts_sig_cond1(:);
        cond1_neg.counts_sig_cond2       = cond1_neg.counts_sig_cond2(:);
        cond1_neg.counts_sig_cond1_cond2 = cond1_neg.counts_sig_cond1_cond2(:);
        cond2_neg.counts_sig_cond2       = cond2_neg.counts_sig_cond2(:);
        cond2_neg.counts_nonsig          = cond2_neg.counts_nonsig(:);
        cond2_neg.counts_sig_cond2       = cond2_neg.counts_sig_cond2(:);
        cond2_neg.counts_sig_cond1       = cond2_neg.counts_sig_cond1(:);
        cond2_neg.counts_sig_cond1_cond2 = cond2_neg.counts_sig_cond1_cond2(:);
        
        figure,
        set(gcf, 'Position', [557   620   754   206])
        s1 = subplot(1,2,1);
        colororder([0 0 1; 1 0 0; 0.5 0 0.5; 1 1 1]) % blue, red, magenta, white
        bar(phase_bin_edges, [cond1_neg.counts_sig_cond1 cond1_neg.counts_sig_cond2 cond1_neg.counts_sig_cond1_cond2 cond1_neg.counts_nonsig], 'stacked')
        xlabel('Heart Cycle Phase [0-2\pi]')
        ylabel('# Units')
        legend({['Only ' cfg.condition(cond1_num).name], ['Only ' cfg.condition(cond2_num).name], ['Both ' cfg.condition(cond2_num).name ' and ' cfg.condition(cond1_num).name], 'None'}, 'Location', 'best')
        title(cfg.condition(cond1_num).name)
        
        s2 = subplot(1,2,2);
        colororder([0 0 1; 1 0 0; 0.5 0 0.5; 1 1 1]) % blue, red, magenta, white
        bar(phase_bin_edges, [cond2_neg.counts_sig_cond1 cond2_neg.counts_sig_cond2 cond2_neg.counts_sig_cond1_cond2 cond2_neg.counts_nonsig], 'stacked')
        xlabel('Heart Cycle Phase [0-2\pi]')
        title(cfg.condition(cond2_num).name)
        
        linkaxes([s1 s2], 'xy')
        
        sgtitle([unqTargets{targNum} ': Neg. von Mises'])
        
        filename = ['SigPhasesVMNeg_' unqTargets{targNum} '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
        export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
        
        
        % plot number of heart-responsive units in the rest and in the task
        pos_num_sig_cond1(targNum)          = 100 * sum(selection.vmpos_sig_cond1_ids) / length(curr_unit_ids);
        pos_num_sig_cond2(targNum)          = 100 * sum(selection.vmpos_sig_cond2_ids) / length(curr_unit_ids);
        pos_num_sig_cond1_cond2(targNum)    = 100 * sum(selection.vmpos_sig_both_ids) / length(curr_unit_ids);
        
        neg_num_sig_cond1(targNum)          = 100 * sum(selection.vmneg_sig_cond1_ids) / length(curr_unit_ids);
        neg_num_sig_cond2(targNum)          = 100 * sum(selection.vmneg_sig_cond2_ids) / length(curr_unit_ids);
        neg_num_sig_cond1_cond2(targNum)    = 100 * sum(selection.vmneg_sig_both_ids) / length(curr_unit_ids);
        
        all_num_nonsig(targNum)             = 100 * sum(selection.all_nonsig_ids) / length(curr_unit_ids);
        lin_num_sig_cond1_or_cond2(targNum) = 100 * sum(lin_sig_pos_cond1 | lin_sig_pos_cond2 | lin_sig_neg_cond1 | lin_sig_neg_cond2) / length(curr_unit_ids);
        
        if targNum == length(unqTargets)
            figure,
            pos_colors = brighten([0 0 0.5; 0.5 0 0; 0.5 0 0.5], 0.5);
            neg_colors = brighten([0 0 0.5; 0.5 0 0; 0.5 0 0.5], -0.25);
            colororder([pos_colors; 1 1 1; 0 0 0; neg_colors]) % blue, red, magenta, white, black
            bar([pos_num_sig_cond1; pos_num_sig_cond2; pos_num_sig_cond1_cond2; ...
                all_num_nonsig; lin_num_sig_cond1_or_cond2; ...
                neg_num_sig_cond1; neg_num_sig_cond2; neg_num_sig_cond1_cond2]', 'stacked')
            set(gca, 'XTickLabel', unqTargets)
            legend({['Pos. von Mises: Only ' cfg.condition(cond1_num).name], ...
                ['Pos. von Mises: Only ' cfg.condition(cond2_num).name], ...
                ['Pos. von Mises: Both ' cfg.condition(cond2_num).name ' and ' cfg.condition(cond1_num).name], ...
                'None', ...
                'Linear Slopes', ...
                ['Neg. von Mises: Only ' cfg.condition(cond1_num).name], ...
                ['Neg. von Mises: Only ' cfg.condition(cond2_num).name], ...
                ['Neg. von Mises: Both ' cfg.condition(cond2_num).name ' and ' cfg.condition(cond1_num).name]
                }, 'Location', 'southoutside')
            ylabel('Percentage of Units')
            
            filename = ['HistogramPos_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
            export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
            close all;
            
        end
        
%         if targNum == length(unqTargets)
%             
%             figure,
%             colororder([0 0 1; 1 0 0; 0.5 0 0.5; 1 1 1; 0 0 0]) % blue, red, magenta, white, black
%             bar([neg_num_sig_cond1; neg_num_sig_cond2; neg_num_sig_cond1_cond2; neg_num_nonsig; lin_num_sig_cond1_or_cond2]', 'stacked')
%             set(gca, 'XTickLabel', unqTargets)
%             legend({['Only ' cfg.condition(cond1_num).name], ...
%                 ['Only ' cfg.condition(cond2_num).name], ...
%                 ['Both ' cfg.condition(cond2_num).name ' and ' cfg.condition(cond1_num).name], ...
%                 'None', ...
%                 'Linear Slope'}, 'Location', 'southoutside')
%             ylabel('Percentage of Units')
%             
%             filename = ['HistogramNeg_' unqTargets{targNum} '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
%             export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
%             close all;
%             
%         end
        
%         % phase for significant fits
%         cos_Phase_cond1   = histc(dt.(cfg.condition(cond1_num).name).cosine.coefs(cos_sig_cond1,2), [0:0.2:2*pi]);
%         vMPos_Phase_cond1 = histc(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(vmpos_sig_cond1,4), [0:0.2:2*pi]);
%         vMNeg_Phase_cond1 = histc(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(vmneg_sig_cond1,4), [0:0.2:2*pi]);
%         
%         cos_Phase_cond2   = histc(dt.(cfg.condition(cond2_num).name).cosine.coefs(cos_sig_cond2,2), [0:0.2:2*pi]);
%         vMPos_Phase_cond2 = histc(dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(vmpos_sig_cond2,4), [0:0.2:2*pi]);
%         vMNeg_Phase_cond2 = histc(dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(vmneg_sig_cond2,4), [0:0.2:2*pi]);
%         
%         if size(cos_Phase_cond1,1) < size(cos_Phase_cond1,2)
%             cos_Phase_cond1 = cos_Phase_cond1';
%         end
%         
%         figure
%         set(gcf, 'Position', [557   620   754   206])
%         colororder(C)
%         s1 = subplot(1,2,1);
%         stairs([0:0.2:2*pi], [cos_Phase_cond1 vMPos_Phase_cond1 vMNeg_Phase_cond1], 'LineWidth', 2)
%         xlabel('Heart-Cycle Phase of Peak')
%         title(cfg.condition(cond1_num).name)
%         legend('Cosine', 'von Mises Pos.', 'von Mises Neg.')
%         xlim([0 2*pi])
%         
%         s2 = subplot(1,2,2);
%         stairs([0:0.2:2*pi], [cos_Phase_cond2 vMPos_Phase_cond2 vMNeg_Phase_cond2], 'LineWidth', 2)
%         xlabel('Heart-Cycle Phase of Peak')
%         title(cfg.condition(cond2_num).name)
%         xlim([0 2*pi])
%         linkaxes([s1 s2], 'xy')
%         
%         sgtitle([cfg.version ' ' unqTargets{targNum}], 'interpreter', 'none')
        
        %% plot significant kappas
        plot_kappa_histograms(dt, {''}, cfg, cond1_num, cond2_num, C, unqTargets{targNum}, 'SigKappas_', basepath_to_save, selection)
        
        plot_kappa_histograms(dt, {'lowIBI_', 'highIBI_'}, cfg, cond1_num, cond2_num, C, unqTargets{targNum}, 'MedSplit_SigKappas_', basepath_to_save, selection)

        %% plot kappas vs phases for significant units
        
        % take kappas
        vMPos_sigKappa_cond1 = log(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,3));
        vMNeg_sigKappa_cond1 = log(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,3));
        
        vMPos_sigKappa_cond2 = log(dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,3));
        vMNeg_sigKappa_cond2 = log(dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,3));
        
        % take phases
        vMPos_sigPhase_cond1 = dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,4);
        vMNeg_sigPhase_cond1 = dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,4);
        
        vMPos_sigPhase_cond2 = dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,4);
        vMNeg_sigPhase_cond2 = dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,4);
        
        figure,
        set(gcf, 'Position', [557   620   754   206])
        subplot(1,2,1)
        scatter(vMPos_sigKappa_cond1,vMPos_sigPhase_cond1, '.', 'MarkerEdgeColor', cfg.condition(1).color)
        hold on
        scatter(vMPos_sigKappa_cond2,vMPos_sigPhase_cond2, '.', 'MarkerEdgeColor', cfg.condition(2).color)
        box on
        ylim([0 2*pi])
        title('Pos. von Mises')
        xlabel('log(\kappa)')
        ylabel('Peak Phase [0-2\pi]')
        legend({cfg.condition.name}, 'Location', 'Best')
        
        subplot(1,2,2)
        scatter(vMNeg_sigKappa_cond1,vMNeg_sigPhase_cond1, '.', 'MarkerEdgeColor', cfg.condition(1).color)
        hold on
        scatter(vMNeg_sigKappa_cond2,vMNeg_sigPhase_cond2, '.', 'MarkerEdgeColor', cfg.condition(2).color)
        box on
        ylim([0 2*pi])
        title('Neg. von Mises')
        
        filename = ['SigKappa_vs_SigPhase_' unqTargets{targNum} '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
        export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
        
        
        %% compare circular distribution parameters
        disp(['Compare Means: ' cfg.condition(cond1_num).name ' vs. ' cfg.condition(cond2_num).name ': Pos. von Mises'])
        circ_wwtest(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(vm_pos_ids,4), dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(vm_pos_ids,4))
        disp(['Compare Concentrations: ' cfg.condition(cond1_num).name ' vs. ' cfg.condition(cond2_num).name ': Pos. von Mises'])
        circ_ktest(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(vm_pos_ids,4)', dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(vm_pos_ids,4)')
        
        disp(['Compare Means: ' cfg.condition(cond1_num).name ' vs. ' cfg.condition(cond2_num).name ': Neg. von Mises'])
        circ_wwtest(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(vm_neg_ids,4), dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(vm_neg_ids,4))
        disp(['Compare Concentrations: ' cfg.condition(cond1_num).name ' vs. ' cfg.condition(cond2_num).name ': Neg. von Mises'])
        circ_ktest(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(vm_neg_ids,4), dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(vm_neg_ids,4))
        
        % compare mean phases by area
        phase_by_area_cond1_vMpos = [phase_by_area_cond1_vMpos; dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(vm_pos_ids,4)];
        phase_by_area_cond2_vMpos = [phase_by_area_cond2_vMpos; dt.(cfg.condition(cond2_num).name).vonMisesPos.coefs(vm_pos_ids,4)];
        
        phase_by_area_cond1_vMneg = [phase_by_area_cond1_vMneg; dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(vm_neg_ids,4)];
        phase_by_area_cond2_vMneg = [phase_by_area_cond2_vMneg; dt.(cfg.condition(cond2_num).name).vonMisesNeg.coefs(vm_neg_ids,4)];
        
        area_id_vMpos            = [area_id_vMpos targNum*ones(1, length(dt.(cfg.condition(cond1_num).name).vonMisesPos.coefs(vm_pos_ids,4)))];
        area_id_vMneg            = [area_id_vMneg targNum*ones(1, length(dt.(cfg.condition(cond1_num).name).vonMisesNeg.coefs(vm_neg_ids,4)))];
        
        clear cond1 cond2
    end
    
    %% 2-way ANOVA for pos. von Mises distributions
    condition_id_vMpos = [zeros(1, length(phase_by_area_cond1_vMpos)) ones(1, length(phase_by_area_cond2_vMpos))];
    [pval_vMpos, stats_vMpos] = circ_hktest([phase_by_area_cond1_vMpos; phase_by_area_cond2_vMpos], [area_id_vMpos area_id_vMpos], condition_id_vMpos, 1, {'Area: VPL, dPul, MD', 'Condition: Rest, Task', 'Area * Condition'});
    
    T = table(stats_vMpos);
    writetable(T, [basepath_to_save filesep 'PosVM_Harrison-Kanji.xls'])
    
    figure, boxplot([phase_by_area_cond1_vMpos; phase_by_area_cond2_vMpos],[area_id_vMpos area_id_vMpos])
    set(gca, 'XTickLabel', {'VPL', 'dPul', 'MD'})
    ylim([0 2*pi])
    xlabel('Area')
    ylabel('Heart-Cycle Phase [0-2\pi]')
    
    figure, boxplot([phase_by_area_cond1_vMpos; phase_by_area_cond2_vMpos],condition_id_vMpos)
    set(gca, 'XTickLabel', {'Rest', 'Task'})
    ylim([0 2*pi])
    xlabel('Condition')
    ylabel('Heart-Cycle Phase [0-2\pi]')
    
    %% 2-way ANOVA for neg. von Mises distributions
    condition_id_vMneg = [zeros(1, length(phase_by_area_cond1_vMneg)) ones(1, length(phase_by_area_cond2_vMneg))];
    [pval_vMneg, stats_vMneg] = circ_hktest([phase_by_area_cond1_vMneg; phase_by_area_cond2_vMneg], [area_id_vMneg area_id_vMneg], condition_id_vMneg, 1, {'Area: VPL, dPul, MD', 'Condition: Rest, Task', 'Area * Condition'});
    
    T = table(stats_vMneg);
    writetable(T, [basepath_to_save filesep 'NegVM_Harrison-Kanji.xls'])
    
    figure, boxplot([phase_by_area_cond1_vMneg; phase_by_area_cond2_vMneg],[area_id_vMneg area_id_vMneg])
    set(gca, 'XTickLabel', {'VPL', 'dPul', 'MD'})
    ylim([0 2*pi])
    xlabel('Area')
    ylabel('Heart-Cycle Phase [0-2\pi]')
    
    figure, boxplot([phase_by_area_cond1_vMneg; phase_by_area_cond2_vMneg],condition_id_vMneg)
    set(gca, 'XTickLabel', {'Rest', 'Task'})
    ylim([0 2*pi])
    xlabel('Condition')
    ylabel('Heart-Cycle Phase [0-2\pi]')
    
end
end

function plot_stability_vs_linear_gof(data, lin_field_name, cfg, cond1_num, cond2_num, targName, file_prefix, basepath_to_save)

stability = data.stability_rating;

figure,
set(gcf, 'Position', [557   620   754   200*length(lin_field_name)])

for lin_field_num = 1:length(lin_field_name)
    
    Rsq_cond1{lin_field_num} = data.(cfg.condition(cond1_num).name).(lin_field_name{lin_field_num}).rsquared;
    Rsq_cond2{lin_field_num} = data.(cfg.condition(cond2_num).name).(lin_field_name{lin_field_num}).rsquared;
    
    [cc_cond1_tmp, p_cond1_tmp] = corrcoef(Rsq_cond1{lin_field_num}, stability);
    cc_cond1{lin_field_num}     = cc_cond1_tmp(1,2);
    p_cond1{lin_field_num}      = p_cond1_tmp(1,2);
    [cc_cond2_tmp, p_cond2_tmp] = corrcoef(Rsq_cond2{lin_field_num}, stability);
    cc_cond2{lin_field_num}     = cc_cond2_tmp(1,2);
    p_cond2{lin_field_num}      = p_cond2_tmp(1,2);
    clear cc_cond1_tmp p_cond1_tmp cc_cond2_tmp p_cond2_tmp
    
    s1(lin_field_num) = subplot(length(lin_field_name),length(cfg.condition),1+2*(lin_field_num-1));
    scatter(log10(Rsq_cond1{lin_field_num}), stability, '.', 'MarkerEdgeColor', cfg.condition(cond1_num).color)
    xlabel('log10(R^2 of Linear Fit)')
    ylabel('Stability')
    title({[cfg.condition(cond1_num).name ': ' lin_field_name{lin_field_num}],
        ['cc = ' num2str(cc_cond1{lin_field_num}) '; p = ' num2str(p_cond1{lin_field_num})]}, 'interpreter', 'none')
    box on
    
    s2(lin_field_num) = subplot(length(lin_field_name),length(cfg.condition),2*lin_field_num);
    scatter(log10(Rsq_cond2{lin_field_num}), stability, '.', 'MarkerEdgeColor', cfg.condition(cond2_num).color)
    title({[cfg.condition(cond2_num).name ': ' lin_field_name{lin_field_num}], ...
        ['cc = ' num2str(cc_cond2{lin_field_num}) '; p = ' num2str(p_cond2{lin_field_num})]}, 'interpreter', 'none')
    box on
    
end

linkaxes([s1 s2], 'xy')
    
sgtitle([cfg.version ' ' targName], 'interpreter', 'none')

filename = [file_prefix targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf');
close(gcf);
    
T = table({[cfg.condition(cond1_num).name ': ' lin_field_name{1}], [cfg.condition(cond2_num).name]}', ...
    [cc_cond1; cc_cond2], [p_cond1; p_cond2], ...
    'VariableNames', {'Condition', 'Pearson''s Rho', 'Pearson''s p-value'});
writetable(T, [basepath_to_save filesep filename '.xls'])

end

function plot_RpeakNum_histograms(ephys_dt, var_prefix, cfg, cond1_num, cond2_num, targName, file_prefix, basepath_to_save)
figure,
set(gcf, 'Position', [557   620   754   206])

M1 = nanmedian(ephys_dt.(cfg.condition(cond1_num).name).NrEvents);
M2 = nanmedian(ephys_dt.(cfg.condition(cond2_num).name).NrEvents);

s1 = subplot(1,2,1);
hist(ephys_dt.(cfg.condition(cond1_num).name).NrEvents);
title([cfg.condition(cond1_num).name '-' var_prefix ': median = ' num2str(M1) ' heart cycles'])
xlabel('Number of Heart Cycles')
ylabel('Number of Units')

s2 = subplot(1,2,2);
hist(ephys_dt.(cfg.condition(cond2_num).name).NrEvents);
title([cfg.condition(cond2_num).name '-' var_prefix ': median = ' num2str(M2) ' heart cycles'])

linkaxes([s1 s2], 'xy')

filename = [file_prefix targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf');
close(gcf)

if ~isempty(rmmissing(ephys_dt.(cfg.condition(cond1_num).name).NrEvents)) && ~isempty(rmmissing(ephys_dt.(cfg.condition(cond2_num).name).NrEvents))
    p_wilcox = signrank(ephys_dt.(cfg.condition(cond1_num).name).NrEvents, ephys_dt.(cfg.condition(cond2_num).name).NrEvents);
else
    p_wilcox = [];
end

T = table(M1, M2, p_wilcox','VariableNames', {[cfg.condition(cond1_num).name ': median NrEvents'], [cfg.condition(cond2_num).name ': median NrEvents'], 'Paired Wilcoxon''s p'});
writetable(T, [basepath_to_save filesep filename '.xls'])

end

function plot_Rsquared_histograms(dt, var_prefix, cfg, cond1_num, cond2_num, C, targName, file_prefix, basepath_to_save)

bin_edges              = 10 .^ (-6:0.2:0);

figure,
set(gcf, 'Position', [557   620   754   200*length(var_prefix)])
colororder(C)

for lin_field_num = 1:length(var_prefix)
    
    cos_Rsq_cond1{lin_field_num}   = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'cosine']).rsquared;
    vMPos_Rsq_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).rsquared;
    vMNeg_Rsq_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).rsquared;
    
    cos_Rsq_cond2{lin_field_num}   = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'cosine']).rsquared;
    vMPos_Rsq_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).rsquared;
    vMNeg_Rsq_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).rsquared;
    
    cos_Rsq_cond1_binned{lin_field_num}   = histc(cos_Rsq_cond1{lin_field_num}, bin_edges);
    vMPos_Rsq_cond1_binned{lin_field_num} = histc(vMPos_Rsq_cond1{lin_field_num}, bin_edges);
    vMNeg_Rsq_cond1_binned{lin_field_num} = histc(vMNeg_Rsq_cond1{lin_field_num}, bin_edges);

    cos_Rsq_cond2_binned{lin_field_num}   = histc(cos_Rsq_cond2{lin_field_num}, bin_edges);
    vMPos_Rsq_cond2_binned{lin_field_num} = histc(vMPos_Rsq_cond2{lin_field_num}, bin_edges);
    vMNeg_Rsq_cond2_binned{lin_field_num} = histc(vMNeg_Rsq_cond2{lin_field_num}, bin_edges);
    
    med_rest(:, lin_field_num) = {median(cos_Rsq_cond1{lin_field_num}); median(vMPos_Rsq_cond1{lin_field_num}); median(vMNeg_Rsq_cond1{lin_field_num})};
    med_task(:, lin_field_num) = {median(cos_Rsq_cond2{lin_field_num}); median(vMPos_Rsq_cond2{lin_field_num}); median(vMNeg_Rsq_cond2{lin_field_num})};
    
    s1(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),1+2*(lin_field_num-1));
    stairs(bin_edges, [cos_Rsq_cond1_binned{lin_field_num} vMPos_Rsq_cond1_binned{lin_field_num} vMNeg_Rsq_cond1_binned{lin_field_num}], 'LineWidth', 2)
    hold on
    plot(med_rest{1,lin_field_num}, 0, '.', 'Color', C(1,:), 'MarkerSize', 15)
    plot(med_rest{2,lin_field_num}, 0, '.', 'Color', C(2,:), 'MarkerSize', 15)
    plot(med_rest{3,lin_field_num}, 0, '.', 'Color', C(3,:), 'MarkerSize', 15)
    set(gca,'XScale','log')
    xlabel('log10(R^2)')
    ylabel('Number of Units')
    title(cfg.condition(cond1_num).name)
    legend('Cosine', 'von Mises Pos.', 'von Mises Neg.', 'Location', 'Best')
    
    s2(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),2*lin_field_num);
    stairs(bin_edges, [cos_Rsq_cond2_binned{lin_field_num} vMPos_Rsq_cond2_binned{lin_field_num} vMNeg_Rsq_cond2_binned{lin_field_num}], 'LineWidth', 2)
    hold on
    plot(med_task{1,lin_field_num}, 0, '.', 'Color', C(1,:), 'MarkerSize', 15)
    plot(med_task{2,lin_field_num}, 0, '.', 'Color', C(2,:), 'MarkerSize', 15)
    plot(med_task{3,lin_field_num}, 0, '.', 'Color', C(3,:), 'MarkerSize', 15)
    set(gca,'XScale','log')
    xlabel('log10(R^2)')
    title(cfg.condition(cond2_num).name)
    
    p_wilcox{1,lin_field_num} = signrank(dt.(cfg.condition(cond1_num).name).([var_prefix{1} 'cosine']).rsquared, dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'cosine']).rsquared);
    p_wilcox{2,lin_field_num} = signrank(dt.(cfg.condition(cond1_num).name).([var_prefix{1} 'vonMisesPos']).rsquared, dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).rsquared);
    p_wilcox{3,lin_field_num} = signrank(dt.(cfg.condition(cond1_num).name).([var_prefix{1} 'vonMisesNeg']).rsquared, dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).rsquared);
    
end

linkaxes([s1 s2], 'xy')

sgtitle([cfg.version ' ' targName], 'interpreter', 'none')

filename = [file_prefix targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
close(gcf)

T = table({'Cosine', 'Pos. von Mises', 'Neg. von Mises'}', ...
    med_rest, med_task, p_wilcox, ...
    'VariableNames', {'Distribution', [cfg.condition(cond1_num).name ': median R^2'], [cfg.condition(cond2_num).name ': median R^2'], 'Paired Wilcoxon''s p'});
writetable(T, [basepath_to_save filesep filename '.xls'])

end

function plot_phase_histograms(dt, var_prefix, cfg, cond1_num, cond2_num, C, targName, file_prefix, basepath_to_save)

phase_bin_edges = [0:0.2:2*pi];

figure
set(gcf, 'Position', [557   620   754   200*length(var_prefix)])
colororder(C)

for lin_field_num = 1:length(var_prefix)

    cos_Phase_cond1{lin_field_num}   = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'cosine']).coefs(:,2);
    vMPos_Phase_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,4);
    vMNeg_Phase_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,4);

    cos_Phase_cond2{lin_field_num}   = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'cosine']).coefs(:,2);
    vMPos_Phase_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,4);
    vMNeg_Phase_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,4);

    cos_Phase_cond1_binned{lin_field_num}   = histc(cos_Phase_cond1{lin_field_num}, phase_bin_edges);
    vMPos_Phase_cond1_binned{lin_field_num} = histc(vMPos_Phase_cond1{lin_field_num}, phase_bin_edges);
    vMNeg_Phase_cond1_binned{lin_field_num} = histc(vMNeg_Phase_cond1{lin_field_num}, phase_bin_edges);

    cos_Phase_cond2_binned{lin_field_num}   = histc(cos_Phase_cond2{lin_field_num}, phase_bin_edges);
    vMPos_Phase_cond2_binned{lin_field_num} = histc(vMPos_Phase_cond2{lin_field_num}, phase_bin_edges);
    vMNeg_Phase_cond2_binned{lin_field_num} = histc(vMNeg_Phase_cond2{lin_field_num}, phase_bin_edges);

    circ_med_rest(:,lin_field_num) = ...
        {mod(circ_median(cos_Phase_cond1{lin_field_num}), 2*pi);
        mod(circ_median(vMPos_Phase_cond1{lin_field_num}), 2*pi);
        mod(circ_median(vMNeg_Phase_cond1{lin_field_num}), 2*pi)};
    circ_med_task(:,lin_field_num) = ...
        {mod(circ_median(cos_Phase_cond2{lin_field_num}), 2*pi);
        mod(circ_median(vMPos_Phase_cond2{lin_field_num}), 2*pi);
        mod(circ_median(vMNeg_Phase_cond2{lin_field_num}), 2*pi)};

    cos_Phase_cond1_binned{lin_field_num} = cos_Phase_cond1_binned{lin_field_num}(:);

    s1(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),1+2*(lin_field_num-1));
    stairs(phase_bin_edges, [cos_Phase_cond1_binned{lin_field_num} vMPos_Phase_cond1_binned{lin_field_num} vMNeg_Phase_cond1_binned{lin_field_num}], 'LineWidth', 2)
    hold on
    plot(circ_med_rest{1,lin_field_num}, 0, '.', 'Color', C(1,:), 'MarkerSize', 15)
    plot(circ_med_rest{2,lin_field_num}, 0, '.', 'Color', C(2,:), 'MarkerSize', 15)
    plot(circ_med_rest{3,lin_field_num}, 0, '.', 'Color', C(3,:), 'MarkerSize', 15)
    hold on

    xlabel('Heart-Cycle Phase of Peak')
    title(cfg.condition(cond1_num).name)
    legend('Cosine', 'von Mises Pos.', 'von Mises Neg.')
    xlim([0 2*pi])

    s2(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),2*lin_field_num);
    stairs(phase_bin_edges, [cos_Phase_cond2_binned{lin_field_num} vMPos_Phase_cond2_binned{lin_field_num} vMNeg_Phase_cond2_binned{lin_field_num}], 'LineWidth', 2)
    hold on
    plot(circ_med_task{1,lin_field_num}, 0, '.', 'Color', C(1,:), 'MarkerSize', 15)
    plot(circ_med_task{2,lin_field_num}, 0, '.', 'Color', C(2,:), 'MarkerSize', 15)
    plot(circ_med_task{3,lin_field_num}, 0, '.', 'Color', C(3,:), 'MarkerSize', 15)
    xlabel('Heart-Cycle Phase of Peak')
    title(cfg.condition(cond2_num).name)
    xlim([0 2*pi])

    p_wwtest{lin_field_num,1} = circ_wwtest(cos_Phase_cond1_binned{lin_field_num}, cos_Phase_cond2_binned{lin_field_num});
    p_wwtest{lin_field_num,2} = circ_wwtest(vMPos_Phase_cond1_binned{lin_field_num}, vMPos_Phase_cond2_binned{lin_field_num});
    p_wwtest{lin_field_num,3} = circ_wwtest(vMNeg_Phase_cond1_binned{lin_field_num}, vMNeg_Phase_cond2_binned{lin_field_num});
    
end
   
linkaxes([s1 s2], 'xy')

sgtitle([cfg.version ' ' targName], 'interpreter', 'none')

filename = [file_prefix targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
close all;

T = table({'Cosine', 'Pos. von Mises', 'Neg. von Mises'}', circ_med_rest, circ_med_task, p_wwtest', ...
    'VariableNames', {'Distribution', [cfg.condition(cond1_num).name ': median phase'], [cfg.condition(cond2_num).name ': median phase'], 'Watson-Williams''s p'});
writetable(T, [basepath_to_save filesep filename '.xls'])
end

function plot_scalingFactor_histograms(dt, var_prefix, cfg, cond1_num, cond2_num, C, targName, file_prefix, basepath_to_save, selection)

bin_centers_a    = -0.195:0.01:0.195;

bin_centers_diff = -0.145:0.01:0.145;

figure
set(gcf, 'Position', [557   620   754   200*length(var_prefix)])
colororder(C(2:end,:))

for lin_field_num = 1:length(var_prefix)
    
    if exist('selection', 'var')
    % select only significant units
        vMPos_a_cond1{lin_field_num} = ...
            dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,1);
        vMNeg_a_cond1{lin_field_num} = ...
            dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,1);
    
        vMPos_a_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,1);
        vMNeg_a_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,1);
    else
        vMPos_a_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,1);
        vMNeg_a_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,1);
    
        vMPos_a_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,1);
        vMNeg_a_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,1);
    end
    
    % bin scaling factors
    vMPos_a_cond1_binned{lin_field_num} = histc(vMPos_a_cond1{lin_field_num}, bin_centers_a);
    vMNeg_a_cond1_binned{lin_field_num} = histc(vMNeg_a_cond1{lin_field_num}, bin_centers_a);
    
    vMPos_a_cond2_binned{lin_field_num} = histc(vMPos_a_cond2{lin_field_num}, bin_centers_a);
    vMNeg_a_cond2_binned{lin_field_num} = histc(vMNeg_a_cond2{lin_field_num}, bin_centers_a);
    
    % plot scaling factor histograms
    s1(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),1+2*(lin_field_num-1));
    stairs(bin_centers_a, [vMPos_a_cond1_binned{lin_field_num} vMNeg_a_cond1_binned{lin_field_num}], 'LineWidth', 2)
    legend('von Mises Pos.', 'von Mises Neg.')
    title(cfg.condition(cond1_num).name)
    xlabel('Scaling Factor')
    ylabel('Unit Counts')
    
    s2(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),2*lin_field_num);
    stairs(bin_centers_a, [vMPos_a_cond2_binned{lin_field_num} vMNeg_a_cond2_binned{lin_field_num}], 'LineWidth', 2)
    title(cfg.condition(cond2_num).name)
    
    sgtitle([targName ': Scaling Factors'])
    
end

linkaxes([s1 s2], 'y')

filename = [file_prefix targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
close(gcf);
    
figure
set(gcf, 'Position', [557   620   754   200*length(var_prefix)])

for lin_field_num = 1:length(var_prefix)

    % compute differences of absolute scaling factors
    vMPos_a_absdiff = abs(vMPos_a_cond1{lin_field_num}) - abs(vMPos_a_cond2{lin_field_num});

    vMNeg_a_absdiff = abs(vMNeg_a_cond1{lin_field_num}) - abs(vMNeg_a_cond2{lin_field_num});
    
    % median difference in scaling factors
    M_pos{lin_field_num} = median(vMPos_a_absdiff);
    M_neg{lin_field_num} = median(vMNeg_a_absdiff);
    
    % binned scaling factor differences
    vMPos_a_absdiff_binned = histc(vMPos_a_absdiff, bin_centers_diff);
    vMNeg_a_absdiff_binned = histc(vMNeg_a_absdiff, bin_centers_diff);
    
    S1(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),1+2*(lin_field_num-1));
    stairs(bin_centers_diff, vMPos_a_absdiff_binned, 'LineWidth', 2)
    hold on
    plot(M_pos{lin_field_num}, 0, '.', 'MarkerSize', 15)
    title('Pos. von Mises: Scaling Factor Difference')
    xlabel('Scaling Factor (abs(Rest) - abs(Task))')
    
    S2(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),2*lin_field_num);
    stairs(bin_centers_diff, vMNeg_a_absdiff_binned, 'LineWidth', 2)
    hold on
    plot(M_neg{lin_field_num}, 0, '.', 'MarkerSize', 15)
    title('Neg. von Mises: Scaling Factor Difference')
    
    sgtitle([targName ': Scaling Factor Differences'])
    
    % check if means are zeros
    p_wilcoxon_paired{1,lin_field_num} = signrank(abs(vMPos_a_cond1{lin_field_num}),abs(vMPos_a_cond2{lin_field_num}));
    
    if ~isempty(vMNeg_a_cond1{lin_field_num}) && ~isempty(vMNeg_a_cond2{lin_field_num})
        p_wilcoxon_paired{2,lin_field_num} = signrank(abs(vMNeg_a_cond1{lin_field_num}),abs(vMNeg_a_cond2{lin_field_num}));
    else
        p_wilcoxon_paired{2,lin_field_num} = NaN;
    end

end

linkaxes([S1 S2], 'xy')

if strfind(file_prefix, 'MedSplit_')
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
close(gcf);

T = table(M_pos', M_neg', p_wilcoxon_paired', ...
    'VariableNames', {'vonMisesPos Scaling Factor Median Diff.', 'vonMisesNeg Scaling Factor Median Diff.', 'Paired Wilcoxon''s p'});
writetable(T, [basepath_to_save filesep filename '.xls'])
end

function plot_kappa_histograms(dt, var_prefix, cfg, cond1_num, cond2_num, C, targName, file_prefix, basepath_to_save, selection)

kappa_bin_edges = exp(-1:2);

figure
set(gcf, 'Position', [557   620   754   200*length(var_prefix)])
colororder(C(2:end,:))

for lin_field_num = 1:length(var_prefix)
    
    if exist('selection', 'var')
        % select only significant units
        vMPos_kappa_cond1{lin_field_num} = ...
            log(dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,3));
        vMNeg_kappa_cond1{lin_field_num} = ...
            log(dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,3));
        
        vMPos_kappa_cond2{lin_field_num} = ...
            log(dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,3));
        vMNeg_kappa_cond2{lin_field_num} = ...
            log(dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,3));
        
    else
        vMPos_kappa_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,3);
        vMNeg_kappa_cond1{lin_field_num} = dt.(cfg.condition(cond1_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,3);
        
        vMPos_kappa_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,3);
        vMNeg_kappa_cond2{lin_field_num} = dt.(cfg.condition(cond2_num).name).([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,3);
    end
    
    vMPos_kappa_cond1_binned{lin_field_num} = histc(vMPos_kappa_cond1{lin_field_num}, kappa_bin_edges);
    vMNeg_kappa_cond1_binned{lin_field_num} = histc(vMNeg_kappa_cond1{lin_field_num}, kappa_bin_edges);
    
    vMPos_kappa_cond2_binned{lin_field_num} = histc(vMPos_kappa_cond2{lin_field_num}, kappa_bin_edges);
    vMNeg_kappa_cond2_binned{lin_field_num} = histc(vMNeg_kappa_cond2{lin_field_num}, kappa_bin_edges);
    
    M_cond1(lin_field_num,:) = ...
        {median(vMPos_kappa_cond1{lin_field_num}); median(vMNeg_kappa_cond1{lin_field_num})};
    M_cond2(lin_field_num,:) = ...
        {median(vMPos_kappa_cond2{lin_field_num}); median(vMNeg_kappa_cond2{lin_field_num})};
    
    vMNeg_kappa_cond1_binned{lin_field_num} = ...
        vMNeg_kappa_cond1_binned{lin_field_num}(:);
    vMNeg_kappa_cond2_binned{lin_field_num} = ...
        vMNeg_kappa_cond2_binned{lin_field_num}(:);
    
    s1(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),1+2*(lin_field_num-1));
    bar(kappa_bin_edges,[vMPos_kappa_cond1_binned{lin_field_num} vMNeg_kappa_cond1_binned{lin_field_num}])
    xlabel('log(\kappa)')
    title(cfg.condition(cond1_num).name)
    legend({'von Mises Pos.', 'von Mises Neg.'})
    xlim([-2 4])
    
    s2(lin_field_num) = subplot(length(var_prefix),length(cfg.condition),2*lin_field_num);
    bar(kappa_bin_edges,[vMPos_kappa_cond2_binned{lin_field_num} vMNeg_kappa_cond2_binned{lin_field_num}])
    xlabel('log(\kappa)')
    title(cfg.condition(cond2_num).name)
    xlim([-2 4])
    
    p_wilcoxon_paired{lin_field_num,1} = signrank(vMPos_kappa_cond1{lin_field_num}, vMPos_kappa_cond2{lin_field_num});
    
    if ~isempty(vMNeg_kappa_cond1{lin_field_num}) && ~isempty(vMNeg_kappa_cond2{lin_field_num})
        p_wilcoxon_paired{2,lin_field_num} = signrank(abs(vMNeg_kappa_cond1{lin_field_num}),abs(vMNeg_kappa_cond2{lin_field_num}));
    else
        p_wilcoxon_paired{2,lin_field_num} = NaN;
    end
    
end

linkaxes([s1 s2], 'xy')

sgtitle([cfg.version ' ' targName], 'interpreter', 'none')

filename = [file_prefix targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
close all;

T = table({'Pos. von Mises', 'Neg. von Mises'}', M_cond1', M_cond2', p_wilcoxon_paired, ...
    'VariableNames', {'Distribution', [cfg.condition(cond1_num).name ': median \kappa'], [cfg.condition(cond2_num).name ': median \kappa'], 'Paired Wilcoxon''s p'});
writetable(T, [basepath_to_save filesep filename '.xls'])
end

function plot_kappas_vs_phase(data, var_prefix, cfg, cond1_num, cond2_num, C, targName, file_prefix, basepath_to_save, selection)

if exist('selection', 'var')
    % take kappas
    vMPos_kappa_cond1 = log(dt.(cfg.condition(cond1_num).name)([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,3));
    vMNeg_kappa_cond1 = log(dt.(cfg.condition(cond1_num).name)([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,3));
    
    vMPos_kappa_cond2 = log(dt.(cfg.condition(cond2_num).name)([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,3));
    vMNeg_kappa_cond2 = log(dt.(cfg.condition(cond2_num).name)([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,3));
    
    % take phases
    vMPos_Phase_cond1 = dt.(cfg.condition(cond1_num).name)([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,4);
    vMNeg_Phase_cond1 = dt.(cfg.condition(cond1_num).name)([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,4);
    
    vMPos_Phase_cond2 = dt.(cfg.condition(cond2_num).name)([var_prefix{lin_field_num} 'vonMisesPos']).coefs(selection.vmpos_sig_cond1_ids | selection.vmpos_sig_cond2_ids | selection.vmpos_sig_both_ids,4);
    vMNeg_Phase_cond2 = dt.(cfg.condition(cond2_num).name)([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(selection.vmneg_sig_cond1_ids | selection.vmneg_sig_cond2_ids | selection.vmneg_sig_both_ids,4);
else
    % take kappas
    vMPos_kappa_cond1 = log(data.(cfg.condition(cond1_num).name)([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,3));
    vMNeg_kappa_cond1 = log(data.(cfg.condition(cond1_num).name)([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,3));
    
    vMPos_kappa_cond2 = log(data.(cfg.condition(cond2_num).name)([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,3));
    vMNeg_kappa_cond2 = log(data.(cfg.condition(cond2_num).name)([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,3));
    
    % take phases
    vMPos_Phase_cond1 = data.(cfg.condition(cond1_num).name)([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,4);
    vMNeg_Phase_cond1 = data.(cfg.condition(cond1_num).name)([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,4);
    
    vMPos_Phase_cond2 = data.(cfg.condition(cond2_num).name)([var_prefix{lin_field_num} 'vonMisesPos']).coefs(:,4);
    vMNeg_Phase_cond2 = data.(cfg.condition(cond2_num).name)([var_prefix{lin_field_num} 'vonMisesNeg']).coefs(:,4);
end

figure,
set(gcf, 'Position', [557   620   754   206])
colororder(C(2:end,:))

subplot(1,2,1)
scatter(vMPos_kappa_cond1,vMPos_Phase_cond1, '.', 'MarkerEdgeColor', cfg.condition(1).color)
hold on
scatter(vMPos_kappa_cond2,vMPos_Phase_cond2, '.', 'MarkerEdgeColor', cfg.condition(2).color)
box on
ylim([0 2*pi])
title('Pos. von Mises')
xlabel('log(\kappa)')
ylabel('Peak Phase [0-2\pi]')
legend({cfg.condition.name}, 'Location', 'Best')

subplot(1,2,2)
scatter(vMNeg_kappa_cond1,vMNeg_Phase_cond1, '.', 'MarkerEdgeColor', cfg.condition(1).color)
hold on
scatter(vMNeg_kappa_cond2,vMNeg_Phase_cond2, '.', 'MarkerEdgeColor', cfg.condition(2).color)
box on
ylim([0 2*pi])
title('Neg. von Mises')

filename = [file_prefix targName '_' cfg.condition(cond1_num).name '_' cfg.condition(cond2_num).name];
export_fig(gcf, [basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
close(gcf)

end
