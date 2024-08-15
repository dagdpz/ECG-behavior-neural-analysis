function ecg_bna_population_correlation_analysis(cfg)

unqTargets = {'VPL', 'dPul', 'MD'};
N_Areas = length(unqTargets);
N_conditions = length(cfg.condition);

var_list = {'pearson_r', 'pearson_p', 'n_cycles', ...
    'FR_pearson_r', 'FR_pearson_p', 'FR_n_cycles', ...
    'RR_pearson_r', 'RR_pearson_p', 'RR_n_cycles'};

var_prefix = {'', 'FR_', 'RR_'};

bar_colors = [244 149 173; 126 221 95; 127 127 127]/255;

basepath_to_save = [cfg.SPK_root_results_fldr filesep 'Population_correlation_analysis'];
if ~exist(basepath_to_save, 'dir')
    mkdir(basepath_to_save)
end

%% load data
%load([cfg.SPK_root_results_fldr filesep 'unit_lists_ECG\unitInfo_after_SNR_exclusion_selected_noLow_amplitude_ccs_any.mat'], 'unit_ids', 'targets')
load([cfg.SPK_root_results_fldr filesep 'unit_lists_ECG\unitInfo_after_SNR_exclusion_stable_noLow_amplitude_ccs_any.mat'], 'unit_ids', 'targets', 'ids_both')

for a = 1: N_Areas
    T=unqTargets{a};
    currTargIds = cellfun(@(x) strcmp(x, unqTargets{a}), targets);
    curr_unit_ids = unit_ids(currTargIds);
    curr_ids_both = ids_both(currTargIds);
    dt.(T) = ecg_bna_load_variables(cfg, curr_unit_ids, 'correlation_analysis', 'data', var_list, curr_ids_both);
end

%% autocorrelation & correlation coefs between FR and RR
for prefixNum = 1:length(var_prefix)
    
    f0 = figure;
    set(f0, 'Position', [22 441 1857 491])
    
    for a = 1: N_Areas
        
        T=unqTargets{a};
        
        for c=1:N_conditions
            
            L=cfg.condition(c).name;
            
            splot_num = (c-1)*N_Areas + a;
            subplot(N_conditions, N_Areas, splot_num)
            plot(dt.(T).cc_lag_list', dt.(T).(L).([var_prefix{prefixNum} 'pearson_r']), 'Color', [0.5 0.5 0.5])
            
            xlim([-12 12])
            ylim([-1 1])
            
            title([T ':' L])
            
        end
        
    end
    
    subplot(N_conditions, N_Areas, 1)
    xlabel('Shift of RR sequence relative to FR, # cardiac cycles')
    ylabel('Correlation Coefficients')
    
    if isempty(var_prefix{prefixNum})
        sgttl = 'Correlation: FR vs. RR';
    elseif contains(var_prefix{prefixNum}, 'FR_')
        sgttl = 'Autocorrelation FR';
    elseif contains(var_prefix{prefixNum}, 'RR_')
        sgttl = 'Autocorrelation RR';
    end
    
    figure(f0)
    mtit(sgttl) %% sgtitle not working in matlab 2015
    save_figure_as(f0, [var_prefix{prefixNum} 'CC_Lags_AllAreas'], basepath_to_save, 1)
    
end

%% normalized correlation between FR and RR
f0 = figure;
set(f0, 'Position', [22 441 1857 491])

for a = 1: N_Areas
    
    T=unqTargets{a};
    
    for c=1:N_conditions
        
        L=cfg.condition(c).name;
        
        figure(f0);
        % compute normalized ccs
        R=dt.(T).(L).('pearson_r');
        P=dt.(T).(L).('pearson_p');
        
        norm_r  = R ./ repmat(R(ceil(end/2),:),size(R,1),1);
        pos_sig = R(ceil(end/2),:) > 0 & P(ceil(end/2),:) < 0.05;
        neg_sig = R(ceil(end/2),:) < 0 & P(ceil(end/2),:) < 0.05;
        non_sig = ~(pos_sig|neg_sig);
        
        splot_num = (c-1)*N_Areas + a;
        subplot(N_conditions, N_Areas, splot_num)
        hold on
        if sum(non_sig)>0
            plot(dt.(T).cc_lag_list, norm_r(:, non_sig), 'Color', bar_colors(3,:))
        end
        if sum(pos_sig)>0
            plot(dt.(T).cc_lag_list, norm_r(:, pos_sig), 'Color', bar_colors(1,:))
        end
        if sum(neg_sig)>0
            plot(dt.(T).cc_lag_list, norm_r(:, neg_sig), 'Color', bar_colors(2,:))
        end
        hline(0)
        hline(1)
        
        ylim([-2 2])
        title([T ': ' L])
        xlabel('Lag, # Cardiac Cycles')
        ylabel('cc normalized by cc(lag 0)')
        box on
    end
end

sgttl = 'Correlation: FR vs. RR';

figure(f0)
mtit(sgttl)
save_figure_as(f0, [var_prefix{prefixNum} 'CC_Lags_normalized_AllAreas'], basepath_to_save, 1)


%% lag scatters 
f1 = figure;
set(f1, 'Position', [22 441 1857 491])

for a = 1: N_Areas
    T=unqTargets{a};
    
    for c=1:N_conditions
        
        L=cfg.condition(c).name;
        
        % plot optimal correlation coefficients
        % extract data
        R=dt.(T).(L).('pearson_r');
        P=dt.(T).(L).('pearson_p');
        [~, ids]  = max(abs(R), [], 1); % find max abs ccs; I will need 'ids' to figure out significance of corresponding ccs
        ccs=R(sub2ind(size(R),ids,1:numel(ids)));
        invalid=isnan(ccs);
        P=P(:,~invalid);R=R(:,~invalid);ccs=ccs(~invalid);ids=ids(~invalid);
        P=P(ceil(end/2),:);
        R=R(ceil(end/2),:);
        sig=sign(R).*(P<0.05);
        
        splot_num = (c-1)*N_Areas + a;
        subplot(N_conditions, N_Areas, splot_num)
        hold on
        
        plot(dt.(T).cc_lag_list(ids(sig==0)),  ccs(sig==0),  'o', 'Color', bar_colors(3, :))
        plot(dt.(T).cc_lag_list(ids(sig==-1)), ccs(sig==-1), 'o', 'Color', bar_colors(2, :))
        plot(dt.(T).cc_lag_list(ids(sig==1)),  ccs(sig==1),  'o', 'Color', bar_colors(1, :))
        xlabel('Lag, # Heart Cycles')
        ylabel('Correlation Coefficient')
        title([T ': ' L])
        xlim([-13 13])
        ylim([-0.5 0.5])
        box on
    end
    
end

figure(f1)
save_figure_as(f1, 'CC_MaxAbsLags_AllAreas', basepath_to_save, 1)


%% histograms of lag of max cc
f10 = figure;
set(f10, 'Position', [22 441 1857 491])

for a = 1: N_Areas
    T=unqTargets{a};
    
    for c=1:N_conditions
        
        L=cfg.condition(c).name;
        
        % plot optimal correlation coefficients
        % extract data
        R=dt.(T).(L).('pearson_r');
        P=dt.(T).(L).('pearson_p');
        [~, ids]  = max(abs(R), [], 1); % find max abs ccs; I will need 'ids' to figure out significance of corresponding ccs
        ccs=R(sub2ind(size(R),ids,1:numel(ids)));
        invalid=isnan(ccs);
        P=P(:,~invalid);R=R(:,~invalid);ccs=ccs(~invalid);ids=ids(~invalid);
        P=P(ceil(end/2),:);
        R=R(ceil(end/2),:);
        sig=sign(R).*(P<0.05);
        
        splot_num = (c-1)*N_Areas + a;
        subplot(N_conditions, N_Areas, splot_num)
        hold on
        
        laglist=dt.(T).cc_lag_list(1,:);
        N=numel(ids);
        
        
        h0=hist(laglist(ids(sig==0)),laglist);
        h1=hist(laglist(ids(sig==-1)),laglist);
        h2=hist(laglist(ids(sig==1)),laglist);
        
        m0=median(laglist(ids(sig==0)));
        m1=median(laglist(ids(sig==-1)));
        m2=median(laglist(ids(sig==1)));
        
        
        plot(laglist,  h0, 'Color', bar_colors(3, :))
        plot(laglist,  h1, 'Color', bar_colors(2, :))
        plot(laglist,  h2, 'Color', bar_colors(1, :))
        plot([m0 m0],  [0 max(h0)], 'Color', bar_colors(3, :))
        plot([m1 m1],  [0 max(h1)], 'Color', bar_colors(2, :))
        plot([m2 m2],  [0 max(h2)], 'Color', bar_colors(1, :))
        
        
        xlabel('Lag, # Heart Cycles')
        ylabel('N units')
        title([T ': ' L])
        xlim([-13 13])
    end
    
end

figure(f10)
save_figure_as(f10, 'Histograms_MaxAbsLags_AllAreas', basepath_to_save, 1)


%% lag scatters task versus rest
for a = 1: N_Areas
    
    T=unqTargets{a};
    L1=cfg.condition(1).name;
    L2=cfg.condition(2).name;
    
    [rc, lag_rest] = max(abs(dt.(T).(L1).pearson_r));
    [tc, lag_task] = max(abs(dt.(T).(L2).pearson_r));
    invalid=isnan(rc)|isnan(tc);
    lag_rest(invalid)=[];
    lag_task(invalid)=[];
    
    f2 = figure;
    set(f2,'Position',[364   319   580   584])
    
    aa=...
        scatter(dt.(T).cc_lag_list(lag_rest),dt.(T).cc_lag_list(lag_task), 20, cfg.area_colors{a}, 'filled', 'o', 'MarkerFaceAlpha', 0.3);
    %         scatterhistogram(dt.(T).cc_lag_list(lag_rest),dt.(T).cc_lag_list(lag_task),...
    %             'Title',T,'HistogramDisplayStyle','smooth',...
    %             'ScatterPlotLocation','SouthEast','Color',cfg.area_colors{a}, 'MarkerAlpha',0.3)%,'Kernel','on','Marker','.'
    xlabel('Rest: Lag for Abs. Max. CC')
    ylabel('Task: Lag for Abs. Max. CC')
    title(T)
    box on
    save_figure_as(f2, ['Scatterhist_Lags_Rest_vs_Task_' T], basepath_to_save, 1)
    
end

%% CC histograms per area, for all lags and conditions
for a = 1: N_Areas
    
    figure,
    set(gcf, 'Position', [22 441 1857 491])
    colormap(bar_colors)
    bin_resolution=0.1;
    T=unqTargets{a};
    
    for c=1:N_conditions
        N_lags = size(dt.(T).cc_lag_list,2);
        L=cfg.condition(c).name;
        y_lim=1;
        
        for lag_num = 1:N_lags
            curr_cc = dt.(T).(L).pearson_r(lag_num, :);
            curr_pp = dt.(T).(L).pearson_p(lag_num, :);
            [~, sig_idx] = bonf_holm(curr_pp, 0.05); % correct for multiple comparisons
            pos_idx = curr_cc > 0;
            neg_idx = curr_cc < 0;
            
            % prepare data for plotting
            counts_pos    = histc(curr_cc(sig_idx & pos_idx)', -1:bin_resolution:1); % significant and positive
            counts_neg    = histc(curr_cc(sig_idx & neg_idx)', -1:bin_resolution:1); % significant and negative
            counts_nonsig = histc(curr_cc(~sig_idx)', -1:bin_resolution:1);
            counts_pos    = counts_pos(:);
            counts_neg    = counts_neg(:);
            counts_nonsig = counts_nonsig(:);
            
            % plot
            splot_num = (c-1)*N_lags + lag_num;
            plot_data = [counts_pos counts_neg counts_nonsig];
            subplot(N_conditions, N_lags, splot_num)
            b = bar(-1:bin_resolution:1, plot_data, 'stacked');
            for ii = 1:size(bar_colors,1)
                b(ii).FaceColor = bar_colors(ii,:);
            end
%             set(b, 'LineWidth',0.000000000000000000000000000000000001);
            xlim([-0.5 0.5])
            title(num2str(dt.(T).cc_lag_list(lag_num)))
            if lag_num > 1
                set(gca,'yticklabel',[]);
            else
                ylabel([T ': ' L ': N neurons']);
            end
            set(gca,'xtick',0,'xticklabel',0);
            xlabel('CC');
            if c == 1 && lag_num == 1
                legend({'Sig.Pos.', 'Sig.Neg.', 'Non-Sig.'})
            end
            
            y_lim=max(y_lim,max(sum(plot_data,2)));
        end
        
        for lag_num = 1:N_lags
            splot_num = (c-1)*N_lags + lag_num;
            subplot(N_conditions, N_lags, splot_num)
            set(gca,'ylim',[0,y_lim]);
        end
    end
    save_figure_as(gcf, ['Histograms_CC_FR_&_RR_' T], basepath_to_save, 1)
end

%% CC histograms for one lag, for areas and conditions
figure,
set(gcf, 'Position', [667 519 930 477])
lag_to_plot=0;
lag=cfg.correlation.lag_list==lag_to_plot;
for a = 1: N_Areas
    
    T=unqTargets{a};
    colormap(bar_colors)
    bin_resolution=0.05;
    
    y_lim=1;
    for c=1:N_conditions
        L=cfg.condition(c).name;
        
        curr_cc = dt.(T).(L).pearson_r(lag, :);
        curr_pp = dt.(T).(L).pearson_p(lag, :);
        [~, sig_idx] = bonf_holm(curr_pp, 0.05); % correct for multiple comparisons
        pos_idx = curr_cc > 0;
        neg_idx = curr_cc < 0;
        
        % prepare data for plotting
        counts_pos    = histc(curr_cc(sig_idx & pos_idx)', -1:bin_resolution:1); % significant and positive
        counts_neg    = histc(curr_cc(sig_idx & neg_idx)', -1:bin_resolution:1); % significant and negative
        counts_nonsig = histc(curr_cc(~sig_idx)', -1:bin_resolution:1);
        
        counts_pos = counts_pos(:);
        counts_neg = counts_neg(:);
        
        % plot
        splot_num = (c-1)*N_Areas + a;
        plot_data = [counts_pos counts_neg counts_nonsig];
        subplot(N_conditions, N_Areas, splot_num)
        b = bar(-1:bin_resolution:1, plot_data, 'stacked');
        set(b, 'FaceColor', 'Flat')
        for ii = 1:size(bar_colors,1)
            b(ii).FaceColor = bar_colors(ii,:);
        end
        xlim([-0.5 0.5])
        xlabel('CC between FR and RR duration')
        title([T ': ' L ', ' num2str(lag_to_plot) ' lag'])
        
        if a == 1 && c == 1
            legend({'Sig.Pos.', 'Sig.Neg.', 'Non-Sig.'})
        end
        
        y_lim=max(y_lim,max(sum(plot_data,2)));
        
    end
    
    for c=1:N_conditions
        splot_num = (c-1)*N_Areas + a;
        subplot(N_conditions, N_Areas, splot_num)
        set(gca,'ylim',[0,y_lim]);
    end
    
end
save_figure_as(gcf, ['Histograms_CC_FR_&_RR_at_lag_' num2str(lag_to_plot)], basepath_to_save, 1)

%% plot number of cardiac cycles used
figure,
set(gcf, 'Position', [667 519 930 477])
for a = 1: N_Areas
    T=unqTargets{a};
    for c=1:N_conditions
        L=cfg.condition(c).name;
        N_cycles=dt.(T).(L).n_cycles;        
        laglist=dt.(T).cc_lag_list;
        zero_lag_idx=find(laglist==0);
        N_cycles=N_cycles(zero_lag_idx,:);   
        N_cycles(isnan(N_cycles))=[];
        splot_num = (c-1)*N_Areas + a;
        subplot(N_conditions,N_Areas,splot_num)        
        hist(N_cycles, 20,'o-', 'Color', [0.5 0.5 0.5])
        title([T ': ' L ', min: ' num2str(min(N_cycles))])
    end
end
save_figure_as(gcf, ['Percentages_RR_Counts_' T], basepath_to_save, 1)


%% plot number of cardiac cycles used
figure,
set(gcf, 'Position', [667 519 930 477])
for a = 1: N_Areas
    T=unqTargets{a};
    for c=1:N_conditions
        L=cfg.condition(c).name;
        N_cycles=dt.(T).(L).n_cycles;
        N_cycles(:,all(isnan(N_cycles),1))=[];
        splot_num = (c-1)*N_Areas + a;
        subplot(N_conditions,N_Areas,splot_num)
        plot(cfg.correlation.lag_list, 100 * N_cycles ./ repmat(N_cycles(ceil(end/2),:),size(N_cycles,1),1), 'o-', 'Color', [0.5 0.5 0.5])
        title([T ': ' L])
    end
end
save_figure_as(gcf, ['Percentages_RR_Counts_' T], basepath_to_save, 1)

end

function save_figure_as(fig_id, filename,basepath_to_save,savePlot)
if savePlot;
    export_fig(fig_id,[basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
    close(fig_id)
end
end