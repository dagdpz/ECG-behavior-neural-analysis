clear all, close all

% This script creates the figure on correlation analysis summary for the
% manuscript Vasileva, Kaduk et al., 2024

unqTargets    = {'VPL', 'dPul', 'MD'};
unqConditions = {'Rest', 'Task'};

bin_resolution = 0.05;
bar_colors = [244 149 173; ... % pink  - pos
    126 221 95; ...    % green - neg
    127 127 127]/255;  % grey  - non-significant
deep_green = [6,64,43]/255;

%% load data
load('Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_Bacchus_TaskRest\correlation_analysis\Bac_20210716_11_VPL_R_correlation.mat')
load('Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_Magnus_TaskRest\Population_correlation_analysis_after_SNR_exclusion_stable_noLow_amplitude_ccs_any\data.mat')

figure,
set(gcf, 'Position', [160 237 1612 726])


%% sequences of FR and RRs in the rest
subplot(2,4,1)
yyaxis left
plot(data.Rest.timeRRstart, data.Rest.FRbyRR_Hz)
ylim([0 70])
yyaxis right
plot(data.Rest.timeRRstart, data.Rest.cycleDurations_s)
ylim([0.3 0.5])
xlim([3280 3290])
xlabel('Time from the recording onset, s')

%% histogram of ccs in example case at lag 0
subplot(2,4,2)
curr_cc = dt.VPL.Rest.pearson_r(ceil(end/2), :);
curr_pp = dt.VPL.Rest.pearson_p(ceil(end/2), :);
[sig_idx,~] = fdr_bky(curr_pp, 0.05);
pos_idx = curr_cc > 0;
neg_idx = curr_cc < 0;
nan_idx = isnan(curr_cc);

median_cc = median(curr_cc);

% prepare data for plotting
counts_pos    = histc(curr_cc(sig_idx & pos_idx),   -1+bin_resolution/2:bin_resolution:1-bin_resolution/2); % significant and positive
counts_neg    = histc(curr_cc(sig_idx & neg_idx),   -1+bin_resolution/2:bin_resolution:1-bin_resolution/2); % significant and negative
counts_nonsig = histc(curr_cc(~sig_idx & ~nan_idx), -1+bin_resolution/2:bin_resolution:1-bin_resolution/2);

counts_pos    = counts_pos(:);
counts_neg    = counts_neg(:);
counts_nonsig = counts_nonsig(:);

plot_data = [counts_pos counts_neg counts_nonsig];

hold on
b = bar([-1+bin_resolution/2:bin_resolution:1-bin_resolution/2]+bin_resolution/2, plot_data, 'stacked');
plot(median_cc, 0, '^', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none')
set(b, 'FaceColor', 'Flat')
for ii = 1:size(bar_colors,1)
    b(ii).FaceColor = bar_colors(ii,:);
end
xlim([-0.5 0.5])
set(gca,'XTick',-0.4:0.1:0.4,'XTickLabel',-0.4:0.1:0.4)
box on
legend({'Sig.Pos.', 'Sig.Neg.', 'Non-Sig.'})
xlabel('CC between FR and RR duration')
ylabel('Unit Counts')

%% histograms of unit FR-RR correlation signs
% loop through monkeys
unqMonkey = {'Bacchus', 'Magnus'};
sbpl_pos  = [3 7];
for curr_monk = 1:2
    load(['Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_' unqMonkey{curr_monk} '_TaskRest\Population_correlation_analysis_after_SNR_exclusion_stable_noLow_amplitude_ccs_any\data.mat'])
    for a = 1:3
        
        T = unqTargets{a};
        
        for c=1:2
            L = unqConditions{c};
            
            curr_cc = dt.(T).(L).pearson_r(ceil(end/2), :);
            curr_pp = dt.(T).(L).pearson_p(ceil(end/2), :);
            [~, sig_idx] = bonf_holm(curr_pp, 0.05); % correct for multiple comparisons
            [sig_idx,~] = fdr_bky(curr_pp, 0.05);
            pos_idx = curr_cc > 0;
            neg_idx = curr_cc < 0;
            nan_idx = isnan(curr_cc);
            
            counts_pos    = histc(curr_cc(sig_idx & pos_idx), -1:bin_resolution:1); % significant and positive
            counts_neg    = histc(curr_cc(sig_idx & neg_idx), -1:bin_resolution:1); % significant and negative
            counts_nonsig = histc(curr_cc(~sig_idx & ~nan_idx), -1:bin_resolution:1);
            
            counts_pos    = counts_pos(:);
            counts_neg    = counts_neg(:);
            counts_nonsig = counts_nonsig(:);
            
            unit_counts_by_area(c,a,:) = permute([sum(counts_pos) sum(counts_neg) sum(counts_nonsig)], [1 3 2]);
            
        end
    end
    
    for c = 1:2
        
        L = unqConditions{c};
        
        % proportion histograms
        subplot(2,4,sbpl_pos(curr_monk)-1+c)
        tmp = permute(unit_counts_by_area(c,:,:), [2 3 1]);
        unit_percentages = 100 * tmp ./ sum(tmp,2);
        b = bar(unit_percentages,'stacked');
        for ii = 1:size(bar_colors,1)
            b(ii).FaceColor = bar_colors(ii,:);
        end
        title(L)
        set(gca,'XTickLabel',unqTargets)
        ylabel('Unit fractions, %')
%         legend({'Sig.Pos.', 'Sig.Neg.', 'Non-Sig.'},'Location','southoutside')
        
    end
    
end

%% example unit with ccs for all the lags
[~,idx_max_rest] = max(abs(data.Rest.pearson_r));
[~,idx_max_task] = max(abs(data.Task.pearson_r));

subplot(2,4,5)
hold on
plot(data.cc_lag_list, data.Rest.pearson_r, 'bo-')
plot(data.cc_lag_list, data.Task.pearson_r, 'rs-')
plot(data.cc_lag_list(idx_max_rest), data.Rest.pearson_r(idx_max_rest), 'Marker', 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
plot(data.cc_lag_list(idx_max_task), data.Task.pearson_r(idx_max_task), 'Marker', 's', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])

box on
xlabel('Lag, # Cardiac Cycles')
ylabel('Correlation Coefficient')
legend({'Rest', 'Task'}, 'Location', 'best')

%% plot for lag of most effective correlation
load('Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_Magnus_TaskRest\Population_correlation_analysis_after_SNR_exclusion_stable_noLow_amplitude_ccs_any\data.mat')

laglist=dt.(T).cc_lag_list(1,:);
R=dt.dPul.Rest.('pearson_r');
P=dt.dPul.Rest.('pearson_p');
[~, ids]  = max(abs(R), [], 1); % find max abs ccs; I will need 'ids' to figure out significance of corresponding ccs
ccs=R(sub2ind(size(R),ids,1:numel(ids)));
invalid=isnan(ccs);
P=P(:,~invalid);
R=R(:,~invalid);
ccs=ccs(~invalid);
ids=ids(~invalid);
P=P(ceil(end/2),:);
R=R(ceil(end/2),:);
sig_idx = fdr_bky(P, 0.05); % figure out significant ones
sig=sign(R).*sig_idx;

% comptute histogram counts
h0=hist(laglist(ids(sig==0)),laglist);
h1=hist(laglist(ids(sig==-1)),laglist);
h2=hist(laglist(ids(sig==1)),laglist);

% compute medians by subgroup
m0(a,c)=median(laglist(ids(sig==0)));
m1(a,c)=median(laglist(ids(sig==-1)));
m2(a,c)=median(laglist(ids(sig==1)));
% compute lag medians only sig.
M(a,c) = median(laglist(ids(sig_idx)));

subplot(2,4,6)
x     = -20:0.125:20;
h1_h2_smo = smooth(interp1(laglist,h1+h2,x,'pchip'),50);
h1_smo    = smooth(interp1(laglist,h1,x,'pchip'),50);
h2_smo    = smooth(interp1(laglist,h2,x,'pchip'),50);
hold on

fill([x fliplr(x) x(1)], [h1_h2_smo; zeros(size(h1_h2_smo)); h1_h2_smo(1)]', deep_green, 'EdgeColor','none')
fill([x fliplr(x) x(1)], [h2_smo; zeros(size(h2_smo)); h2_smo(1)]', bar_colors(2,:), 'EdgeColor',deep_green,'FaceAlpha',1)
fill([x fliplr(x) x(1)], [h1_smo; zeros(size(h1_smo)); h1_smo(1)]', bar_colors(1,:), 'EdgeColor','none','FaceAlpha',0.7)
% plot(M(a,c),  0, '^', 'MarkerFaceColor', deep_green, 'MarkerEdgeColor', deep_green)         % median overall
plot(m1(a,c), 0, '^', 'MarkerFaceColor', bar_colors(2,:), 'MarkerEdgeColor', bar_colors(2,:)/2) % median sig pos
plot(m2(a,c), 0, '^', 'MarkerFaceColor', bar_colors(1,:), 'MarkerEdgeColor', bar_colors(1,:)/2) % median sig neg
box on

xlabel('Lag for Most Effective Correlation, # cardiac cycles')
ylabel('N units')
xlim([-20 20])
% ylim([-0.5 0.5])
box on
legend({'All Sig.', 'Sig. Neg.', 'Sig. Pos.', 'Median Neg.', 'Median Pos.'}, 'Location', 'best')
