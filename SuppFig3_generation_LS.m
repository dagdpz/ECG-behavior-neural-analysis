% to plot the whol sessions ECG bpm
function SuppFig3_generation_LS
clear all
load('Y:\Projects\Pulv_bodysignal\HeartRate_Changes_all_Blocks\Bpm_both_monkeys_taskrest.mat')
foi='mean_bpm'; %'median_bpm'; %'std_bpm'
cols=[1 0 0; 0 0 1];

Alltogether=[Task.Magnus.(foi), Rest.Magnus.(foi), Task.Bacchus.(foi), Rest.Bacchus.(foi)];
y_lims=[floor(min(Alltogether)) ceil(max(Alltogether))];


%% Violins - BPM
h=figure;
subplot(1,2,1)
% Combine data
Magnus_all_HR = [Task.Magnus.(foi), Rest.Magnus.(foi)];
% Create group labels
group = [ones(1, numel(Task.Magnus.(foi))), ones(1, numel(Rest.Magnus.(foi)))*2];
violinplot(Magnus_all_HR, group, 'ViolinColor', cols);
ylabel('Mean HR (bpm)');
%[~,ptt] = ttest(Task.Magnus.(foi), Rest.Magnus.(foi));
[ptt, ~, stats]  = signrank(Task.Magnus.(foi), Rest.Magnus.(foi));
%[~,ptt] = ttest(Task.Magnus.(foi), Rest.Magnus.(foi));
if ptt<0.001
    p='<0.001';
else
    p=num2str(round(ptt*1000)/1000);
end
text(1,median(Task.Magnus.(foi)),[num2str(round(median(Task.Magnus.(foi)))) ';' num2str(round(prctile(Task.Magnus.(foi),2.5))) '-' num2str(round(prctile(Task.Magnus.(foi),97.5)))]);
text(2,median(Rest.Magnus.(foi)),[num2str(round(median(Rest.Magnus.(foi)))) ';' num2str(round(prctile(Rest.Magnus.(foi),2.5))) '-' num2str(round(prctile(Rest.Magnus.(foi),97.5)))]);
text(1.5,median(Magnus_all_HR),[num2str(round(median(Magnus_all_HR))) ';' num2str(round(prctile(Magnus_all_HR,2.5))) '-' num2str(round(prctile(Magnus_all_HR,97.5)))]);
title(['Magnus, Wilcoxon signed rank p=' p ',zval=' num2str(stats.zval)]);
set(gca,'xticklabels',{'Task','Rest'},'ylim',y_lims);
axis square
[~, pLM] = lillietest(Task.Magnus.(foi) -Rest.Magnus.(foi));

subplot(1,2,2)
% Combine data
Bacchus_all_HR = [Task.Bacchus.(foi), Rest.Bacchus.(foi)];
% Create group labels
group = [ones(1, numel(Task.Bacchus.(foi))), ones(1, numel(Rest.Bacchus.(foi)))*2];
violinplot(Bacchus_all_HR, group, 'ViolinColor', cols);
ylabel('Mean HR (bpm)');
[ptt, ~, stats] = signrank(Task.Bacchus.(foi), Rest.Bacchus.(foi));
%[~,ptt] = ttest(Task.Magnus.(foi), Rest.Magnus.(foi));
if ptt<0.001
    p='<0.001';
else
    p=num2str(round(ptt*1000)/1000);
end
text(1,median(Task.Bacchus.(foi)),[num2str(round(median(Task.Bacchus.(foi)))) ';' num2str(round(prctile(Task.Bacchus.(foi),2.5))) '-' num2str(round(prctile(Task.Bacchus.(foi),97.5)))]);
text(2,median(Rest.Bacchus.(foi)),[num2str(round(median(Rest.Bacchus.(foi)))) ';' num2str(round(prctile(Rest.Bacchus.(foi),2.5))) '-' num2str(round(prctile(Rest.Bacchus.(foi),97.5)))]);
text(1.5,median(Bacchus_all_HR),[num2str(round(median(Bacchus_all_HR))) ';' num2str(round(prctile(Bacchus_all_HR,2.5))) '-' num2str(round(prctile(Bacchus_all_HR,97.5)))]);
title(['Bacchus,  Wilcoxon signed rank  p=' p ',zval=' num2str(stats.zval)]);
set(gca,'xticklabels',{'Task','Rest'},'ylim',y_lims);
axis square
[~, pLB] = lillietest(Task.Bacchus.(foi) -Rest.Bacchus.(foi));
mtit(['HR over all: ' num2str(round(median([Bacchus_all_HR, Magnus_all_HR]))) ';' num2str(round(prctile([Bacchus_all_HR, Magnus_all_HR],2.5))) '-' num2str(round(prctile([Bacchus_all_HR, Magnus_all_HR],97.5)))]);

results_file='Y:\Projects\Pulv_bodysignal\Figures\Suppfig. 3-Heartrate Violins';
wanted_size=[50 30];
set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    %
export_fig(h, results_file, '-pdf'); %% how come this does not export most plots ??

%% Histograms -BPM
h=figure;
subplot(1,2,1)
% Plot first histogram
bins=100:10:180;
h1 = histogram(Task.Magnus.(foi), bins,'Normalization', 'probability', ...
    'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;

% Plot second histogram
h2 = histogram(Rest.Magnus.(foi), bins,'Normalization', 'probability', ...
    'FaceColor', 'blue', 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Customize plot
xlabel('Heart Rate (bpm)');
ylabel('Probability');
title({'Monkey M';'Heart Rate Distribution:'});
legend({'Task', 'Rest'});
axis square

subplot(1,2,2)
% Plot first histogram
h2 = histogram(Task.Bacchus.(foi), bins,'Normalization', 'probability', ...
    'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
h1 = histogram(Rest.Bacchus.(foi), bins,'Normalization', 'probability', ...
    'FaceColor', 'blue', 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Customize plot
xlabel('Heart Rate (bpm)');
ylabel('Probability');
title({'Monkey B';'Heart Rate Distribution:'});
legend({'Task', 'Rest'});
axis square

results_file='Y:\Projects\Pulv_bodysignal\Figures\Suppfig. 3-Heartrate Histograms';
wanted_size=[50 30];
set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    
export_fig(h, results_file, '-pdf'); 

%% MAX ECG

Magnus = load('Y:\Projects\Pulv_bodysignal\ECG\ECG_Magnus_TaskRest_generalTrig_ECG\Magnus_ecgMaxMin_perSite.mat');
Bacchus = load('Y:\Projects\Pulv_bodysignal\ECG\ECG_Bacchus_TaskRest_generalTrig_ECG\Bacchus_ecgMaxMin_perSite.mat');

h=figure;
subplot(121)
% Combine data
Magnus_MaxTaskRest = [Magnus.ecgMaxMin.Rest_max, Magnus.ecgMaxMin.Task_max];
% Create group labels
group = [ones(1, length(Magnus.ecgMaxMin)), repmat(2, 1, length(Magnus.ecgMaxMin))];
violinplot(Magnus_MaxTaskRest, group);
ylabel('Max Task & Rest ecg');
title('Magnus')
axis square

subplot(122)
% Combine data
Bacchus_MaxTaskRest = [Bacchus.ecgMaxMin.Rest_max, Bacchus.ecgMaxMin.Task_max];
% Create group labels
group = [ones(1, length(Bacchus.ecgMaxMin)), repmat(2, 1, length(Bacchus.ecgMaxMin))];
violinplot(Bacchus_MaxTaskRest, group);
ylabel('Max Task & Rest ecg');
title('Bacchus')
axis square

results_file='Y:\Projects\Pulv_bodysignal\Figures\Suppfig. 3-max ecg distribution';
wanted_size=[50 30];
set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    %
export_fig(h, results_file, '-pdf'); 

clear all
clc
end