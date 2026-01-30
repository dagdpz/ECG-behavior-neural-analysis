close all, 
clear all
clc

% load('Y:\Data\Magnus_phys_mat_from_TDT\20230614\MagTDT2023-06-14_02.mat')
load('Y:\Data\Bacchus_phys_mat_from_TDT\20210720\BacTDT2021-07-20_03.mat')

%%
Rpeaks = [];
for t = 1:length(TDT_trial)
    %     sampleECG = cat(3,sampleECG,TDT_trial(t).ECG1(1:nsamples));
    %     plot(time, TDT_trial(t).ECG1(1:nsamples))
    [pks,locs] = findpeaks(TDT_trial(t).ECG1, 'MinPeakDistance',500,'MinPeakHeight',1.5*10^-3);
    Rpeaks{t}= locs;
end

wnd = 0.550;
nsamples = floor(wnd*TDT_trial(1).ECG1_samplingrate);
sampleECG = [];
for t = 1:length(TDT_trial)
    numRpeaks = length(Rpeaks{t});
    for peak = 1:numRpeaks
        if  (Rpeaks{1, t}(peak)+nsamples) <= length(TDT_trial(t).ECG1) && ( Rpeaks{1, t}(peak)-nsamples)>=1
            sampleECG{t,peak} = TDT_trial(t).ECG1( Rpeaks{1, t}(peak)-nsamples : Rpeaks{1, t}(peak)+nsamples);
        else
            continue;
            
        end
    end
end
idx = cell2mat(cellfun(@(x) ~isempty(x),sampleECG,'uniformoutput',false));
validSampl = sampleECG(idx);
concatECG = [];
for i = 1:length(validSampl)
concatECG(i,:) = validSampl{i,1};
end
timeAll = -wnd:1/TDT_trial(1).ECG1_samplingrate:wnd-1/TDT_trial(1).ECG1_samplingrate;
strtIDX = find(timeAll>-0.3);
time2plot = timeAll(strtIDX(1)-1:end);
lineProps={'color',[0.1 0.1 0.1]};
meanECG = mean(concatECG(:,strtIDX(1)-1:end),1);
stdECG = std(concatECG(:,strtIDX(1)-1:end),[],1);

h = figure;
h.WindowState = 'maximized';
shadedErrorBar(time2plot,meanECG,stdECG,lineProps,1);
xline(0);
title('Analysis Window (0.5 sec)','Interpreter','none');
xlabel('Time from R-Peak (s)');
ax = gca;
ax.YAxis.Visible = 'off';   % removes the y-axis fully
yl = ylim;  
rectangle('Position', [-0.25, yl(1), 0.25 - (-0.25), yl(2) - yl(1)], ...
          'EdgeColor', [0.5,0.5,0.5], 'LineWidth', 1);
figname = fullfile(['Y:\Projects\Pulv_bodysignal\LFP\Figures',filesep,'Fig1_PanelB']);
export_fig(h,[figname,'.pdf']);

