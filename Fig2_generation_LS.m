function Fig2_generation_LS
clc
%%
Site_ID = 'Bac_20220315_Site_15';%'Bac_20211007_Site_01';%'Bac_20210730_Site_03';
project = 'Pulv_bodysignal';
version = 'ECG_Bacchus_complete';
driver_path = 'Y:';%'/home/shamim/fileserver';
%showfiltered='filtered';

cfg = [];
cfg.ecg.timestep = 1; %??

cfg.process_ecg=0;
cfg.project = project;

cfg.fldr.root     = [driver_path,filesep,'Projects',filesep,cfg.project];
cfg.results_folder = [driver_path,filesep,'Projects',filesep,cfg.project];
ecg_bna_location     =which('ecg_bna_define_folders');
github_folder        =ecg_bna_location(1:strfind(ecg_bna_location,['ECG-behavior-neural-analysis' filesep 'ecg_bna_define_folders'])-1);
cfg.version = version;
run([github_folder filesep 'Settings' filesep cfg.project filesep 'ECG_bna' filesep cfg.version '.m']);
cfg = ecg_bna_define_folders(cfg);
%
% filteredornot={'unfiltered','filtered'};
% for flt=1:numel(filteredornot)
%     showfiltered=filteredornot{flt};



%% 'LFP'
%     base_path=['Y:\Projects\Pulv_bodysignal\LFP\' ECG_Bacchus_Complete '\'];
%     load([base_path '\Per_Site_Smpl\' Site_ID '.mat'])
%     LFPEV=triggered_site_data;
base_path=['Y:\Projects\Pulv_bodysignal\LFP\' version];
load([base_path '\Per_Site\' Site_ID '.mat'])
LFP=triggered_site_data;
%% 'MUA'
base_path=['Y:\Projects\Pulv_bodysignal\MUA\' version] ;
load([base_path '\Per_Site\' Site_ID '.mat'])
MUA=triggered_site_data;

c=2; %task?
e=1;

%% time
sr=1017; %hz
sl=1/sr*10;
time=[-sl*26:sl:sl*26];

%% frequencies !
% % potential cutoff -> probably startfreq should go into settings at somepoint
% startfreq=3.5;
% % if strcmp(showfiltered,'unfiltered')
% %     startfreq=0;
% % end
% %colsbp=jet(length(cfg.lfp.frequency_bands));
% %freqName = cfg.lfp.freqName(fbx:end);
% %colsbp=colsbp(fbx:end,:);
% %freqb=num2cell(strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz'),2);
% fx=find(cfg.lfp.foi>=startfreq,1,'first');
% freq = cfg.lfp.foi(fx:end);

startfreq=cfg.lfp.foi(1);
fbx=find(cfg.lfp.frequency_bands(:,1)>=startfreq,1,'first');
frequency_bands=cfg.lfp.frequency_bands(fbx:end,:);
%frequency_bands=cfg.lfp.frequency_bands; %(fbx:end,:);
freq = cfg.lfp.foi;


%% start plotting
h = figure('units','normalized','position',[0 0 1 1]);
methods= {'observed','surrogate','normalized'};
titles= {'observed','surrogate','normalized'};
cols= {[0 0 1],[1 0 0],[0.2 0.7 0.1]};
for mt = 1: numel(methods)
    M=methods{mt};
    %% LFP evoked
    spLFP(mt)=subplot(4, 3, mt);
    TP=LFP.condition(c).event(e);
    TPG=TP.(M);
    TPS=TP.significance;
    hold on;
    
    % error bar
    lineprops={'color',cols{mt}};
    switch M
        case 'surrogate'
            conf=squeeze([TPG.evoked.mean-TPG.evoked.conf95(1,:,:);-1*(TPG.evoked.mean-TPG.evoked.conf95(2,:,:))]);
        case 'observed'
            %conf=[TPG.evoked.mean-TPG.evoked.conf95(1,:,:);-1*(TPG.evoked.mean-TPG.evoked.conf95(2,:,:))];
            conf=squeeze(TPG.evoked.sterr)';
        case 'normalized'
            conf=squeeze(LFP.condition(c).event(e).surrogate.evoked.sterr)'./LFP.condition(c).event(e).surrogate.evoked.std;
    end
    meantp=squeeze(TPG.evoked.mean)';
    shadedErrorBar(time,meantp,conf,lineprops,1);
    %
    %         % add surrogate
    %         if strcmp(M,'observed')
    %             lineprops={'color',cols{2}};
    %             meantp=squeeze(TP.surrogate.evoked.mean)';
    %             surrogate_conf=[TP.surrogate.evoked.mean-TP.surrogate.evoked.conf95(1,:,:);-1*(TP.surrogate.evoked.mean-TP.surrogate.evoked.conf95(2,:,:))];
    %             shadedErrorBar(time,meantp,squeeze(surrogate_conf),lineprops,1);
    %         end
    
    xlabel('Time(s)');
    
    if strcmp(M,'normalized')
        ylabel('normalized LFP');
    else
        ylabel('Voltage (mV)');
    end
    title(titles{mt},'Interpreter', 'none');
    if strcmp(M,'surrogate')
        ylm = get(spLFP(1),'Ylim');
        set(gca,'Ylim',ylm);
    else
        ylm = get(gca,'Ylim');
    end
    line([0 0], ylm, 'color', 'k');
    if strcmp(M,'normalized')
        significance = double(squeeze(TPS.evoked));
        significance(significance==0)=NaN;
        significance=significance.*ylm(1);
        plot(time,significance','linewidth',3);
    end
    set(gca, 'xlim', [time(1) time(end)]);
    axis square
    
    %% POW
    TP=LFP.condition(c).event(e);
    TPG=TP.(M);
    TPS=TP.significance;
    spP(mt)=subplot(4, 3, mt+numel(methods));
    hold on;
    %toplot=squeeze(TPG.pow.mean(:,fx:end,:));
    toplot=squeeze(TPG.pow.mean);
    
    image(time, 1:numel(freq), toplot,'CDataMapping','scaled');
    set(gca,'YDir','normal');
    line([0 0], ylim, 'color', 'k');
    
    if strcmp(M,'normalized')
        %significance = double(squeeze(TPS.pow(:,fx:end,:)));
        significance = double(squeeze(TPS.pow));
        % plotting a contour around the significant parts:
        
        ecg_bna_draw_outlines(significance==1,'r')
        ecg_bna_draw_outlines(significance==-1,'b')
        %contour(time,1:numel(freq),significance,1,'linecolor','k')
    end
    nonnan=toplot;nonnan(isnan(nonnan))=[];
    collimP{mt}=[min(nonnan(:)) max(nonnan(:))];
    
    % horizontal lines to separate frequency bands
    fbandstart = unique(frequency_bands(:))';
    fbandstart_idx = zeros(size(fbandstart));
    for f = fbandstart
        f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
        line([time(1) time(end)], [f_idx f_idx], 'color', 'k', 'linestyle', '--');
        fbandstart_idx(fbandstart == f) = f_idx;
    end
    
    set(gca,'TickDir','out')
    set(gca, 'ytick', fbandstart_idx);
    set(gca, 'yticklabel', fbandstart);
    set(gca, 'ylim', [0.5,numel(freq) + 0.5]);
    set(gca, 'xlim', [time(1) time(end)]);
    ylabel('Frequency (Hz)');
    axis square
    
    %% ITPC
    TP=LFP.condition(c).event(e);
    TPG=TP.(M);
    TPS=TP.significance;
    spI(mt)=subplot(4, 3, mt+2*numel(methods));
    hold on;
    %toplot=squeeze(TPG.itpc.mean(:,fx:end,:));
    toplot=squeeze(TPG.itpc.mean);
    
    image(time, 1:numel(freq), toplot,'CDataMapping','scaled');
    set(gca,'YDir','normal');
    line([0 0], ylim, 'color', 'k');
    
    if strcmp(M,'normalized')
        %significance = double(squeeze(TPS.itpc(:,fx:end,:)));
        significance = double(squeeze(TPS.itpc));
        contour(time,1:numel(freq),significance,1,'linecolor','k')
    end
    nonnan=toplot;nonnan(isnan(nonnan))=[];
    collimI{mt}=[min(nonnan(:)) max(nonnan(:))];
    
    % horizontal lines to separate frequency bands
    fbandstart = unique(frequency_bands(:))';
    fbandstart_idx = zeros(size(fbandstart));
    for f = fbandstart
        f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
        line([time(1) time(end)], [f_idx f_idx], 'color', 'k', 'linestyle', '--');
        fbandstart_idx(fbandstart == f) = f_idx;
    end
    
    set(gca,'TickDir','out')
    set(gca, 'ytick', fbandstart_idx);
    set(gca, 'yticklabel', fbandstart);
    set(gca, 'ylim', [0.5,numel(freq) + 0.5]);
    set(gca, 'xlim', [time(1) time(end)]);
    ylabel('Frequency (Hz)');
    axis square
    
    %% MUA
    TP=MUA.condition(c).event(e);
    TPG=TP.(M);
    TPS=TP.significance;
    spMUA(mt)=subplot(4, 3, mt+3*numel(methods));
    hold on;
    
    % error bar
    lineprops={'color',cols{mt}};
    TPGM=TPG.evoked.mean;
    switch M
        case 'surrogate'
            conf=squeeze([TPGM-TPG.evoked.conf95(1,:,:);-1*(TPGM-TPG.evoked.conf95(2,:,:))]);
        case 'observed'
            %conf=[TPGM-TPG.evoked.conf95(1,:,:);-1*(TPGM-TPG.evoked.conf95(2,:,:))];
            conf=squeeze(TPG.evoked.sterr)';
        case 'normalized'
            
            conf=squeeze(MUA.condition(c).event(e).surrogate.evoked.sterr)'./MUA.condition(c).event(e).surrogate.evoked.std;
            %conf=squeeze(TPG.evoked.sterr)';
    end
    meantp=squeeze(TPGM)';
    shadedErrorBar(time,meantp,conf,lineprops,1);
    
    xlabel('Time(s)');
    if strcmp(M,'normalized')
        ylabel('normalized MUA');
    else
        ylabel('Voltage (mV)');
    end
    %         if strcmp(M,'observed')
    %             lineprops={'color',cols{2}};
    %             meantp=squeeze(TP.surrogate.evoked.mean)';
    %             surrogate_conf=[TP.surrogate.evoked.mean-TP.surrogate.evoked.conf95(1,:,:);-1*(TP.surrogate.evoked.mean-TP.surrogate.evoked.conf95(2,:,:))];
    %             shadedErrorBar(time, meantp,squeeze(surrogate_conf),lineprops,1);
    %         end
    if strcmp(M,'surrogate')
        ylm = get(spMUA(1),'Ylim');
        set(gca,'Ylim',ylm);
    else
        ylm = get(gca,'Ylim');
    end
    line([0 0], ylm, 'color', 'k');
    if strcmp(M,'normalized')
        significance = double(squeeze(TPS.mua));
        significance(significance==0)=NaN;
        significance=significance.*ylm(1);
        plot(time,significance','linewidth',3);
    end
    set(gca, 'xlim', [time(1) time(end)]);
    axis square
end

spP_pos=[1,2];
for s=spP_pos
    maxP=max([collimP{spP_pos}]);
    cm=jet(255);
    cm=cm(128:end,:);
    colormap(spP(s),cm);
    cb = colorbar('peer',spP(s));
    set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
    set(spP(s),'CLim',[0 maxP]);
end

spI_pos=[1,2];
for s=spI_pos
    maxI=max([collimI{spI_pos}]);
    cm=jet(255);
    cm=cm(128:end,:);
    colormap(spI(s),cm);
    cb = colorbar('peer',spI(s));
    set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
    set(spI(s),'CLim',[0 maxI]);
end

maxI=max([collimI{3}]);
cm=jet(255);
cm=cm(128:end,:);
colormap(spI(3),cm);
cb = colorbar('peer',spI(3));
set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
set(spI(3),'CLim',[0 maxI]);

cm=jet(128);
colormap(spP(3),cm);
cb = colorbar('peer',spP(3));
set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
%set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
%set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
M=max(abs([min([collimP{3}]),max([collimP{3}])]));
set(spP(3),'CLim',[-M M]);
%results_file=['Y:\Projects\Pulv_bodysignal\Figures\' Site_ID '_' showfiltered];
results_file=['Y:\Projects\Pulv_bodysignal\Figures\' Site_ID];
wanted_size=[50 30];
set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    %
export_fig(h, results_file, '-pdf'); %% how come this does not export most plots ??
close all
end






