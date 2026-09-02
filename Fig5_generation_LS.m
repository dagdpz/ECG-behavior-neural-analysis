function Fig5_generation_LS(figstoplot, withunits, version, project, driver_path)
if nargin < 1 || isempty(figstoplot), figstoplot = [8]; end
if nargin < 2 || isempty(withunits), withunits = 'all'; end
if nargin < 3 || isempty(version), version = 'ECG_Bacchus_complete'; end
if nargin < 4 || isempty(project), project = 'Pulv_bodysignal'; end
if nargin < 5 || isempty(driver_path), driver_path = 'Y:'; end

%%
cfg.ecg.timestep = 1; %??

cfg.project = project;
cfg.results_folder = [driver_path,filesep,'Projects',filesep,cfg.project];
ecg_bna_location     =which('ecg_bna_define_folders');
github_folder        =ecg_bna_location(1:strfind(ecg_bna_location,['ECG-behavior-neural-analysis' filesep 'ecg_bna_define_folders'])-1);
cfg.version = version;
run([github_folder filesep 'Settings' filesep cfg.project filesep 'ECG_bna' filesep cfg.version '.m']);
% cfg.fldr.root=
% cfg = ecg_bna_define_folders(cfg);



%% Bacchus
% LFP
file=['Y:\Projects\Pulv_bodysignal\LFP\ECG_Bacchus_complete\Bacchus_Rpeak_Triggered_target_wise_Grand_grand_avg_sessions_targets' withunits '.mat'];
D.LFP.B=load(file,'tar');
% MUA
file=['Y:\Projects\Pulv_bodysignal\MUA\ECG_Bacchus_complete\Grand_avg_mua_targets_' withunits '.mat'];
D.MUA.B=load(file,'tar');
%% Magnus
% LFP
file=['Y:\Projects\Pulv_bodysignal\LFP\ECG_Magnus_complete\Magnus_Rpeak_Triggered_target_wise_Grand_grand_avg_sessions_targets' withunits '.mat'];
D.LFP.M=load(file,'tar');
% MUA
file=['Y:\Projects\Pulv_bodysignal\MUA\ECG_Magnus_complete\Grand_avg_mua_targets_' withunits '.mat'];
D.MUA.M=load(file,'tar');

sr=1017.2526041; %hz
sl=1/sr*10;
lfp_time=[-sl*26:sl:sl*26];
ticklabelsx=[-0.2,0,0.2];
for tl=1:numel(ticklabelsx)
    tt=ticklabelsx(tl);
    tlstepsize=numel(lfp_time)/(lfp_time(end)-lfp_time(1));
    ticksx(tl)=tt*tlstepsize-lfp_time(1)*tlstepsize;
end



E='R';
monkeys={'B','M'};
targets={'VPL','dPul','MD'};
datatype={'LFP','MUA'};
datafield={'evoked','evoked'};

%% colormaps

logscale_127=exp((0:126)'/127);
logscale_127=(logscale_127-1)/(max(logscale_127)-1);
%logscale_127=(0:126)'/127;
bluered_colormap=[logscale_127*([255 255 255]-[0 0 255])+repmat([0 0 255],size(logscale_127)); 255 255 255;...
    flipud(logscale_127*([255 255 255]-[255 0 0])+repmat([255 0 0],size(logscale_127)))]/255;


%% FIG 4
if ismember(4,figstoplot)
    h = figure('units','normalized','position',[0 0 1 1]);
    y_lim={[-15 15] [-3 5] [-3 5] [-4 3]};
    ylabels={'Normalized LFP','Normalized MUA'};
    for m=1:numel(monkeys)
        M=monkeys{m};
        for y=1:numel(datatype)
            Y=datatype{y};
            tartmp=D.(Y).(M).tar;
            curr_ylim=y_lim{(m-1)*2+y};
            
            for t=1:numel(targets)
                T=targets{t};
                [~,Ti]=ismember({tartmp.target},T);
                subplot(2,6,(m-1)*6+(y-1)*3+t)
                xlim([min(lfp_time) max(lfp_time)]);
                ylim(curr_ylim);
                ylabel(ylabels{y});
                xlabel('Time relative to R peak (s)');
                line([0 0], curr_ylim,'color','k')
                axis square
                hold on
                title([T ' (' num2str(tartmp(Ti==1).nSites) ' sites)'],'interpreter','none')
                
                for c=1:2
                    %cname=tartmp(Ti==1).con(c).cond_name;
                    col=cfg.condition(c).color;
                    f=datafield{y};
                    mean_tp=  smoothit(tartmp(Ti==1).con(c).(E).(f));
                    EB_tp=  smoothit(tartmp(Ti==1).con(c).(E).([f '_sterr']));
                    
                    lineProps={'color',col};
                    shadedErrorBar(lfp_time,mean_tp,EB_tp,lineProps,1);
                end
            end
        end
    end
    results_file=['Y:\Projects\Pulv_bodysignal\Figures\Fig4-mualfp_grandaverage' '_' withunits];
    wanted_size=[50 30];
    set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    %
    export_fig(h, results_file, '-pdf'); %% how come this does not export most plots ??
end

%% FIG5
if ismember(5,figstoplot)
    h = figure('units','normalized','position',[0 0 1 1]);
    y_lim={[0 100] [0 100] [0 50] [0 5]};
    ylabels={'sig. LFP (% of sites)','sig. MUA (% of sites)'};
    for m=1:numel(monkeys)
        M=monkeys{m};
        for y=1:numel(datatype)
            Y=datatype{y};
            tartmp=D.(Y).(M).tar;
            
            for t=1:numel(targets)
                if y==2 && m==2 && t==1
                    curr_ylim=[0 50];
                else
                    curr_ylim=y_lim{(m-1)*2+y};
                end
                
                T=targets{t};
                [~,Ti]=ismember({tartmp.target},T);
                subplot(2,6,(m-1)*6+(y-1)*3+t)
                xlim([min(lfp_time) max(lfp_time)]);
                ylim(curr_ylim);
                ylabel(ylabels{y});
                xlabel('Time relative to R peak (s)');
                line([0 0], curr_ylim,'color','k')
                axis square
                hold on
                title([T ' (' num2str(tartmp(Ti==1).nSites) ' sites)'],'interpreter','none')
                
                for c=1:2
                    %cname=tartmp(Ti==1).con(c).cond_name;
                    col=cfg.condition(c).color;
                    f=datafield{y};
                    Sig=  smoothit(tartmp(Ti==1).con(c).(E).([f '_sig']));
                    p=  tartmp(Ti==1).con(c).(E).([f '_p']);
                    if p<0.001
                        P='p<0.001';
                    else
                        P=['p=' num2str(round(p*1000)/1000)];
                    end
                    plot(lfp_time,Sig,'color',col,'linewidth',1.5);
                    % peak latency
                    [pkv,pki]=max(Sig(lfp_time>0));
                    pkt=lfp_time(pki+sum(lfp_time<=0));
                    sepy=diff(curr_ylim)/70;
                    plot(pkt,pkv+sepy,'color',col,'Marker','v');
                    text(pkt-0.02,pkv+sepy,[num2str(round(pkt*100)*10) 'ms'],'color',col,'HorizontalAlignment', 'right');
                    text(lfp_time(2),curr_ylim(1)+diff(curr_ylim)/20*c,num2str(P),'color',col);
                end
            end
        end
    end
    results_file=['Y:\Projects\Pulv_bodysignal\Figures\Fig5-%significant' '_' withunits];
    wanted_size=[50 30];
    set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    %
    export_fig(h, results_file, '-pdf');
end







%% frequencies !

% potential cutoff -> probably startfreq should go into settings at somepoint
%
startfreq=3.5;
colsbp=cool(length(cfg.lfp.frequency_bands));
fbx=find(cfg.lfp.frequency_bands(:,1)>=startfreq,1,'first');
frequency_bands=cfg.lfp.frequency_bands(fbx:end,:);
%freqName = cfg.lfp.freqName(fbx:end);
colsbp=colsbp(fbx:end,:);
freqb=num2cell(strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz'),2);
fx=find(cfg.lfp.foi>=startfreq,1,'first');
freq = cfg.lfp.foi(fx:end);





if ismember(6,figstoplot)
    h = figure('units','normalized','position',[0 0 1 1]);
    y_lim={[-15 15] [-3 3] [-3 5] [-4 3]};
    for m=1:numel(monkeys)
        M=monkeys{m};
        tartmp=D.LFP.(M).tar;
        %curr_ylim=y_lim{(m-1)*2+y};
        for t=1:numel(targets)
            T=targets{t};
            [~,Ti]=ismember({tartmp.target},T);
            for c=1:2
                %% POWER
                sph(t,m,c,1)=subplot(6,4,(t-1)*8+(m-1)*4+(c-1)*2 +1);
                ax=sph(t,m,c,1);
                cname=tartmp(Ti==1).con(c).cond_name;
                if m==1 && t==1
                    title(['POW-' cname])
                end
                if c==1
                    ylabel({['Monkey ' M],[T, ',' num2str(tartmp(Ti==1).nSites) ' sites'],'Frequency (Hz)'});
                else
                    ylabel('Frequency (Hz)');
                end
                if m==numel(monkeys) && t==numel(targets)
                    xlabel('Time relative to R peak (s)');
                end
                con=tartmp(Ti==1).con(c).(E);
                toplot=con.pow;
                sigplot=con.pow_popsig;
                
                %xlim([min(lfp_time) max(lfp_time)]);
                xlim([0 size(toplot,2)]);
                axis square
                hold on
                
                %image(lfp_time-sl/2,1:size(toplot,1),toplot,'CDataMapping','scaled');
                image(toplot,'CDataMapping','scaled');
                set(gca,'YDir','normal');
                significance = double(sigplot);
                %                 contour(lfp_time-sl/2,(1:size(toplot,1)),significance==1,1,'linecolor','r')
                %                 contour(lfp_time-sl/2,(1:size(toplot,1)),significance==-1,1,'linecolor','b')
                
                ecg_bna_draw_outlines(significance==1,'r')
                ecg_bna_draw_outlines(significance==-1,'b')
                
                % horizontal lines to separate frequency bands
                fbandstart = unique(frequency_bands(:))';
                fbandstart_idx = zeros(size(fbandstart));
                for f = fbandstart
                    f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                    line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                    fbandstart_idx(fbandstart == f) = f_idx;
                end
                line([ticksx(ticklabelsx==0) ticksx(ticklabelsx==0)], ylim, 'color', 'k');
                
                set(ax,'TickDir','out')
                set(ax, 'ytick', fbandstart_idx);
                set(ax, 'yticklabel', fbandstart);
                set(ax, 'xtick', ticksx);
                set(ax, 'xticklabel', ticklabelsx);
                
                
                set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                box on
                
                %% POWER sig
                sph(t,m,c,2)=subplot(6,4,(t-1)*8+(m-1)*4+(c-1)*2 +2);
                ax=sph(t,m,c,2);
                hold on;
                cname=tartmp(Ti==1).con(c).cond_name;
                if m==1 && t==1
                    title(['Sig. POW-' cname])
                end
                if m==numel(monkeys) && t==numel(targets)
                    xlabel('Time relative to R peak (s)');
                end
                toplot=abs(con.powbp_sig);
                toplot=smoothit(toplot);
                ylabel('Proportion of sites (%)');
                
                y_lim=[0 ceil(max(toplot(:))/10)*10];
                sepy=diff(y_lim)/20;
                set(ax,'ylim',y_lim);
                for bpf=1:size(toplot,1)
                    col=colsbp(bpf,:);
                    Sig=toplot(bpf,:);
                    [pkv,pki]=max(Sig(lfp_time>=0));
                    pkt=lfp_time(pki+sum(lfp_time<0));
                    if pkv>(y_lim(end)-sepy)
                        plot(pkt,pkv-sepy,'color',col,'Marker','^');
                        text(pkt-0.02,pkv-sepy,[num2str(round(pkt*100)*10) 'ms'],'color',col,'HorizontalAlignment', 'right');
                    else
                        plot(pkt,pkv+sepy,'color',col,'Marker','v');
                        text(pkt-0.02,pkv+sepy,[num2str(round(pkt*100)*10) 'ms'],'color',col,'HorizontalAlignment', 'right');
                    end
                end
                %
                %                 for bpf=1:size(toplot,1)
                %                     col=colsbp(bpf,:);
                %                     Sig=toplot(bpf,:);
                %                     [pkv,pki]=max(Sig(lfp_time>=0));
                %                     pkt=lfp_time(pki+sum(lfp_time<0));
                %                     sepy=5;
                %                     plot(pkt,pkv+sepy,'color',col,'Marker','v');
                %                     text(pkt-0.02,pkv+sepy,[num2str(round(pkt*100)*10) 'ms'],'color',col,'HorizontalAlignment', 'right');
                %                 end
                
                xlim([min(lfp_time) max(lfp_time)]);
                axis square
                set(ax,'ColorOrder',colsbp);
                plot(lfp_time,toplot')
                if m==1 && t==1 && c==1
                    legend(freqb,'fontsize',3);
                end
            end
        end
    end
    
    %% adjusting limits
    for m=1:numel(monkeys)
        for t=1:numel(targets)
            y_lim=[0,max([get(sph(t,m,1,2),'ylim') get(sph(t,m,2,2),'ylim')])];
            c_lim=max(abs([get(sph(t,m,1,1),'clim'), get(sph(t,m,2,1),'clim')]));
            for c=1:2
                % sig %
                ax=sph(t,m,c,2);
                set(ax,'ylim',y_lim);
                subplot(ax);
                line([0 0], y_lim, 'color', 'k');
                % POW
                ax=sph(t,m,c,1);
                subplot(ax);
                %cm=jet(128);
                cm=customcolormap_preset('red-yellow-blue',256);
                %cm=bluered_colormap;
                colormap(ax,cm);
                cb = colorbar;
                set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                set(get(cb,'title'),'string', 'n. Pow', 'fontsize',8);
                set(ax,'clim',[-c_lim c_lim]);
            end
        end
    end
    
    results_file=['Y:\Projects\Pulv_bodysignal\Figures\Fig6-POWER_' withunits];
    wanted_size=[50 30];
    set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    %
    export_fig(h, results_file, '-pdf'); %% how come this does not export most plots ??
end



if ismember(6.5,figstoplot)
    h = figure('units','normalized','position',[0 0 1 1]);
    y_lim={[-15 15] [-3 3] [-3 5] [-4 3]};
    for m=1:numel(monkeys)
        M=monkeys{m};
        tartmp=D.LFP.(M).tar;
        %curr_ylim=y_lim{(m-1)*2+y};
        for t=1:numel(targets)
            T=targets{t};
            [~,Ti]=ismember({tartmp.target},T);
            for c=1:2
                %% POWER
                sph(t,m,c,1)=subplot(6,4,(t-1)*8+(m-1)*4+(c-1)*2 +1);
                ax=sph(t,m,c,1);
                cname=tartmp(Ti==1).con(c).cond_name;
                if m==1 && t==1
                    title(['POW-' cname])
                end
                if c==1
                    ylabel({['Monkey ' M],[T, ',' num2str(tartmp(Ti==1).nSites) ' sites'],'Frequency (Hz)'});
                else
                    ylabel('Frequency (Hz)');
                end
                if m==numel(monkeys) && t==numel(targets)
                    xlabel('Time relative to R peak (s)');
                end
                con=tartmp(Ti==1).con(c).(E);
                toplot=con.pow;
                sigplot=con.pow_popsig;
                
                %xlim([min(lfp_time) max(lfp_time)]);
                xlim([0 size(toplot,2)]);
                axis square
                hold on
                
                %image(lfp_time-sl/2,1:size(toplot,1),toplot,'CDataMapping','scaled');
                image(toplot,'CDataMapping','scaled');
                set(gca,'YDir','normal');
                significance = double(sigplot);
                %                 contour(lfp_time-sl/2,(1:size(toplot,1)),significance==1,1,'linecolor','r')
                %                 contour(lfp_time-sl/2,(1:size(toplot,1)),significance==-1,1,'linecolor','b')
                
                ecg_bna_draw_outlines(significance==1,'r')
                ecg_bna_draw_outlines(significance==-1,'b')
                
                % horizontal lines to separate frequency bands
                fbandstart = unique(frequency_bands(:))';
                fbandstart_idx = zeros(size(fbandstart));
                for f = fbandstart
                    f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                    line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                    fbandstart_idx(fbandstart == f) = f_idx;
                end
                line([ticksx(ticklabelsx==0) ticksx(ticklabelsx==0)], ylim, 'color', 'k');
                
                set(ax,'TickDir','out')
                set(ax, 'ytick', fbandstart_idx);
                set(ax, 'yticklabel', fbandstart);
                set(ax, 'xtick', ticksx);
                set(ax, 'xticklabel', ticklabelsx);
                
                
                set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                box on
                
                %% POWER sig
                sph(t,m,c,2)=subplot(6,4,(t-1)*8+(m-1)*4+(c-1)*2 +2);
                ax=sph(t,m,c,2);
                hold on;
                cname=tartmp(Ti==1).con(c).cond_name;
                if m==1 && t==1
                    title(['Sig. POW-' cname])
                end
                if m==numel(monkeys) && t==numel(targets)
                    xlabel('Time relative to R peak (s)');
                end
                toplot1=con.powbp_sig_pos;
                toplot1=smoothit(toplot1);
                toplot2=con.powbp_sig_neg;
                toplot2=smoothit(toplot2);
                ylabel('Proportion of sites (%)');
                
                y_lim=[floor(min(toplot2(:))/10)*10 ceil(max(toplot1(:))/10)*10];
                sepy=diff(y_lim)/20;
                set(ax,'ylim',y_lim);
                for bpf=1:size(toplot1,1)
                    col=colsbp(bpf,:);
                    Sig=toplot1(bpf,:);
                    [pkv,pki]=max(Sig(lfp_time>=0));
                    pkt=lfp_time(pki+sum(lfp_time<0));
                    if pkv>(y_lim(end)-sepy)
                        plot(pkt,pkv-sepy,'color',col,'Marker','^');
                        text(pkt-0.02,pkv-sepy,[num2str(round(pkt*100)*10) 'ms'],'color',col,'HorizontalAlignment', 'right');
                    else
                        plot(pkt,pkv+sepy,'color',col,'Marker','v');
                        text(pkt-0.02,pkv+sepy,[num2str(round(pkt*100)*10) 'ms'],'color',col,'HorizontalAlignment', 'right');
                    end
                end
                %
                %                 for bpf=1:size(toplot,1)
                %                     col=colsbp(bpf,:);
                %                     Sig=toplot(bpf,:);
                %                     [pkv,pki]=max(Sig(lfp_time>=0));
                %                     pkt=lfp_time(pki+sum(lfp_time<0));
                %                     sepy=5;
                %                     plot(pkt,pkv+sepy,'color',col,'Marker','v');
                %                     text(pkt-0.02,pkv+sepy,[num2str(round(pkt*100)*10) 'ms'],'color',col,'HorizontalAlignment', 'right');
                %                 end
                
                xlim([min(lfp_time) max(lfp_time)]);
                axis square
                set(ax,'ColorOrder',colsbp);
                plot(lfp_time,toplot1')
                plot(lfp_time,toplot2')
                if m==1 && t==1 && c==1
                    legend(freqb,'fontsize',3);
                end
            end
        end
    end
    
    %% adjusting limits
    for m=1:numel(monkeys)
        for t=1:numel(targets)
            y_lim=[min([get(sph(t,m,1,2),'ylim') get(sph(t,m,2,2),'ylim')]),max([get(sph(t,m,1,2),'ylim') get(sph(t,m,2,2),'ylim')])];
            c_lim=max(abs([get(sph(t,m,1,1),'clim'), get(sph(t,m,2,1),'clim')]));
            for c=1:2
                % sig %
                ax=sph(t,m,c,2);
                set(ax,'ylim',y_lim);
                subplot(ax);
                line([0 0], y_lim, 'color', 'k');
                % POW
                ax=sph(t,m,c,1);
                subplot(ax);
                %cm=jet(128);
                cm=customcolormap_preset('red-yellow-blue',256);
                %cm=bluered_colormap;
                colormap(ax,cm);
                cb = colorbar;
                set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                set(get(cb,'title'),'string', 'n. Pow', 'fontsize',8);
                set(ax,'clim',[-c_lim c_lim]);
            end
        end
    end
    
    results_file=['Y:\Projects\Pulv_bodysignal\Figures\Fig6,5-POWER_' withunits];
    wanted_size=[50 30];
    set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    %
    export_fig(h, results_file, '-pdf'); %% how come this does not export most plots ??
end



if ismember(6.8,figstoplot)
    h = figure('units','normalized','position',[0 0 1 1]);
    y_lim={[-15 15] [-3 3] [-3 5] [-4 3]};
    for m=1:numel(monkeys)
        M=monkeys{m};
        tartmp=D.LFP.(M).tar;
        %curr_ylim=y_lim{(m-1)*2+y};
        for t=1:numel(targets)
            T=targets{t};
            [~,Ti]=ismember({tartmp.target},T);
            for c=1:2
                %% POWER
                sph(t,m,c,1)=subplot(6,4,(t-1)*8+(m-1)*4+(c-1)*2 +1);
                ax=sph(t,m,c,1);
                cname=tartmp(Ti==1).con(c).cond_name;
                if m==1 && t==1
                    title(['POW-' cname])
                end
                if c==1
                    ylabel({['Monkey ' M],[T, ',' num2str(tartmp(Ti==1).nSites) ' sites'],'Frequency (Hz)'});
                else
                    ylabel('Frequency (Hz)');
                end
                if m==numel(monkeys) && t==numel(targets)
                    xlabel('Time relative to R peak (s)');
                end
                con=tartmp(Ti==1).con(c).(E);
                toplot=con.pow;
                sigplot=con.pow_popsig;
                
                %xlim([min(lfp_time) max(lfp_time)]);
                xlim([0 size(toplot,2)]);
                axis square
                hold on
                
                %image(lfp_time-sl/2,1:size(toplot,1),toplot,'CDataMapping','scaled');
                image(toplot,'CDataMapping','scaled');
                set(gca,'YDir','normal');
                significance = double(sigplot);
                %                 contour(lfp_time-sl/2,(1:size(toplot,1)),significance==1,1,'linecolor','r')
                %                 contour(lfp_time-sl/2,(1:size(toplot,1)),significance==-1,1,'linecolor','b')
                
                ecg_bna_draw_outlines(significance==1,'r')
                ecg_bna_draw_outlines(significance==-1,'b')
                
                % horizontal lines to separate frequency bands
                fbandstart = unique(frequency_bands(:))';
                fbandstart_idx = zeros(size(fbandstart));
                for f = fbandstart
                    f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                    line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                    fbandstart_idx(fbandstart == f) = f_idx;
                end
                line([ticksx(ticklabelsx==0) ticksx(ticklabelsx==0)], ylim, 'color', 'k');
                
                set(ax,'TickDir','out')
                set(ax, 'ytick', fbandstart_idx);
                set(ax, 'yticklabel', fbandstart);
                set(ax, 'xtick', ticksx);
                set(ax, 'xticklabel', ticklabelsx);
                
                
                set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                box on
                
                %% POWER sig
                sph(t,m,c,2)=subplot(6,4,(t-1)*8+(m-1)*4+(c-1)*2 +2);
                ax=sph(t,m,c,2);
                hold on;
                cname=tartmp(Ti==1).con(c).cond_name;
                if m==1 && t==1
                    title(['bandpassed POW-' cname])
                end
                if m==numel(monkeys) && t==numel(targets)
                    xlabel('Time relative to R peak (s)');
                end
                toplot1=con.powbp;
                toplot1=smoothit(toplot1);
                ylabel('Average normalized power (%)');
                
                y_lim=[floor(min(toplot1(:))) ceil(max(toplot1(:)))];
                sepy=diff(y_lim)/20;
                set(ax,'ylim',y_lim);
                for bpf=1:size(toplot1,1)
                    col=colsbp(bpf,:);
                    Sig=toplot1(bpf,:);
                    [pkv,pki]=max(Sig(lfp_time>=0));
                    pkt=lfp_time(pki+sum(lfp_time<0));
                    if pkv>(y_lim(end)-sepy)
                        plot(pkt,pkv-sepy,'color',col,'Marker','^');
                        text(pkt-0.02,pkv-sepy,[num2str(round(pkt*100)*10) 'ms'],'color',col,'HorizontalAlignment', 'right');
                    else
                        plot(pkt,pkv+sepy,'color',col,'Marker','v');
                        text(pkt-0.02,pkv+sepy,[num2str(round(pkt*100)*10) 'ms'],'color',col,'HorizontalAlignment', 'right');
                    end
                end
                %
                %                 for bpf=1:size(toplot,1)
                %                     col=colsbp(bpf,:);
                %                     Sig=toplot(bpf,:);
                %                     [pkv,pki]=max(Sig(lfp_time>=0));
                %                     pkt=lfp_time(pki+sum(lfp_time<0));
                %                     sepy=5;
                %                     plot(pkt,pkv+sepy,'color',col,'Marker','v');
                %                     text(pkt-0.02,pkv+sepy,[num2str(round(pkt*100)*10) 'ms'],'color',col,'HorizontalAlignment', 'right');
                %                 end
                
                xlim([min(lfp_time) max(lfp_time)]);
                axis square
                set(ax,'ColorOrder',colsbp);
                plot(lfp_time,toplot1')
                if m==1 && t==1 && c==1
                    legend(freqb,'fontsize',3);
                end
            end
        end
    end
    
    %% adjusting limits
    for m=1:numel(monkeys)
        for t=1:numel(targets)
            y_lim=[min([get(sph(t,m,1,2),'ylim') get(sph(t,m,2,2),'ylim')]),max([get(sph(t,m,1,2),'ylim') get(sph(t,m,2,2),'ylim')])];
            c_lim=max(abs([get(sph(t,m,1,1),'clim'), get(sph(t,m,2,1),'clim')]));
            for c=1:2
                % sig %
                ax=sph(t,m,c,2);
                set(ax,'ylim',y_lim);
                subplot(ax);
                line([0 0], y_lim, 'color', 'k');
                % POW
                ax=sph(t,m,c,1);
                subplot(ax);
                %cm=jet(128);
                cm=customcolormap_preset('red-yellow-blue',256);
                %cm=bluered_colormap;
                colormap(ax,cm);
                cb = colorbar;
                set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                set(get(cb,'title'),'string', 'n. Pow', 'fontsize',8);
                set(ax,'clim',[-c_lim c_lim]);
            end
        end
    end
    
    results_file=['Y:\Projects\Pulv_bodysignal\Figures\Fig6,8-POWER_' withunits];
    wanted_size=[50 30];
    set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    %
    export_fig(h, results_file, '-pdf'); %% how come this does not export most plots ??
end


if ismember(6,figstoplot)
    h = figure('units','normalized','position',[0 0 1 1]);
    y_lim={[-15 15] [-3 3] [-3 5] [-4 3]};
    for m=1:numel(monkeys)
        M=monkeys{m};
        tartmp=D.LFP.(M).tar;
        %curr_ylim=y_lim{(m-1)*2+y};
        for t=1:numel(targets)
            T=targets{t};
            [~,Ti]=ismember({tartmp.target},T);
            for c=1:2
                %% POWER
                sph(t,m,c,1)=subplot(6,4,(t-1)*8+(m-1)*4+(c-1)*2 +1);
                cname=tartmp(Ti==1).con(c).cond_name;
                if m==1 && t==1
                    title(['POW-' cname])
                end
                if c==1
                    ylabel({['Monkey ' M],[T, ',' num2str(tartmp(Ti==1).nSites) ' sites'],'Frequency (Hz)'});
                else
                    ylabel('Frequency (Hz)');
                end
                if m==numel(monkeys) && t==numel(targets)
                    xlabel('Time relative to R peak (s)');
                end
                con=tartmp(Ti==1).con(c).(E);
                toplot=con.pow_poppval';
                
                %                toplot=abs(toplot);
                %                 toplot(toplot>0.1)=0.1;
                %                 toplot=1-toplot;
                
                sigplot=con.pow_popsig;
                
                ax=sph(t,m,c,1);
                %xlim([min(lfp_time) max(lfp_time)]);
                xlim([0 size(toplot,2)]);
                axis square
                hold on
                
                %image(lfp_time-sl/2,1:size(toplot,1),toplot,'CDataMapping','scaled');
                image(toplot,'CDataMapping','scaled');
                set(gca,'YDir','normal');
                significance = double(sigplot);
                %                 contour(lfp_time-sl/2,(1:size(toplot,1)),significance==-1,1,'linecolor','b')
                %                 contour(lfp_time-sl/2,(1:size(toplot,1)),significance==1,1,'linecolor','r')
                %
                ecg_bna_draw_outlines(significance==1,'r')
                ecg_bna_draw_outlines(significance==-1,'b')
                
                % horizontal lines to separate frequency bands
                fbandstart = unique(frequency_bands(:))';
                fbandstart_idx = zeros(size(fbandstart));
                for f = fbandstart
                    f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                    line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                    fbandstart_idx(fbandstart == f) = f_idx;
                end
                
                line([0 0], ylim, 'color', 'k');
                set(ax,'TickDir','out')
                set(ax, 'ytick', fbandstart_idx);
                set(ax, 'yticklabel', fbandstart);
                set(ax, 'xtick', ticksx);
                set(ax, 'xticklabel', ticklabelsx);
                set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                box on
                
                %% POWER sig
                sph(t,m,c,2)=subplot(6,4,(t-1)*8+(m-1)*4+(c-1)*2 +2);
                cname=tartmp(Ti==1).con(c).cond_name;
                if m==1 && t==1
                    title(['Sig. POW-' cname])
                end
                if m==numel(monkeys) && t==numel(targets)
                    xlabel('Time relative to R peak (s)');
                end
                toplot=abs(con.powbp_sig);
                
                ax=sph(t,m,c,2);
                xlim([min(lfp_time) max(lfp_time)]);
                axis square
                hold on;
                set(ax,'ColorOrder',colsbp);
                smoothed=smoothit(toplot);
                plot(lfp_time,smoothed')
                if m==1 && t==1 && c==1
                    legend(freqb,'fontsize',3);
                end
            end
        end
    end
    
    %% adjusting limits
    for m=1:numel(monkeys)
        for t=1:numel(targets)
            y_lim=[0,max([get(sph(t,m,1,2),'ylim') get(sph(t,m,2,2),'ylim')])];
            c_lim=max(abs([get(sph(t,m,1,1),'clim'), get(sph(t,m,2,1),'clim')]));
            for c=1:2
                % sig %
                ax=sph(t,m,c,2);
                set(ax,'ylim',y_lim);
                subplot(ax);
                line([0 0], y_lim, 'color', 'k');
                % POW
                ax=sph(t,m,c,1);
                subplot(ax);
                cm=hot(128);
                %cm=bluered_colormap;
                colormap(ax,cm);
                cb = colorbar;
                set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                %                 set(get(cb,'title'),'string', '1-p', 'fontsize',8);
                %                 set(ax,'clim',[0.9 1]);
            end
        end
    end
    
    results_file=['Y:\Projects\Pulv_bodysignal\Figures\Fig6-pvals_POWER_' withunits];
    wanted_size=[50 30];
    set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    %
    export_fig(h, results_file, '-pdf'); %% how come this does not export most plots ??
end


if ismember(7,figstoplot)
    h = figure('units','normalized','position',[0 0 1 1]);
    y_lim={[-15 15] [-3 3] [-3 5] [-4 3]};
    for m=1:numel(monkeys)
        M=monkeys{m};
        tartmp=D.LFP.(M).tar;
        %curr_ylim=y_lim{(m-1)*2+y};
        for t=1:numel(targets)
            T=targets{t};
            [~,Ti]=ismember({tartmp.target},T);
            for c=1:2
                %% ITPC
                sph(t,m,c,1)=subplot(6,4,(t-1)*8+(m-1)*4+(c-1)*2 +1);
                ax=sph(t,m,c,1);
                cname=tartmp(Ti==1).con(c).cond_name;
                if m==1 && t==1
                    title(['ITPC-' cname])
                end
                if c==1
                    ylabel({['Monkey ' M],[T, ',' num2str(tartmp(Ti==1).nSites) ' sites'],'Frequency (Hz)'});
                else
                    ylabel('Frequency (Hz)');
                end
                if m==numel(monkeys) && t==numel(targets)
                    xlabel('Time (s)');
                end
                con=tartmp(Ti==1).con(c).(E);
                toplot=fig5_log_plot_scale(con.itpc);
                
                
                sigplot=con.itpc_popsig;
                
                xlim([0 size(toplot,2)]);
                
                %xlim([min(lfp_time) max(lfp_time)]);
                axis square
                hold on
                
                %image(lfp_time-sl/2,1:size(toplot,1),toplot,'CDataMapping','scaled');
                image(toplot,'CDataMapping','scaled');
                set(gca,'YDir','normal');
                significance = double(sigplot);
                %                 contour(lfp_time-sl/2,(1:size(toplot,1)),significance==1,1,'linecolor','r')
                %                 contour(lfp_time-sl/2,(1:size(toplot,1)),significance==-1,1,'linecolor','b')
                
                
                ecg_bna_draw_outlines(significance==1,'r')
                ecg_bna_draw_outlines(significance==-1,'b')
                
                % horizontal lines to separate frequency bands
                fbandstart = unique(frequency_bands(:))';
                fbandstart_idx = zeros(size(fbandstart));
                for f = fbandstart
                    f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                    line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                    fbandstart_idx(fbandstart == f) = f_idx;
                end
                set(sph(t,m,1,1),'clim',[min(toplot(:)), max(toplot(:))]);
                line([ticksx(ticklabelsx==0) ticksx(ticklabelsx==0)], ylim, 'color', 'k');
                set(ax,'TickDir','out')
                set(ax, 'ytick', fbandstart_idx);
                set(ax, 'yticklabel', fbandstart);
                set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                set(ax, 'xtick', ticksx);
                set(ax, 'xticklabel', ticklabelsx);
                box on
                
                %% ITPC sig
                sph(t,m,c,2)=subplot(6,4,(t-1)*8+(m-1)*4+(c-1)*2 +2);
                ax=sph(t,m,c,2);
                hold on;
                cname=tartmp(Ti==1).con(c).cond_name;
                if m==1 && t==1
                    title(['Sig. ITPC-' cname])
                end
                if m==numel(monkeys) && t==numel(targets)
                    xlabel('Time relative to R peak (s)');
                end
                ylabel('Proportion of sites (%)');
                toplot=abs(con.itpcbp_sig);
                toplot=smoothit(toplot);
                y_lim=[0 ceil(max(toplot(:))/10)*10];
                set(ax,'ylim',y_lim);
                for bpf=1:size(toplot,1)
                    col=colsbp(bpf,:);
                    Sig=toplot(bpf,:);
                    [pkv,pki]=max(Sig(lfp_time>=0));
                    pkt=lfp_time(pki+sum(lfp_time<0));
                    sepy=diff(y_lim)/20;
                    if pkv>(y_lim(end)-sepy)
                        plot(pkt,pkv-sepy,'color',col,'Marker','^');
                        text(pkt-0.02,pkv-sepy,[num2str(round(pkt*100)*10) 'ms'],'color',col,'HorizontalAlignment', 'right');
                    else
                        plot(pkt,pkv+sepy,'color',col,'Marker','v');
                        text(pkt-0.02,pkv+sepy,[num2str(round(pkt*100)*10) 'ms'],'color',col,'HorizontalAlignment', 'right');
                    end
                end
                
                xlim([min(lfp_time) max(lfp_time)]);
                axis square
                set(ax,'ColorOrder',colsbp);
                plot(lfp_time,toplot')
                if m==1 && t==1 && c==1
                    legend(freqb,'fontsize',3);
                end
            end
        end
    end
    
    %% adjusting limits
    for m=1:numel(monkeys)
        for t=1:numel(targets)
            y_lim=[0,max([get(sph(t,m,1,2),'ylim') get(sph(t,m,2,2),'ylim')])];
            c_lim=max([get(sph(t,m,1,1),'clim'), get(sph(t,m,2,1),'clim')]);
            for c=1:2
                % sig %
                ax=sph(t,m,c,2);
                set(ax,'ylim',y_lim);
                subplot(ax);
                line([0 0], y_lim, 'color', 'k');
                % itpc
                ax=sph(t,m,c,1);
                subplot(ax);
                
                apply_fig5_itpc_colorbar(ax, c_lim, 'n. ITPC');
            end
        end
    end
    
    results_file=['Y:\Projects\Pulv_bodysignal\Figures\Fig7-ITPC_' withunits];
    wanted_size=[50 30];
    set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    %
    export_fig(h, results_file, '-pdf'); %% how come this does not export most plots ??
end


if ismember(7.5,figstoplot)
    h = figure('units','normalized','position',[0 0 1 1]);
    y_lim={[-15 15] [-3 3] [-3 5] [-4 3]};
    for m=1:numel(monkeys)
        M=monkeys{m};
        tartmp=D.LFP.(M).tar;
        %curr_ylim=y_lim{(m-1)*2+y};
        for t=1:numel(targets)
            T=targets{t};
            [~,Ti]=ismember({tartmp.target},T);
            for c=1:2
                %% ITPC
                sph(t,m,c,1)=subplot(6,4,(t-1)*8+(m-1)*4+(c-1)*2 +1);
                ax=sph(t,m,c,1);
                cname=tartmp(Ti==1).con(c).cond_name;
                if m==1 && t==1
                    title(['ITPC-' cname])
                end
                if c==1
                    ylabel({['Monkey ' M],[T, ',' num2str(tartmp(Ti==1).nSites) ' sites'],'Frequency (Hz)'});
                else
                    ylabel('Frequency (Hz)');
                end
                if m==numel(monkeys) && t==numel(targets)
                    xlabel('Time (s)');
                end
                con=tartmp(Ti==1).con(c).(E);
                toplot=fig5_log_plot_scale(con.itpc);
                
                
                sigplot=con.itpc_popsig;
                
                xlim([0 size(toplot,2)]);
                
                %xlim([min(lfp_time) max(lfp_time)]);
                axis square
                hold on
                
                %image(lfp_time-sl/2,1:size(toplot,1),toplot,'CDataMapping','scaled');
                image(toplot,'CDataMapping','scaled');
                set(gca,'YDir','normal');
                significance = double(sigplot);
                %                 contour(lfp_time-sl/2,(1:size(toplot,1)),significance==1,1,'linecolor','r')
                %                 contour(lfp_time-sl/2,(1:size(toplot,1)),significance==-1,1,'linecolor','b')
                
                
                ecg_bna_draw_outlines(significance==1,'r')
                ecg_bna_draw_outlines(significance==-1,'b')
                
                % horizontal lines to separate frequency bands
                fbandstart = unique(frequency_bands(:))';
                fbandstart_idx = zeros(size(fbandstart));
                for f = fbandstart
                    f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                    line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                    fbandstart_idx(fbandstart == f) = f_idx;
                end
                set(sph(t,m,1,1),'clim',[min(toplot(:)), max(toplot(:))]);
                line([ticksx(ticklabelsx==0) ticksx(ticklabelsx==0)], ylim, 'color', 'k');
                set(ax,'TickDir','out')
                set(ax, 'ytick', fbandstart_idx);
                set(ax, 'yticklabel', fbandstart);
                set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                set(ax, 'xtick', ticksx);
                set(ax, 'xticklabel', ticklabelsx);
                box on
                
                %% ITPC sig
                sph(t,m,c,2)=subplot(6,4,(t-1)*8+(m-1)*4+(c-1)*2 +2);
                ax=sph(t,m,c,2);
                hold on;
                cname=tartmp(Ti==1).con(c).cond_name;
                if m==1 && t==1
                    title(['Sig. ITPC-' cname])
                end
                if m==numel(monkeys) && t==numel(targets)
                    xlabel('Time relative to R peak (s)');
                end
                ylabel('Proportion of sites (%)');
                
                
                toplot1=con.itpcbp_sig_pos;
                toplot1=smoothit(toplot1);
                toplot2=con.itpcbp_sig_neg;
                toplot2=smoothit(toplot2);
                ylabel('Proportion of sites (%)');
                
                y_lim=[floor(min(toplot2(:))/10)*10 ceil(max(toplot1(:))/10)*10];
                set(ax,'ylim',y_lim);
                for bpf=1:size(toplot1,1)
                    col=colsbp(bpf,:);
                    Sig=toplot1(bpf,:);
                    [pkv,pki]=max(Sig(lfp_time>=0));
                    pkt=lfp_time(pki+sum(lfp_time<0));
                    sepy=diff(y_lim)/20;
                    if pkv>(y_lim(end)-sepy)
                        plot(pkt,pkv-sepy,'color',col,'Marker','^');
                        text(pkt-0.02,pkv-sepy,[num2str(round(pkt*100)*10) 'ms'],'color',col,'HorizontalAlignment', 'right');
                    else
                        plot(pkt,pkv+sepy,'color',col,'Marker','v');
                        text(pkt-0.02,pkv+sepy,[num2str(round(pkt*100)*10) 'ms'],'color',col,'HorizontalAlignment', 'right');
                    end
                end
                
                xlim([min(lfp_time) max(lfp_time)]);
                axis square
                set(ax,'ColorOrder',colsbp);
                plot(lfp_time,toplot1')
                plot(lfp_time,toplot2')
                if m==1 && t==1 && c==1
                    legend(freqb,'fontsize',3);
                end
            end
        end
    end
    
    %% adjusting limits
    for m=1:numel(monkeys)
        for t=1:numel(targets)
            y_lim=[min([get(sph(t,m,1,2),'ylim') get(sph(t,m,2,2),'ylim')]),max([get(sph(t,m,1,2),'ylim') get(sph(t,m,2,2),'ylim')])];
            c_lim=max([get(sph(t,m,1,1),'clim'), get(sph(t,m,2,1),'clim')]);
            for c=1:2
                % sig %
                ax=sph(t,m,c,2);
                set(ax,'ylim',y_lim);
                subplot(ax);
                line([0 0], y_lim, 'color', 'k');
                % itpc
                ax=sph(t,m,c,1);
                subplot(ax);
                
                apply_fig5_itpc_colorbar(ax, c_lim, 'n. ITPC');
            end
        end
    end
    
    results_file=['Y:\Projects\Pulv_bodysignal\Figures\Fig7,5-ITPC_' withunits];
    wanted_size=[50 30];
    set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    %
    export_fig(h, results_file, '-pdf'); %% how come this does not export most plots ??
end






if ismember(7,figstoplot)
    h = figure('units','normalized','position',[0 0 1 1]);
    y_lim={[-15 15] [-3 3] [-3 5] [-4 3]};
    for m=1:numel(monkeys)
        M=monkeys{m};
        tartmp=D.LFP.(M).tar;
        %curr_ylim=y_lim{(m-1)*2+y};
        for t=1:numel(targets)
            T=targets{t};
            [~,Ti]=ismember({tartmp.target},T);
            for c=1:2
                %% ITPC
                sph(t,m,c,1)=subplot(6,4,(t-1)*8+(m-1)*4+(c-1)*2 +1);
                cname=tartmp(Ti==1).con(c).cond_name;
                if m==1 && t==1
                    title(['ITPC-' cname])
                end
                if c==1
                    ylabel({['Monkey ' M],[T, ',' num2str(tartmp(Ti==1).nSites) ' sites'],'Frequency (Hz)'});
                else
                    ylabel('Frequency (Hz)');
                end
                if m==numel(monkeys) && t==numel(targets)
                    xlabel('Time (s)');
                end
                con=tartmp(Ti==1).con(c).(E);
                toplot=con.itpc_poppval';
                
                %                toplot=abs(toplot);
                %                 toplot(toplot>0.1)=0.1;
                %                 toplot=1-toplot;
                
                sigplot=con.itpc_popsig;
                
                ax=sph(t,m,c,1);
                %xlim([min(lfp_time) max(lfp_time)]);
                xlim([0 size(toplot,2)]);
                axis square
                hold on
                
                %image(lfp_time-sl/2,1:size(toplot,1),toplot,'CDataMapping','scaled');
                image(toplot,'CDataMapping','scaled');
                set(gca,'YDir','normal');
                significance = double(sigplot);
                
                ecg_bna_draw_outlines(significance==1,'r')
                ecg_bna_draw_outlines(significance==-1,'b')
                %                 contour(lfp_time-sl/2,(1:size(toplot,1)),significance==-1,1,'linecolor','b')
                %                 contour(lfp_time-sl/2,(1:size(toplot,1)),significance==1,1,'linecolor','r')
                %
                % horizontal lines to separate frequency bands
                fbandstart = unique(frequency_bands(:))';
                fbandstart_idx = zeros(size(fbandstart));
                for f = fbandstart
                    f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                    line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                    fbandstart_idx(fbandstart == f) = f_idx;
                end
                
                line([0 0], ylim, 'color', 'k');
                set(ax,'TickDir','out')
                set(ax, 'ytick', fbandstart_idx);
                set(ax, 'yticklabel', fbandstart);
                set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                set(ax, 'xtick', ticksx);
                set(ax, 'xticklabel', ticklabelsx);
                box on
                
                %% ITPC sig
                sph(t,m,c,2)=subplot(6,4,(t-1)*8+(m-1)*4+(c-1)*2 +2);
                cname=tartmp(Ti==1).con(c).cond_name;
                if m==1 && t==1
                    title(['Sig. ITPC-' cname])
                end
                if m==numel(monkeys) && t==numel(targets)
                    xlabel('Time (s)');
                end
                toplot=abs(con.itpcbp_sig);
                
                ax=sph(t,m,c,2);
                xlim([min(lfp_time) max(lfp_time)]);
                axis square
                hold on;
                set(ax,'ColorOrder',colsbp);
                smoothed=smoothit(toplot);
                plot(lfp_time,smoothed')
                if m==1 && t==1 && c==1
                    legend(freqb,'fontsize',3);
                end
            end
        end
    end
    
    %% adjusting limits
    for m=1:numel(monkeys)
        for t=1:numel(targets)
            y_lim=[0,max([get(sph(t,m,1,2),'ylim') get(sph(t,m,2,2),'ylim')])];
            c_lim=max([get(sph(t,m,1,1),'clim'), get(sph(t,m,2,1),'clim')]);
            for c=1:2
                % sig %
                ax=sph(t,m,c,2);
                set(ax,'ylim',y_lim);
                subplot(ax);
                line([0 0], y_lim, 'color', 'k');
                % ITPC
                ax=sph(t,m,c,1);
                subplot(ax);
                cm=hot(128);
                colormap(ax,cm);
                cb = colorbar;
                set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                
                %
                %                 set(get(cb,'title'),'string', '1-p', 'fontsize',8);
                %                 set(ax,'clim',[0.9 1]);
            end
        end
    end
    
    results_file=['Y:\Projects\Pulv_bodysignal\Figures\Fig7-pvals_ITPC_' withunits];
    wanted_size=[50 30];
    set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    %
    export_fig(h, results_file, '-pdf'); %% how come this does not export most plots ??
end


if ismember(8,figstoplot)
    
    sr=1017.2526041; %hz
    sl=1/sr*10;
    lfp_time=[-sl*11:sl:sl*41];
    ticklabelsx=[0,0.2,0.4];
    for tl=1:numel(ticklabelsx)
        tt=ticklabelsx(tl);
        tlstepsize=numel(lfp_time)/(lfp_time(end)-lfp_time(1));
        ticksx(tl)=tt*tlstepsize-lfp_time(1)*tlstepsize;
    end
    
    
    h = figure('units','normalized','position',[0 0 1 1]);
    clear sph
    E='Cue';
    c=2; %% task?
    for m=1:numel(monkeys)
        M=monkeys{m};
        tartmp=D.LFP.(M).tar;
        for t=1:numel(targets)
            T=targets{t};
            [~,Ti]=ismember({tartmp.target},T);
            
            %% POWER
            sph(t,m,1)=subplot(6,4,(t-1)*8+(m-1)*4 +1);
            ax=sph(t,m,1);
            cname=tartmp(Ti==1).con(c).cond_name;
            if m==1 && t==1
                title('POW-Cue')
            end
            ylabel({['Monkey ' M],[T, ',' num2str(tartmp(Ti==1).nSites) ' sites'],'Frequency (Hz)'});
            if m==numel(monkeys) && t==numel(targets)
                xlabel('Time relative to R peak (s)');
            end
            con=tartmp(Ti==1).con(c).(E);
            toplot=con.pow;
            sigplot=con.pow_popsig;
            
            xlim([0 size(toplot,2)]);
            axis square
            hold on
            
            image(toplot,'CDataMapping','scaled');
            set(gca,'YDir','normal');
            significance = double(sigplot);
            ecg_bna_draw_outlines(significance==1,'r')
            ecg_bna_draw_outlines(significance==-1,'b')
            
            % horizontal lines to separate frequency bands
            fbandstart = unique(frequency_bands(:))';
            fbandstart_idx = zeros(size(fbandstart));
            for f = fbandstart
                f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                fbandstart_idx(fbandstart == f) = f_idx;
            end
            line([ticksx(ticklabelsx==0) ticksx(ticklabelsx==0)], ylim, 'color', 'k');
            
            set(ax, 'TickDir','out')
            set(ax, 'ytick', fbandstart_idx);
            set(ax, 'yticklabel', fbandstart);
            set(ax, 'xtick', ticksx);
            set(ax, 'xticklabel', ticklabelsx);
            set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
            box on
            
            
            %% ITPC
            sph(t,m,2)=subplot(6,4,(t-1)*8+(m-1)*4 +2);
            ax=sph(t,m,2);
            if m==1 && t==1
                title('ITPC-Cue')
            end
            ylabel('Frequency (Hz)');
            if m==numel(monkeys) && t==numel(targets)
                xlabel('Time relative to R peak (s)');
            end
            con=tartmp(Ti==1).con(c).(E);
            toplot=fig5_log_plot_scale(con.itpc);
            sigplot=con.itpc_popsig;
            xlim([0 size(toplot,2)]);
            axis square
            hold on
            
            image(toplot,'CDataMapping','scaled');
            set(gca,'YDir','normal');
            significance = double(sigplot);
            ecg_bna_draw_outlines(significance==1,'r')
            ecg_bna_draw_outlines(significance==-1,'b')
            
            % horizontal lines to separate frequency bands
            fbandstart = unique(frequency_bands(:))';
            fbandstart_idx = zeros(size(fbandstart));
            for f = fbandstart
                f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                fbandstart_idx(fbandstart == f) = f_idx;
            end
            line([ticksx(ticklabelsx==0) ticksx(ticklabelsx==0)], ylim, 'color', 'k');
            
            set(ax, 'TickDir','out')
            set(ax, 'ytick', fbandstart_idx);
            set(ax, 'yticklabel', fbandstart);
            set(ax, 'xtick', ticksx);
            set(ax, 'xticklabel', ticklabelsx);
            
            
            set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
            box on
            
            %% evoked LFP
            sph(t,m,3)=subplot(6,4,(t-1)*8+(m-1)*4 +3);
            ax=sph(t,m,3);
            if m==1 && t==1
                title('LFP-Cue')
            end
            ylabel('LFP');
            if m==numel(monkeys) && t==numel(targets)
                xlabel('Time relative to R peak (s)');
            end
            con=D.LFP.(M).tar(Ti==1).con(c).(E);
            
            xlim([min(lfp_time) max(lfp_time)]);
            
            axis square
            hold on
            
            col=cfg.condition(c).color;
            mean_tp=  smoothit(con.evoked);
            EB_tp=  smoothit(con.evoked_sterr);
            lineProps={'color',col};
            shadedErrorBar(lfp_time,mean_tp,EB_tp,lineProps,1);
            box on
            
            %% evoked MUA
            sph(t,m,4)=subplot(6,4,(t-1)*8+(m-1)*4 +4);
            ax=sph(t,m,4);
            if m==1 && t==1
                title('MUA-Cue')
            end
            ylabel('MUA');
            if m==numel(monkeys) && t==numel(targets)
                xlabel('Time relative to R peak (s)');
            end
            con=D.MUA.(M).tar(Ti==1).con(c).(E);
            xlim([min(lfp_time) max(lfp_time)]);
            axis square
            hold on
            
            col=cfg.condition(c).color;
            mean_tp=  smoothit(con.evoked);
            EB_tp=  smoothit(con.evoked_sterr);
            lineProps={'color',col};
            shadedErrorBar(lfp_time,mean_tp,EB_tp,lineProps,1);
            box on
        end
    end
    
    
    
    %% adjusting limits
    for m=1:numel(monkeys)
        for t=1:numel(targets)
            % LFP & MUA
            for s=3:4
                ax=sph(t,m,s);
                subplot(ax);
                y_lim=get(ax,'ylim');
                plot(ax,[0 0], y_lim,'color','k');
            end
            
            c_lim_pow=max(abs([get(sph(t,m,1),'clim'), get(sph(t,m,1),'clim')]));
            c_lim_itpc=max([get(sph(t,m,2),'clim')]);
            
            % POW
            ax=sph(t,m,1);
            subplot(ax);
            cm=customcolormap_preset('red-yellow-blue',256);
            colormap(ax,cm);
            cb = colorbar;
            set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
            set(get(cb,'title'),'string', 'n. Pow', 'fontsize',8);
            set(ax,'clim',[-c_lim_pow c_lim_pow]);
            
            
            
            % ITPC
            ax=sph(t,m,2);
            subplot(ax);
            apply_fig5_itpc_colorbar(ax, c_lim_itpc, 'n. ITPC');
        end
    end
    
    results_file=['Y:\Projects\Pulv_bodysignal\Figures\Fig8-CUE_' withunits];
    wanted_size=[50 30];
    set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    %
    export_fig(h, results_file, '-pdf'); %% how come this does not export most plots ??
end

close all
end

function out=smoothit(in)

cfg.lfp.smoothWin = 5;
win = 1:cfg.lfp.smoothWin;
win=win-(numel(win)+1)/2;
half_win = ceil(size(win,2)/2)-1;
gaussian_kernel=normpdf(win,0,numel(win)/6);
gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);



jnk=[];
concat_input = cat(2,(in(:,half_win:-1:1)),in);
concat_input = cat(2,concat_input, (in(:,end:-1:end-half_win+1)));
for k=1:size(concat_input,1)
    jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');
end
out = jnk(:,half_win+1:end-half_win);
end

function y = fig5_log_plot_scale(x)
% |x| <= 1: unchanged; |x| > 1: sign(x) * (log2(|x|) + 1)
y = x;
pos = x > 1;
neg = x < -1;
y(pos) = log2(x(pos)) + 1;
y(neg) = -(log2(abs(x(neg))) + 1);
end

function [tick_pos, tick_labels] = fig5_itpc_colorbar_ticks(c_lim)
% ITPC ticks: negative side clipped at -1 in plot space; labels = pre-log integers.
orig_ticks = [-1, 0, 1];
if c_lim > 1
    max_power = floor(c_lim - 1);
    orig_ticks = [orig_ticks, 2.^(1:max_power)];
end
tick_pos = fig5_log_plot_scale(orig_ticks);
valid = tick_pos >= -1 - 1e-6 & tick_pos <= c_lim + 1e-6;
tick_pos = tick_pos(valid);
orig_ticks = orig_ticks(valid);
[tick_pos, sort_ix] = sort(tick_pos);
tick_labels = arrayfun(@num2str, orig_ticks(sort_ix), 'UniformOutput', false);
end

function apply_fig5_itpc_colorbar(ax, c_lim, cb_title)
offset = round(128/c_lim*(1));
cm = customcolormap_preset('red-yellow-blue',256);
cm = cm((128-offset):end,:);
colormap(ax, cm);
cb = colorbar;
[tick_pos, tick_labels] = fig5_itpc_colorbar_ticks(c_lim);
set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
set(get(cb,'title'),'string', cb_title, 'fontsize', 8);
set(ax,'clim',[-1, c_lim]);
cb.Ticks = tick_pos;
cb.TickLabels = tick_labels;
end


%%

% %% time axis management
% tmp_event=LFP.M.tar(1).con(c).(E);
% tfr_time=tmp_event.tfr_time;
% lfp_time=tmp_event.time;
% onset_s = find(tfr_time <= 0, 1, 'last'); % state onset time
% finish_s = length(tfr_time); % start end sample
% ticksamples_tfr=[1 onset_s finish_s];
% onset_s = find(lfp_time <= 0, 1, 'last'); % state onset time
% finish_s = length(lfp_time); % start end sample
% ticksamples_lfp=[1 onset_s finish_s];
%
% %% plot
% tfr_events.onset        = find(tmp_event.tfr_time == 0);
% tfr_events.name         ={tmp_event.event_name};
% tfr_events.ticksamples  = sort(ticksamples_tfr);
% tfr_events.startsamples = ticksamples_tfr(1:3:end);
% tfr_events.endsamples   = ticksamples_tfr(3:3:end);
% tfr_events.ticks        =floor(tmp_event.tfr_time(tfr_events.ticksamples)*100)/100;
%
% lfp_events.onset        = find(tmp_event.time == 0);
% lfp_events.name         ={tmp_event.event_name};
% lfp_events.ticksamples  = sort(ticksamples_lfp);
% lfp_events.startsamples = ticksamples_lfp(1:3:end);
% lfp_events.endsamples   = ticksamples_lfp(3:3:end);
% lfp_events.ticks        =floor(tmp_event.time(lfp_events.ticksamples)*100)/100;

%% frequencies !

% potential cutoff -> probably startfreq should go into settings at somepoint
%
% startfreq=3.5;
% %colsbp=jet(length(cfg.lfp.frequency_bands));
% fbx=find(cfg.lfp.frequency_bands(:,1)>=startfreq,1,'first');
% frequency_bands=cfg.lfp.frequency_bands(fbx:end,:);
% %freqName = cfg.lfp.freqName(fbx:end);
% %colsbp=colsbp(fbx:end,:);
% %freqb=num2cell(strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz'),2);
% fx=find(cfg.lfp.foi>=startfreq,1,'first');
% freq = cfg.lfp.foi(fx:end);
%
% methods= {'real','shuffled','normalized'};
% cols= {[0 0 1],[1 0 0],[0 1 0],};