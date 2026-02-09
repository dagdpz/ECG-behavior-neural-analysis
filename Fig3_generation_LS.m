function Fig3_generation_LS
clc
%%

%sessions = {'20210720','20210730', '20220315'};

targets = {'VPL','dPul', 'MD'};
% project = 'Pulv_bodysignal';
% version = 'ECG_Bacchus_unfiltered';
% driver_path = 'Y:';%'/home/shamim/fileserver';
%showfiltered='filtered';

% cfg = [];
% cfg.ecg.timestep = 1; %??
%
% cfg.process_ecg=0;
% cfg.project = project;
% cfg.results_folder = [driver_path,filesep,'Projects',filesep,cfg.project];
% ecg_bna_location     =which('ecg_bna_define_folders');
% github_folder        =ecg_bna_location(1:strfind(ecg_bna_location,['ECG-behavior-neural-analysis' filesep 'ecg_bna_define_folders'])-1);
% cfg.version = version;
% run([github_folder filesep 'Settings' filesep cfg.project filesep 'ECG_bna' filesep cfg.version '.m']);
% cfg = ecg_bna_define_folders(cfg);

monkeys = {'Bac','Bac', 'Bac'};
monkey = 'Bacchus';
sessions = {'20210720','20210730', '20220211'};
version = 'ECG_Bacchus_complete';


% Site_ID = 'Bac_20220315_Site_15';%'Bac_20211007_Site_01';%'Bac_20210730_Site_03';
% filteredornot={'unfiltered','filtered'};
%for flt=1:numel(filteredornot)
%showfiltered='filtered';




sit=0;
for s=1:numel(sessions)
    ses=sessions{s};
    base_path=['Y:\Projects\Pulv_bodysignal\LFP\' version];
    site_files= dir([base_path '\Per_Site\*' ses '*.mat']);
    
    sit_ses=sit;
    for sf=1:numel(site_files)
        Site_ID=site_files(sf).name(1:end-4);
        load([base_path '\Per_Site\' Site_ID '.mat'])
        isvalid= ~isempty(triggered_site_data.condition(1).event) &&...
            ~isempty(triggered_site_data.condition(2).event) &&...
            ~isempty(strfind(triggered_site_data.target,targets{s}));
        if isvalid
            sit=sit+1;
            AllLFP(sit)=triggered_site_data;
            AllLFP(sit)=ecg_bna_add_sterr_posthoc(AllLFP(sit));
            %base_path='Y:\Projects\Pulv_bodysignal\MUA\ECG_Bacchus_TaskRest_generalTrig_MUA\';
            load(['Y:\Projects\Pulv_bodysignal\MUA\' version '\Per_Site\' Site_ID '.mat'])
            AllMUA(sit)=triggered_site_data;
            AllMUA(sit)=ecg_bna_add_sterr_posthoc(AllMUA(sit));
        end
    end
    sitespersession(s)=sit-sit_ses;
end

%% find depths

% reading all the electrode depth and store it perSite:
if exist(['Y:\Projects\Pulv_bodysignal\MUA\ECG_',monkey,'_TaskRest_generalTrig_MUA',filesep,'allSites_info.mat'],'file')
    load(['Y:\Projects\Pulv_bodysignal\MUA\ECG_',monkey,'_TaskRest_generalTrig_MUA',filesep,'allSites_info.mat']);
else
    % this part not so clear to me
    cd(ephys_path)
    all_perSite_raw_data = dir(['sites_',monkey,'_*.mat']);
    for si = 1:length(all_perSite_raw_data)
        load(all_perSite_raw_data(si).name)
        siteFields = fieldnames(sites);
        site_info(si) = rmfield(sites,{siteFields{8:end-2}}); %% remove everuything exepth required info
    end
    clear sites;
    save(fullfile(['Y:\Projects\Pulv_bodysignal\MUA\ECG_',monkey,'_TaskRest_generalTrig_MUA',filesep,'allSites_info.mat']),'site_info');
end

%% add x/y/z to sites:
for s=1:numel(AllLFP)
    site_ID=AllLFP(s).site_ID;
    ix=ismember({site_info.site_ID},site_ID);
    AllLFP(s).x=site_info(ix).grid_x*0.7;
    AllLFP(s).y=site_info(ix).grid_y*0.7;
    AllLFP(s).z=round(site_info(ix).electrode_depth*10)/10;
end

nrows=max(sitespersession);
ncolumns=numel(targets)*2;
e=1;
M='normalized';
%c: to loop through!
%% time
sr=1017; %hz
sl=1/sr*10;
time=[-sl*26:sl:sl*26];


%% start plotting
h = figure('units','normalized','position',[0 0 1 1]);
cols= {[0 0 1],[1 0 0]};

for ss=1:numel(sessions)
    cs=ismember({AllLFP.session},[monkeys{ss} '_' sessions{ss}]);
    LFP=AllLFP(cs);
    MUA=AllMUA(cs);
    [~,cs]=sort([LFP.z]);
    LFP=LFP(cs);
    MUA=MUA(cs);
    for s=1:numel(LFP)
        
        %% LFP evoked
        sph(s)=subplot(ncolumns, nrows, (ss-1)*2+(s-1)*ncolumns+1);
        hold on;
        if s==1
           us_loc=strfind(LFP(s).target,'_');
           title(LFP(s).target(1:us_loc-1));
        end
        line([time(1) time(end)], [0 0], 'color', 'k');
        ymax=[];
        ymin=[];
        for c=1:2
            col=cols{c};
            TP=LFP(s).condition(c).event(e);
            TPG=TP.(M);
            conf=smoothit(squeeze(TPG.evoked.sterr)');
            
            % error bar
            lineprops={'color',col};
            meantp=smoothit(squeeze(TPG.evoked.mean)');
            shadedErrorBar(time,meantp,conf,lineprops,1);
            ymin=min([ymin meantp-conf]);
            ymax=max([ymax meantp+conf]);
        end
        
        ylm = [ymin ymax];
        line([0 0], ylm, 'color', 'k');
        ax=sph(s);
        %xyz=['x/y/z:' num2str(LFP(s).x) '/' num2str(LFP(s).y) '/' num2str(LFP(s).z)];
        xyz=['depth: ' num2str(LFP(s).z-22) 'mm'];
        if strcmp(M,'normalized')
            ylabel({xyz,'normalized LFP'});
        else
            ylabel({xyz,'Voltage (mV)'});
        end
        xlabel('Time(s)');
        %title(titles{mt},'Interpreter', 'none');
        set(ax, 'xlim', [time(1) time(end)]);
        set(ax, 'ylim', ylm);
        axis square
        
        for c=1:2
            col=cols{c};
            TP=LFP(s).condition(c).event(e);
            TPS=TP.significance;
            significance = double(squeeze(TPS.evoked));
            significance(significance==0)=NaN;
            if any(~isnan(significance))
                sigpos=significance;sigpos(sigpos==-1)=NaN;sigpos=sigpos.*diff(ylm)/40*c+ylm(1);
                signeg=significance;signeg(signeg==1)=NaN;signeg=abs(signeg).*diff(ylm)/40*c+ylm(1);
                
%             significance=significance.*diff(ylm)/40*c+ylm(1);
%             plot(ax,time,significance','linewidth',2,'color',col);
            plot(ax,time,sigpos','linewidth',2,'color',col);
            plot(ax,time,signeg','linewidth',2,'color',col/2);
            
            end
        end
        
        
        %% MUA
        sph(s)=subplot(ncolumns, nrows, (ss-1)*2+(s-1)*ncolumns+2);
        hold on;
        if s==1
           title(MUA(s).session,'interpreter','none');
        end
        line([time(1) time(end)], [0 0], 'color', 'k');
        ymax=[];
        ymin=[];
        for c=1:2
            col=cols{c};
            TP=MUA(s).condition(c).event(e);
            TPG=TP.(M);
            TPS=TP.significance;
            conf=smoothit(squeeze(TPG.evoked.sterr)');
            
            % error bar
            lineprops={'color',col};
            TPGM=TPG.evoked.mean;
            meantp=smoothit(squeeze(TPGM)');
            shadedErrorBar(time,meantp,conf,lineprops,1);
            ymin=min([ymin meantp-conf]);
            ymax=max([ymax meantp+conf]);
        end
        
        ylm = [ymin ymax];
        line([0 0], ylm, 'color', 'k');
            ax=sph(s);                    
        if strcmp(M,'normalized')
            ylabel('normalized MUA');
        else
            ylabel('Voltage (mV)');
        end
        xlabel('Time(s)');        
        set(ax, 'xlim', [time(1) time(end)]);
        set(ax, 'ylim', ylm);
        axis square
        
        for c=1:2
            col=cols{c};
            TP=MUA(s).condition(c).event(e);
            TPS=TP.significance;
            significance = double(squeeze(TPS.evoked));
            significance(significance==0)=NaN;
            if any(~isnan(significance))
                sigpos=significance;sigpos(sigpos==-1)=NaN;sigpos=sigpos.*diff(ylm)/40*c+ylm(1);
                signeg=significance;signeg(signeg==1)=NaN;signeg=abs(signeg).*diff(ylm)/40*c+ylm(1);
                
%             significance=significance.*diff(ylm)/40*c+ylm(1);
%             plot(ax,time,significance','linewidth',2,'color',col);
            plot(ax,time,sigpos','linewidth',2,'color',col);
            plot(ax,time,signeg','linewidth',2,'color',col/2);
            end
        end
        
    end
end

results_file=['Y:\Projects\Pulv_bodysignal\Figures\Fig3_close_electrodes'];
wanted_size=[50 30];
set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    %
export_fig(h, results_file, '-pdf'); 
end
% close all
% end


function out=smoothit(in)

cfg.lfp.smoothWin = 5;
win = 1:cfg.lfp.smoothWin;
win=win-(numel(win)+1)/2;
half_win = ceil(size(win,2)/2)-1;
gaussian_kernel=normpdf(win,0,1);
gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);



jnk=[];
concat_input = cat(2,(in(:,half_win:-1:1)),in);
concat_input = cat(2,concat_input, (in(:,end:-1:end-half_win+1)));
for k=1:size(concat_input,1)
    jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');
end
out = jnk(:,half_win+1:end-half_win);
end

