function ecg_bna_plotting_MUA_LFP_per_elec_depth

% ==================================================================
% just to plot the LFP and MUA per electrode depth
% in the squared axis all target in one plot >> 27.07.25
% =================================================================
close all,
clear all
clc

monkey = 'Bacchus';%'Magnus' or 'Bacchus'

Fin = {'lfp','mua'};
for fin = 1%1:length(Fin)
    withunits = 'all_sites';
    
    if strcmp(Fin{fin},'mua')
        path_to_save = ['Y:\Projects\Pulv_bodysignal\MUA\ECG_',monkey,'_TaskRest_generalTrig_MUA\Per_Session_per_depth'];
        perSitedata_path = ['Y:\Projects\Pulv_bodysignal\MUA\ECG_',monkey,'_TaskRest_generalTrig_MUA\Per_Site'];
        ephys_path = ['Y:\Projects\Pulv_bodysignal\ephys\ECG_TaskRest_',monkey,'_MUA'];
    elseif strcmp(Fin{fin},'lfp')
        path_to_save = ['Y:\Projects\Pulv_bodysignal\LFP\ECG_',monkey,'_TaskRest_generalTrig_finalRun\Per_Session_per_depth'];
        perSitedata_path = ['Y:\Projects\Pulv_bodysignal\LFP\ECG_',monkey,'_TaskRest_generalTrig_finalRun\4hz filtered\Per_Site'];
        ephys_path = ['Y:\Projects\Pulv_bodysignal\ephys\ECG_TaskRest_',monkey,'_MUA']; % it has both mua and LFP so its similar data
    end
    %%
    tic
    % reading all the electrode depth and store it perSite:
    if exist(['Y:\Projects\Pulv_bodysignal\MUA\ECG_',monkey,'_TaskRest_generalTrig_MUA',filesep,'allSites_info.mat'],'file')
        load(['Y:\Projects\Pulv_bodysignal\MUA\ECG_',monkey,'_TaskRest_generalTrig_MUA',filesep,'allSites_info.mat']);
    else
        
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
    toc
    %%
    targets = {'VPL', 'dPul', 'MD'};% unique({sites.target});
    cond = {'Rest','Task'};
    
    cfg.monkey = monkey;
    cfg.lfp.foi = logspace(log10(2),log10(120),60); %
    cfg.lfp.smoothWin = 5;
    cfg.lfp.frequency_bands = [2 4; 4 8; 8 14; 14 30; 30 50; 70 120];
    
    % Smoothing Kernel here:
    win = 1:cfg.lfp.smoothWin; win=win-(numel(win)+1)/2;
    half_win = ceil(size(win,2)/2)-1;
    gaussian_kernel=normpdf(win,0,numel(win)/6);
    gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);
    
    % cfg.analyse_states = {'R',    'Rpeak',1,-0.25, 0.25;...
    %                       'R_low', 'Rpeak_lowIBI',1,-0.25, 0.25;...
    %                       'R_high', 'Rpeak_highIBI',1,-0.25, 0.25;...
    %                       'Cue',  'state',4,-0.10, 0.4};
    cfg.analyse_states = {'R',    'Rpeak',1,-0.25, 0.25};
    %%
    tic
    cd(perSitedata_path)
    % read the files
    all_mua_data = dir('*.mat');
    % storing all the Triggered parameters of the sites in the results folder:
    sites = struct;
    s=0;
    for f = 1:length(all_mua_data)
        load(all_mua_data(f).name);
        session = triggered_site_data.session;
        site_ID = triggered_site_data.site_ID;
        switch withunits
            case 'w_units'
                if ~ismember(site_ID,cfg.site_IDS)
                    continue
                end
            case 'wo_units'
                if ismember(site_ID,cfg.site_IDS)
                    continue
                end
            case 'all_sites'
        end
        s=s+1;
        sites(s).all_conditions_present = 1;
        sites(s).site_ID = site_ID;
        sites(s).target = triggered_site_data.target;
        
        for c = 1:length(triggered_site_data.condition)
            con=triggered_site_data.condition(c);
            if isempty(con.event)
                sites(s).all_conditions_present = 0;
                continue;
            end
            sites(s).condition(c).condition_name = con.label;
            for e=1:size(cfg.analyse_states,1)
                event=con.event(e);
                %% think about event present?
                
                nTriggers = event.real.ntriggers;
                time      = event.time;
                eval([Fin{fin},'_time  = event.tfr_time;']);
                if strcmp(Fin{fin}, 'lfp')
                    eval([Fin{fin},'       = squeeze(event.real.',Fin{fin},'.mean);']); % not normalized LFP
                else
                    eval([Fin{fin},'       = squeeze(event.normalized.',Fin{fin},'.mean);']);
                end
                
                
                eval([Fin{fin},'_sig   =squeeze(event.significance.',Fin{fin},');']);
                
                sites(s).condition(c).event(e).nTriggers = nTriggers;
                eval(['sites(s).condition(c).event(e).',Fin{fin},' = ',Fin{fin},';']);
                eval(['sites(s).condition(c).event(e).',Fin{fin},'_sig = ',Fin{fin},'_sig;']);
                
                sites(s).condition(c).event(e).time = time;
                eval(['sites(s).condition(c).event(e).tfr_time = ',Fin{fin},'_time;']);
            end
        end
    end
    sites = sites([sites.all_conditions_present]==1);
    all_sites_ID = (ismember({site_info.site_ID},{sites.site_ID}));
    site_info = site_info(all_sites_ID);
    site_info = rmfield(site_info,'perturbation_site');
    for s = 1: length(sites)
        sites(s).site_info = site_info(s);
    end
    %     dum = (setdiff({site_info(all_sites_ID).site_ID},{sites.site_ID}))';
    %     % to check if there's andy order diff in the list of site_IDs
    
    toc
    %% Computing the Target-wise averaging of total available sites
    % potential cutoff -> probably startfreq should go into settings at somepoint
    startfreq=3.5;
    colsbp=jet(length(cfg.lfp.frequency_bands));
    fbx=find(cfg.lfp.frequency_bands(:,1)>=startfreq,1,'first');
    frequency_bands=cfg.lfp.frequency_bands(fbx:end,:);
    %freqName = cfg.lfp.freqName(fbx:end);
    colsbp=colsbp(fbx:end,:);
    freqb=num2cell(strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz'),2);
    fx=find(cfg.lfp.foi>=startfreq,1,'first');
    freq = cfg.lfp.foi(fx:end);
    
    %%
    tic
    for t = 1: length(targets)
        sites_for_this_target=arrayfun(@(x) any(strfind(x.target,targets{t})),sites);
        target_sites = sites(sites_for_this_target);
        condition=vertcat(target_sites.condition);
        
        for c = 1:size(condition,2)
            events=cat(1,condition(:,c).event);
            
            
            for e=1:size(cfg.analyse_states,1)
                E               = [];
                eval(['E.',Fin{fin},'           = [];']);
                eval(['E.',Fin{fin},'_time      = [];']);
                eval(['E.',Fin{fin},'_std       = [];']);
                eval(['E.',Fin{fin},'_25        = [];']);
                eval(['E.',Fin{fin},'_75        = [];']);
                eval(['E.',Fin{fin},'_sig       = [];']);
                eval(['E.',Fin{fin},'_sig_signed= [];']);
                
                E.nTriggers     = mean([events(:,e).nTriggers]);
                eval(['E.',Fin{fin},'           = (cat(3,events(:,e).',Fin{fin},'));']);
                eval(['E.',Fin{fin},'_25        = prctile(cat(3,events(:,e).',Fin{fin},'),25,3);']);
                eval(['E.',Fin{fin},'_75        = prctile(cat(3,events(:,e).',Fin{fin},'),75,3);']);
                eval(['E.',Fin{fin},'_time      = events(1,e).time;']);
                eval(['E.',Fin{fin},'_std       = std(cat(3,events(:,e).',Fin{fin},'),0,3);']);
                eval(['E.',Fin{fin},'_sig       = mean(cat(3,events(:,e).',Fin{fin},'_sig),3);']);
                eval(['E.',Fin{fin},'_sig_signed= mean(cat(3,events(:,e).',Fin{fin},'_sig).*sign(cat(3,events(:,e).',Fin{fin},')),3);']);
                
                eval(['E.',Fin{fin},'_sig       = E.',Fin{fin},'_sig*100;']);
                eval(['E.',Fin{fin},'_sig_signed= E.',Fin{fin},'_sig_signed*100;']);
                
                
                eventname = cfg.analyse_states{e,1};
                tar(t).con(c).(eventname)=E;
            end
            tar(t).con(c).cond_name = cond{c};
        end
        tar(t).target = targets{t};
        tar(t).nSites = length(target_sites);
        tar(t).site_info = arrayfun(@(x) (x.site_info),target_sites);
    end
    toc
    %%
    for t = 1: length(targets)
        all_target_session = unique(arrayfun(@(x) (num2str(x.date)),tar(t).site_info,'UniformOutput',false));
        for ses = 1: length(all_target_session)
            session = all_target_session{ses};
            sites_for_this_session = find(contains({tar(t).site_info.site_ID},session));
            [~,depth_order_idx] = sort([tar(t).site_info(sites_for_this_session).electrode_depth]);
            h(ses) = figure;
            h(ses).WindowState = 'maximized';
            for si = 1: length(sites_for_this_session)
                sp = round(sqrt(length(sites_for_this_session)))+1;
                subplot(sp,sp,si)
                
                idx = sites_for_this_session(depth_order_idx(si));
                
                eval(['site_rest = squeeze(tar(t).con(1).R.',Fin{fin},'(:,:,idx));']);
                eval(['site_task = squeeze(tar(t).con(2).R.',Fin{fin},'(:,:,idx));']);
                eval(['plot_time = tar(t).con(1).R.',Fin{fin},'_time;']);
                
                data_rest_smth = smoothit(site_rest',half_win,gaussian_kernel);
                data_task_smth = smoothit(site_task',half_win,gaussian_kernel);
                
                plot(plot_time,data_rest_smth,'color','b','LineWidth',1.5),
                hold on,
                plot(plot_time,data_task_smth,'color','r','LineWidth',1.5),
                axis square,
%                 xlabel('Time(s)');
                line([0 0], ylim, 'color', 'k');
                title(['ch: ', num2str(tar(t).site_info(idx).channel) ,' depth : ', num2str(tar(t).site_info(idx).electrode_depth)],'FontSize',6)
%                 if strcmp(Fin{fin},'mua')
%                     ylabel('MUA evoked (A.U.)','FontSize',6);
%                 else
%                     ylabel('LFP evoked (V)','FontSize',6);
%                 end
            end
            sgtitle([monkey,'-',targets{t},'-',Fin{fin},'- Per electrode depth - session ', session],'fontsize', 10 ,'Interpreter', 'none');
            result_file = fullfile([path_to_save, filesep,monkey,'-',targets{t},'-',Fin{fin},'- Per electrode depth - session ', session]);
            export_fig(h(ses),[result_file,'.pdf']);
        end
        close all
    end
end
end