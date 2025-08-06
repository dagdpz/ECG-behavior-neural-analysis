function ecg_bna_LFP_and_MUA_sig_plotAllTargets
% ==================================================================
% just to plot the LFP and MUA sig in the squared axis all target in one
% plot >> 26.07.25
% =================================================================
close all,
clear all
clc

monkey = 'Magnus';%'Magnus' or 'Bacchus'

Fin = {'lfp','mua'};
for fin = 1:length(Fin)
    withunits = 'all_sites';
    
    if strcmp(Fin{fin},'mua')
        path_to_save = ['Y:\Projects\Pulv_bodysignal\MUA\ECG_',monkey,'_TaskRest_generalTrig_MUA\grand_avg_all\'];
        perSitedata_path = ['Y:\Projects\Pulv_bodysignal\MUA\ECG_',monkey,'_TaskRest_generalTrig_MUA\Per_Site'];
    elseif strcmp(Fin{fin},'lfp')
        path_to_save = ['Y:\Projects\Pulv_bodysignal\LFP\ECG_',monkey,'_TaskRest_generalTrig_finalRun\not filtered\grand_average_wo combine_hemispheres_LFP real_New GrandAvg version\'];
        perSitedata_path = ['Y:\Projects\Pulv_bodysignal\LFP\ECG_',monkey,'_TaskRest_generalTrig_finalRun\not filtered\Per_Site'];
        
    end
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
    end
    %%
    
    for e=1:size(cfg.analyse_states,1)
        h(e) = figure;
        h(e).WindowState = 'maximized';
        for t = 1: length(targets)
            
            
            data_rest = squeeze(tar(t).con(1).(eventname).([Fin{fin},'_sig'])); % grand adverage of MD target, Rest Condition , across grand_avg(2).nSites;
            data_task = squeeze(tar(t).con(2).(eventname).([Fin{fin},'_sig'])); % grand adverage of MD target, Task Condition ,across grand_avg(2).nSites;
            eval(['plot_time = tar(t).con(1).(eventname).',Fin{fin},'_time;']);
            
            data_rest_smth = smoothit(data_rest',half_win,gaussian_kernel);
            data_task_smth = smoothit(data_task',half_win,gaussian_kernel);
            % ==============================================================
            % plotting the significance
            
            subplot(1,3,t)
            plot(plot_time,data_rest_smth','color','b','LineWidth',1.5) , hold on
            plot(plot_time,data_task_smth','color','r','LineWidth',1.5) , hold on
            line([0 0], ylim, 'color', 'k');
            
            if strcmp(Fin{fin},'mua') && strcmp(monkey,'Magnus')
                y1 = [0 50];
            elseif strcmp(Fin{fin},'mua') && strcmp(monkey,'Bacchus')
                y1 = [0 100];
            elseif strcmp(Fin{fin},'lfp') && strcmp(monkey,'Magnus')
                y1 = [0 100];
            elseif strcmp(Fin{fin},'lfp') && strcmp(monkey,'Bacchus')
                y1 = [0 100];
            end
            
            axis square;
            set(gca,'ylim',y1);
            xlabel('Time(s)');
            line([0 0], y1, 'color', 'k');
            title([monkey,'-',targets{t},'-',Fin{fin},' - nSite: ', num2str(tar(t).nSites)],'fontsize', 12,'Interpreter', 'none');
            if strcmp(Fin{fin},'mua')
                ylabel('sig. MUA [% of sites]');
            else
                ylabel('sig. LFP [% of sites]');
            end
            results_file = fullfile([path_to_save,filesep,monkey,'- allTarget -',Fin{fin},'-',cond{2},' vs. ',cond{1},'_Sig_site_percentage']);
            
            
            
        end
        export_fig(h(e),[results_file,'.pdf']);
    end
end
close all
end