function tar = ecg_bna_compute_grand_avg_mua(cfg,withunits)

reprocess=1;
reprocess2=1;
fileName = fullfile([cfg.fldr.MUA_root filesep 'Grand_avg_mua_sites_',withunits,'.mat']);
fileName2 = fullfile([cfg.fldr.MUA_root filesep 'Grand_avg_mua_targets_',withunits,'.mat']);

if cfg.combine_hemispheres
    cfg.targets=unique(cellfun(@(x) x(1:strfind(x,'_')-1),cfg.targets,'uniformoutput',false));
end
targets = cfg.targets;%unique({sites.target}); %%cfg.monkey

if reprocess
    data_path = cfg.fldr.MUA_sites;
%     data_path = 'Y:\Projects\Pulv_bodysignal\LFP\ECG_Magnus_TaskRest_test_generalTrig_Reref_shamim\sample session check';
    cd(data_path)
    % read the files
    all_mua_data = dir('*.mat');
    % storing all the Triggered parameters of the sites in the results folder:
    sites = struct;
    s=0;
    for f = 1:length(all_mua_data)
        load(all_mua_data(f).name)
        session = triggered_site_data.session;
        site_ID = triggered_site_data.site_ID;
        switch withunits
            case 'Cue_responsive'
                if ~ismember(site_ID,cfg.site_IDS)
                    continue
                end
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
                
                nTriggers = event.observed.ntriggers;
                time      = event.time;
                evoked_time  = event.tfr_time;
                evoked       = squeeze(event.normalized.evoked.mean)';
                evoked_observed       = squeeze(event.observed.evoked.mean)';
                evoked_observed_std       = squeeze(event.observed.evoked.std)';
                evoked_surrogate       = squeeze(event.surrogate.evoked.mean)';
                evoked_surrogate_std       = squeeze(event.surrogate.evoked.std)';
                
%                 p=event.prevspost.p;
%                 h=event.prevspost.h;
                
                evoked_sig     =squeeze(abs(event.significance.evoked))'; %% weird dimension thing
                                
                sites(s).condition(c).event(e).nTriggers = nTriggers;
                sites(s).condition(c).event(e).evoked = evoked;               
                sites(s).condition(c).event(e).evoked_sig = evoked_sig;                  
                sites(s).condition(c).event(e).any_sig = any(evoked_sig);   
                sites(s).condition(c).event(e).evoked_observed = evoked_observed;               
                sites(s).condition(c).event(e).evoked_observed_std = evoked_observed_std;               
                sites(s).condition(c).event(e).evoked_surrogate = evoked_surrogate;               
                sites(s).condition(c).event(e).evoked_surrogate_std = evoked_surrogate_std;        
%                 sites(s).condition(c).event(e).p = p;       
%                 sites(s).condition(c).event(e).h = h;                      
                
                sites(s).condition(c).event(e).time = time;
                sites(s).condition(c).event(e).tfr_time = evoked_time;
                
                %% some comparison of before and after trigger
            end
        end
    end
    
    %% this part to remove sites with not all conditions
    % out_mask = cellfun(@isempty, {sites.site_ID});
    % sites = sites(~out_mask);
    sites = sites([sites.all_conditions_present]==1);
    save(fileName,'sites')
    
else
    load(fileName);
end

%% Computing the Target-wise averaging of total available sites
% for t = 1: length(targets)
%     sites_for_this_target=arrayfun(@(x) any(strfind(x.target,targets{t})),sites);
%     target_sites = sites(sites_for_this_target);
%     condition=vertcat(target_sites.condition);
%     
%     for c = 1:size(condition,2)
%         events=cat(1,condition(:,c).event);
%         
%         concat.nTriggers=round(mean(reshape([events.nTriggers],size(events)),1));
%         
%         concat.evoked          = [];
%         concat.evoked_25       = [];
%         concat.evoked_75       = [];
%         concat.evokedsig       = [];
%         concat.evoked_time     = [];
%         concat.evoked_std      = [];
%         
%         for e=1:size(cfg.analyse_states,1)
%             E.nTriggers     = mean(cat(3,events(:,e).nTriggers),3);
%             E.evoked           = mean(cat(3,events(:,e).evoked),3);
%             E.evoked_25        = prctile(cat(3,events(:,e).evoked),25,3);
%             E.evoked_75        = prctile(cat(3,events(:,e).evoked),75,3);
%             E.evoked_time      = events(1,e).time;
%             E.evoked_std       = std(cat(3,events(:,e).evoked),0,3);                        
%             E.evoked_sig       = mean(cat(3,events(:,e).evoked_sig),3);            
%             E.evoked_sig_signed= mean(cat(3,events(:,e).evoked_sig).*sign(cat(3,events(:,e).evoked)),3);
%             
%             E.evoked_sig       = E.evoked_sig*100;                
%             E.evoked_sig_signed= E.evoked_sig_signed*100;     
%             
%             NaNseparator=100/25;
%             concat.evoked          = cat(2, concat.evoked,        E.evoked,        nan(size(E.evoked, 1),        NaNseparator));
%             concat.evoked_25       = cat(2, concat.evoked_25,     E.evoked_25,	  nan(size(E.evoked_25, 1),     NaNseparator));
%             concat.evoked_75       = cat(2, concat.evoked_75,     E.evoked_75,	  nan(size(E.evoked_75, 1),     NaNseparator));
%             
%             concat.evokedsig       = cat(2, concat.evokedsig,     E.evoked_sig,    nan(size(E.evoked_sig, 1),    NaNseparator));
%             concat.evoked_std      = cat(2, concat.evoked_std,    E.evoked_std,	  nan(size(E.evoked_std, 1),    NaNseparator));
%             concat.evoked_time     = [concat.evoked_time, E.evoked_time, nan(1, NaNseparator)];
%             
%             eventname=cfg.analyse_states{e,1};
%             tar(t).con(c).(eventname)=E;
%         end
%         tar(t).con(c).cond_name = cfg.condition(c).name;
%         tar(t).con(c).concat=concat;
%     end
%     tar(t).target = targets{t};
%     tar(t).nSites = length(target_sites);
% end

if reprocess2
    for t = 1: length(targets)
        sites_for_this_target=arrayfun(@(x) any(strfind(x.target,targets{t})),sites);
        target_sites = sites(sites_for_this_target);
        condition=vertcat(target_sites.condition);
        
        taskevents=cat(1,condition(:,2).event);
        restevents=cat(1,condition(:,1).event);
        
        
%         cueresponse_sigbin=[taskevents(:,4).any_sig];
%         cueresponse_ttest=[taskevents(:,4).p]<0.05;
        
        
        
        for c = 1:size(condition,2)
            events=cat(1,condition(:,c).event);
            
            concat.nTriggers=round(mean(reshape([events.nTriggers],size(events)),1));
            
            concat.evoked          = [];
            concat.evoked_abs      = [];
            concat.evoked_abs_25   = [];
            concat.evoked_abs_75   = [];
            concat.evoked_25       = [];
            concat.evoked_75       = [];
            concat.evokedsig       = [];
            concat.evoked_time     = [];
            concat.evoked_std      = [];
            concat.evoked_abs_std      = [];
            
            for e=1:size(cfg.analyse_states,1)
                
                sig_sites=[taskevents(:,e).any_sig] | [restevents(:,e).any_sig];
                %S=sig_sites;
                S=true(size(sig_sites));
                %S=cueresponse_sigbin;
                %S=cueresponse_ttest;
                
                
                %             evoked=zscore(cat(3,events(S,e).evoked),0,2);
                %             evoked_observed=zscore(cat(3,events(S,e).evoked_observed),0,2);
                %             evoked_surrogate=zscore(cat(3,events(S,e).evoked_surrogate),0,2);
                evoked=cat(3,events(S,e).evoked);
                evoked_observed=cat(3,events(S,e).evoked_observed);
                evoked_surrogate=cat(3,events(S,e).evoked_surrogate);
                
                
                E.evoked_p= ecg_bna_population_significant_n_bins(cat(3,events(:,e).evoked_sig));
                E.nTriggers     = mean(cat(3,events(S,e).nTriggers),3);
                E.evoked           = mean(evoked,3);
                E.evoked_observed      = mean(evoked_observed,3);
                E.evoked_surrogate  = mean(evoked_surrogate,3);
                E.evoked_abs       = mean(abs(evoked),3);
                E.evoked_abs_sterr = sterr(abs(evoked),3);
                E.evoked_abs_25    = prctile(abs(evoked),25,3);
                E.evoked_abs_75    = prctile(abs(evoked),75,3);
                E.evoked_25        = prctile(evoked,25,3);
                E.evoked_75        = prctile(evoked,75,3);
                E.evoked_time      = events(1,e).time;
                E.evoked_std       = std(evoked,0,3);
                E.evoked_observed_std  = std(evoked_observed,0,3);
                E.evoked_sterr     = sterr(evoked,3);
                E.evoked_observed_sterr= sterr(evoked_observed,3);
                E.evoked_observed_abs  = mean(abs(evoked_observed),3);
                E.evoked_observed_abs_std  = std(abs(evoked_observed),0,3);
                E.evoked_observed_abs_sterr  = sterr(abs(evoked_observed),3);
                E.evoked_sig       = mean(cat(3,events(S,e).evoked_sig),3);
                E.evoked_abs_std   = std(abs(evoked),0,3);
                E.evoked_sig_signed= mean(cat(3,events(S,e).evoked_sig).*sign(cat(3,events(S,e).evoked)),3);
                [E.evoked_popsig, E.evoked_poppval]=ecg_bna_population_significance_normalized_parametric(evoked);
                   
                
                E.evoked_sig       = E.evoked_sig*100;
                E.evoked_sig_signed= E.evoked_sig_signed*100;
                
                NaNseparator=100/25;
                concat.evoked          = cat(2, concat.evoked,        E.evoked,        nan(size(E.evoked, 1),        NaNseparator));
                concat.evoked_25       = cat(2, concat.evoked_25,     E.evoked_25,	  nan(size(E.evoked_25, 1),     NaNseparator));
                concat.evoked_75       = cat(2, concat.evoked_75,     E.evoked_75,	  nan(size(E.evoked_75, 1),     NaNseparator));
                concat.evoked_abs      = cat(2, concat.evoked_abs,    E.evoked_abs,    nan(size(E.evoked_abs, 1),    NaNseparator));
                concat.evoked_abs_25   = cat(2, concat.evoked_abs_25, E.evoked_abs_25, nan(size(E.evoked_abs_25, 1), NaNseparator));
                concat.evoked_abs_75   = cat(2, concat.evoked_abs_75, E.evoked_abs_75, nan(size(E.evoked_abs_75, 1), NaNseparator));
                
                concat.evokedsig       = cat(2, concat.evokedsig,     E.evoked_sig,    nan(size(E.evoked_sig, 1),    NaNseparator));
                concat.evoked_std      = cat(2, concat.evoked_std,    E.evoked_std,	  nan(size(E.evoked_std, 1),    NaNseparator));
                concat.evoked_abs_std  = cat(2, concat.evoked_abs_std,E.evoked_abs_std,nan(size(E.evoked_abs_std, 1),NaNseparator));
                concat.evoked_time     = [concat.evoked_time, E.evoked_time, nan(1, NaNseparator)];
                
                eventname=cfg.analyse_states{e,1};
                tar(t).con(c).(eventname)=E;
            end
            tar(t).con(c).cond_name = cfg.condition(c).name;
            tar(t).con(c).concat=concat;
        end
        tar(t).target = targets{t};
        tar(t).nSites = sum(S);
    end
    save(fileName2,'tar')
else
    load(fileName2);
end
%% plotting the results:
cond = {cfg.condition.name};
plot_names={'MUA','MUA significance','MUA signed_signficance'};

% Smoothing Kernel here:
win = 1:cfg.lfp.smoothWin; win=win-(numel(win)+1)/2;
half_win = ceil(size(win,2)/2)-1;
gaussian_kernel=normpdf(win,0,numel(win)/6);
gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);


% figure 1 - overview
for t = 1: length(targets)
    if tar(t).nSites == 0
        continue;
    end
    for c = 1:length(cond)
        evoked_time=tar(t).con(c).concat.evoked_time;
        tfr_events.onset        = find(evoked_time == 0);
        tfr_events.name         = cfg.analyse_states(:,1);
        tfr_events.ticksamples  = sort([find(diff([0 ~isnan(evoked_time)])==1) find(diff(~isnan(evoked_time))==-1)]);
        tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
        tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
        tfr_events.ticks        = round(evoked_time(tfr_events.ticksamples)*100)/100;
        
        evoked=tar(t).con(c).concat.evoked;
        evoked_std=tar(t).con(c).concat.evoked_std;
        %evoked_se=sterr(evoked,3,0);
        percentile25=tar(t).con(c).concat.evoked_25;
        percentile75=tar(t).con(c).concat.evoked_75;
        percentile25(isnan(percentile25))=0;
        percentile75(isnan(percentile75))=0;
        
        % create figure
        h(c,t) = figure('units','normalized','position',[0 0 1 1]);
        
        %========================== LFP evoked Potential ==================== %
        
        smoothed_mean=smoothit(evoked,half_win,gaussian_kernel);
        smoothed_std=smoothit(evoked_std,half_win,gaussian_kernel);
        smoothed_25=smoothit(percentile25,half_win,gaussian_kernel);
        smoothed_75=smoothit(percentile75,half_win,gaussian_kernel);
            
        sp=1;
        sph(c,sp,t)=subplot(2,1,sp);
        ax=sph(c,sp,t);
        hold on;
        lineProps={'color',[0 0 1]};
        %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,[smoothed_75-smoothed_mean;smoothed_mean-smoothed_25 ],lineProps,1);        
        shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
        %plot(repmat(time,size(tar(t).con(c).evoked_avg,1),1)', squeeze(smoothed_mean)')
        
        xlabel('Time(s)'); ylabel('MUA evoked Potential');
        title(plot_names{sp},'fontsize',10,'interpreter','none');
        add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)        
        
        results_file{c} = fullfile(cfg.fldr.MUA_root, [cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-','avg of ',num2str(tar(t).nSites),' sites ',withunits]);
        mtit(h(c,t),[ cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-avg of ',num2str(tar(t).nSites),' sites ' ,withunits],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
        
        export_fig(h(c,t),[results_file{c},'.pdf']);
    end
end
close all,
clear results_file sph h collim


% same as figure 1, but separately for each event
plot_names={'MUA','MUA significance','MUA absolute','MUA signed_signficance'};
for e=1:size(cfg.analyse_states,1)
    E=cfg.analyse_states{e,1};
    
    for t = 1: length(targets)
        for c = 1:length(cond)
            % create figure
            h(c,t) = figure('units','normalized','position',[0 0 1 1]);
            if tar(t).nSites == 0
                continue;
            end
            
            evoked_time                = tar(t).con(c).(E).evoked_time;
            tfr_events.onset        = find(evoked_time == 0);
            tfr_events.name         = cfg.analyse_states(e,1);
            tfr_events.ticksamples  = [1,numel(evoked_time)]; 
            tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
            tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
            tfr_events.ticks        = round(evoked_time(tfr_events.ticksamples)*100)/100;
            
            %========================== MUA evoked Potential ==================== %
            evoked=tar(t).con(c).(E).evoked;
            evoked_abs=tar(t).con(c).(E).evoked_abs;
            evoked_std=tar(t).con(c).(E).evoked_std;
            evoked_sterr=tar(t).con(c).(E).evoked_sterr;
            evoked_abs_std=tar(t).con(c).(E).evoked_abs_std;
            evoked_abs_sterr=tar(t).con(c).(E).evoked_abs_sterr;
            %evoked_se=sterr(evoked,3,0);
            percentile25=tar(t).con(c).(E).evoked_25;
            percentile75=tar(t).con(c).(E).evoked_75;
            percentile25(isnan(percentile25))=0;
            percentile75(isnan(percentile75))=0;
           
            % Smoothing of the  LFP evoked Potential here:
            smoothed_mean=smoothit(evoked,half_win,gaussian_kernel);    
            smoothed_std=smoothit(evoked_std,half_win,gaussian_kernel);
            smoothed_sterr=smoothit(evoked_sterr,half_win,gaussian_kernel);
            smoothed_25=smoothit(percentile25,half_win,gaussian_kernel);
            smoothed_75=smoothit(percentile75,half_win,gaussian_kernel);
            
            %========================== evoked  ==================== %
            sp=1;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            lineProps={'color',[0 0 1]};
            
            %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
            shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_sterr,lineProps,1);
            %plot(repmat(time,size(tar(t).con(c).evoked_avg,1),1)', squeeze(smoothed_mean)')
            xlabel('Time(s)'); ylabel('MUA evoked Potential');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)            
            
            
            %========================== evoked significance ==================== %
            evoked_sig=tar(t).con(c).(E).evoked_sig;
            smoothed_mean=smoothit(evoked_sig,half_win,gaussian_kernel);          
            
            sp=2;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
            
            xlabel('Time(s)'); ylabel('sig. MUA [% of sites]');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)   
            
            
            %========================== evoked abs  ==================== %
            smoothed_mean=smoothit(evoked_abs,half_win,gaussian_kernel);    
            smoothed_std=smoothit(evoked_abs_std,half_win,gaussian_kernel);
            smoothed_sterr=smoothit(evoked_abs_sterr,half_win,gaussian_kernel);
            sp=3;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            lineProps={'color',[0 0 1]};
            
            %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
            shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_sterr,lineProps,1);
            %plot(repmat(time,size(tar(t).con(c).evoked_avg,1),1)', squeeze(smoothed_mean)')
            xlabel('Time(s)'); ylabel('MUA Absolute');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)            
            
            
            %========================== evoked significance ==================== %
            evoked_sig=tar(t).con(c).(E).evoked_sig_signed;
            smoothed_mean=smoothit(evoked_sig,half_win,gaussian_kernel);          
            
            sp=4;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
            
            xlabel('Time(s)'); ylabel('signed sig. MUA [% of sites]');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)   
            
            
            
            results_file{c,t} = fullfile(cfg.fldr.MUA_root, [cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-','avg of ',num2str(tar(t).nSites),' sites ',withunits]);
            mtit(h(c,t),[ cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-avg of ',num2str(tar(t).nSites),' sites ' ,withunits],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
            export_fig(h(c,t),[results_file{c,t},'.pdf']);
        end
        
    end
    
end
close all,

clear results_file sph h collim



% same as figure 1, but separately for each event
plot_names={'MUA','MUA significance','MUA absolute','MUA signed_signficance'};
for e=1:size(cfg.analyse_states,1)
    E=cfg.analyse_states{e,1};
    
    for t = 1: length(targets)
        for c = 1:length(cond)
            % create figure
            h(c,t) = figure('units','normalized','position',[0 0 1 1]);
            if tar(t).nSites == 0
                continue;
            end
            
            evoked_time                = tar(t).con(c).(E).evoked_time;
            tfr_events.onset        = find(evoked_time == 0);
            tfr_events.name         = cfg.analyse_states(e,1);
            tfr_events.ticksamples  = [1,numel(evoked_time)]; 
            tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
            tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
            tfr_events.ticks        = round(evoked_time(tfr_events.ticksamples)*100)/100;
            
            %========================== MUA evoked Potential ==================== %
            evoked=tar(t).con(c).(E).evoked_observed;
            evoked_abs=tar(t).con(c).(E).evoked_observed_abs;
            evoked_std=tar(t).con(c).(E).evoked_observed_std;
            evoked_sterr=tar(t).con(c).(E).evoked_observed_sterr;
            evoked_abs_std=tar(t).con(c).(E).evoked_observed_abs_std;
            evoked_abs_sterr=tar(t).con(c).(E).evoked_observed_abs_sterr;
            %evoked_se=sterr(evoked,3,0);
            percentile25=tar(t).con(c).(E).evoked_25;
            percentile75=tar(t).con(c).(E).evoked_75;
            percentile25(isnan(percentile25))=0;
            percentile75(isnan(percentile75))=0;
           
            % Smoothing of the  LFP evoked Potential here:
            smoothed_mean=smoothit(evoked,half_win,gaussian_kernel);    
            smoothed_std=smoothit(evoked_std,half_win,gaussian_kernel);
            smoothed_sterr=smoothit(evoked_sterr,half_win,gaussian_kernel);
            smoothed_25=smoothit(percentile25,half_win,gaussian_kernel);
            smoothed_75=smoothit(percentile75,half_win,gaussian_kernel);
            
            %========================== evoked  ==================== %
            sp=1;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            lineProps={'color',[0 0 1]};
            
            %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
            shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_sterr,lineProps,1);
            %plot(repmat(time,size(tar(t).con(c).evoked_avg,1),1)', squeeze(smoothed_mean)')
            xlabel('Time(s)'); ylabel('MUA evoked Potential');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)            
            
            
            %========================== evoked significance ==================== %
            evoked_sig=tar(t).con(c).(E).evoked_sig;
            smoothed_mean=smoothit(evoked_sig,half_win,gaussian_kernel);          
            
            sp=2;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
            
            xlabel('Time(s)'); ylabel('sig. MUA [% of sites]');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)   
            
            
            %========================== evoked abs  ==================== %
            smoothed_mean=smoothit(evoked_abs,half_win,gaussian_kernel);    
            smoothed_std=smoothit(evoked_abs_std,half_win,gaussian_kernel);
            smoothed_sterr=smoothit(evoked_abs_sterr,half_win,gaussian_kernel);
            sp=3;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            lineProps={'color',[0 0 1]};
            
            %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
            shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_sterr,lineProps,1);
            %plot(repmat(time,size(tar(t).con(c).evoked_avg,1),1)', squeeze(smoothed_mean)')
            xlabel('Time(s)'); ylabel('MUA Absolute');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)            
            
            
            %========================== evoked significance ==================== %
            evoked_sig=tar(t).con(c).(E).evoked_sig_signed;
            smoothed_mean=smoothit(evoked_sig,half_win,gaussian_kernel);          
            
            sp=4;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
            
            xlabel('Time(s)'); ylabel('signed sig. MUA [% of sites]');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)   
            
            
            
            results_file{c,t} = fullfile(cfg.fldr.MUA_root, ['observed_' cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-','avg of ',num2str(tar(t).nSites),' sites ',withunits]);
            mtit(h(c,t),[ cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-avg of ',num2str(tar(t).nSites),' sites ' ,withunits],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
            export_fig(h(c,t),[results_file{c,t},'.pdf']);
        end
        
    end
    
end
close all,


% MUA significance summary figure
for e=1:size(cfg.analyse_states,1)
    E=cfg.analyse_states{e,1};
    
    h(e) = figure('units','normalized','position',[0 0 1 1]);
    k=0;
    for t = 1: length(targets)
        for c = 1:length(cond)
            k=k+1;
            % create figure
            if tar(t).nSites == 0
                continue;
            end
            
            evoked_time                = tar(t).con(c).(E).evoked_time;
            evoked_p                   = tar(t).con(c).(E).evoked_p;
            if evoked_p<0.001
                ptoplot='<0.001';
            else
                ptoplot=num2str(round(evoked_p*1000)/1000);
            end
            tfr_events.onset        = find(evoked_time == 0);
            tfr_events.name         = cfg.analyse_states(e,1);
            tfr_events.ticksamples  = [1,numel(evoked_time)]; 
            tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
            tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
            tfr_events.ticks        = round(evoked_time(tfr_events.ticksamples)*100)/100;      
            
            
            %========================== evoked significance ==================== %
            evoked_sig=tar(t).con(c).(E).evoked_sig;
            smoothed_mean=smoothit(evoked_sig,half_win,gaussian_kernel);          
            
            sp=2;
            sph(c,sp,t)=subplot(3,2,k);
            ax=sph(c,sp,t);
            hold on;
            plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
            
            xlabel('Time(s)'); ylabel('sig. MUA [% of sites]');
            title([targets{t} '-' cond{c} '-' num2str(tar(t).nSites) 'sites-p=' ptoplot],'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)   
            
        end
    end
    results_file = fullfile(cfg.fldr.MUA_root, [cfg.monkey,'-',E,'-sigNbinstesting-',withunits]);
    export_fig(h(e),[results_file,'.pdf']);
end
close all,


clear results_file sph h collim

clc
% 
% % same as figure 1, but now for significances
% for e=1:size(cfg.analyse_states,1)
%     E=cfg.analyse_states{e,1};
%     
%     for t = 1: length(targets)
%         for c = 1:length(cond)
%             k=t+(c-1)*length(targets);
%             bins=tar(t).con(c).(E).tfr_time;
%             % create figure
%             h(c,t) = figure('units','normalized','position',[0 0 1 1]);
%             if tar(t).nSites == 0
%                 continue;
%             end
%             
%             % tfr_time not defined well!
%             evoked_time=tar(t).con(c).(E).tfr_time;
%             %evoked_time=tar(t).con(c).concat.evoked_time;
%             tfr_events.onset        = find(evoked_time == 0);
%             tfr_events.name         = cfg.analyse_states(e,1);
%             tfr_events.ticksamples  = [1,numel(evoked_time)]; %sort([find(diff([0 ~isnan(tfr_time)])==1) find(diff([~isnan(tfr_time)])==-1)]);
%             tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
%             tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
%             tfr_events.ticks        = round(evoked_time(tfr_events.ticksamples)*100)/100;
%             
%             evoked=tar(t).con(c).(E).evoked_sig;
%             evoked_std=tar(t).con(c).(E).evoked_std;
%             %evoked_se=sterr(evoked,3,0);
%             percentile25=tar(t).con(c).(E).evoked_25;
%             percentile75=tar(t).con(c).(E).evoked_75;
%             percentile25(isnan(percentile25))=0;
%             percentile75(isnan(percentile75))=0;
%            
%             
%             
%             %========================== LFP evoked Potential ==================== %
%             % Smoothing of the  LFP evoked Potential here:
%             jnk = [];
%             concat_input = cat(2,(evoked(:,half_win:-1:1)),(evoked(:,:)));
%             concat_input = cat(2,concat_input, (evoked(:,end:-1:end-half_win+1)));
%             for k=1:size(evoked,1)
%                 jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
%             end
%             clear concat_input
%             smoothed_mean = jnk(:,half_win+1:end-half_win);
%             concat_input = cat(2,(evoked_std(:,half_win:-1:1)),(evoked_std(:,:)));
%             concat_input = cat(2,concat_input, (evoked_std(:,end:-1:end-half_win+1)));
%             for k=1:size(evoked,1)
%                 jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
%             end
%             clear concat_input
%             smoothed_std = jnk(:,half_win+1:end-half_win);
%             concat_input = cat(2,(percentile25(:,half_win:-1:1)),(percentile25(:,:)));
%             concat_input = cat(2,concat_input, (percentile25(:,end:-1:end-half_win+1)));
%             for k=1:size(evoked,1)
%                 jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
%             end
%             clear concat_input
%             smoothed_25 = jnk(:,half_win+1:end-half_win);
%             concat_input = cat(2,(percentile75(:,half_win:-1:1)),(percentile75(:,:)));
%             concat_input = cat(2,concat_input, (percentile75(:,end:-1:end-half_win+1)));
%             for k=1:size(evoked,1)
%                 jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
%             end
%             clear concat_input
%             smoothed_75 = jnk(:,half_win+1:end-half_win);
%             
%             sp=5;
%             sph(c,sp,t)=subplot(3,2,sp);
%             ax=sph(c,sp,t);
%             hold on;
%             plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
%             
%             %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
%             %plot(repmat(time,size(tar(t).con(c).evoked_avg,1),1)', squeeze(smoothed_mean)')
%             xlabel('Time(s)'); ylabel('sig. LFP [% of sites]');
%             title(plot_names{sp},'fontsize',10,'interpreter','none');
%             add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
%             
%             
%             results_file{c,t} = fullfile(cfg.analyse_evoked_folder, [cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-','sig% of ',num2str(tar(t).nSites),' sites ',withunits]);
%             mtit(h(c,t),[cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-sig% of ',num2str(tar(t).nSites),' sites ' ,withunits],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
%             
%         end
%         
%     end
%     
%     for c = 1:length(cond)
%         for t=1:length(targets)
%         figure(h(c,t));
%         for sp=1:2
%             subplot(sph(c,sp,t));
%             set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
%         end
%         export_fig(h(c,t),[results_file{c,t},'.pdf']);
%         end
%     end
% end
% close all,
% 
% clear results_file sph h collim
% clc
% 
% 
% % same as figure 1, but now for SIGNED significances
% for e=1:size(cfg.analyse_states,1)
%     E=cfg.analyse_states{e,1};
%     
%     for t = 1: length(targets)
%         for c = 1:length(cond)
%             k=t+(c-1)*length(targets);
%             bins=tar(t).con(c).(E).tfr_time;
%             % create figure
%             h(c,t) = figure('units','normalized','position',[0 0 1 1]);
%             if tar(t).nSites == 0
%                 continue;
%             end
%             
%             % tfr_time not defined well!
%             evoked_time=tar(t).con(c).(E).tfr_time;
%             %evoked_time=tar(t).con(c).concat.evoked_time;
%             tfr_events.onset        = find(evoked_time == 0);
%             tfr_events.name         = cfg.analyse_states(e,1);
%             tfr_events.ticksamples  = [1,numel(evoked_time)]; %sort([find(diff([0 ~isnan(tfr_time)])==1) find(diff([~isnan(tfr_time)])==-1)]);
%             tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
%             tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
%             tfr_events.ticks        = round(evoked_time(tfr_events.ticksamples)*100)/100;
%             
%             powbp=tar(t).con(c).(E).powbp_sig_signed;
%             itpcbp=tar(t).con(c).(E).itpcbp_sig_signed;
%             evoked=tar(t).con(c).(E).evoked_sig_signed;
%             evoked_std=tar(t).con(c).(E).evoked_std;
%             %evoked_se=sterr(evoked,3,0);
%             percentile25=tar(t).con(c).(E).evoked_25;
%             percentile75=tar(t).con(c).(E).evoked_75;
%             percentile25(isnan(percentile25))=0;
%             percentile75(isnan(percentile75))=0;
%             
%             
%             %========================== LFP evoked Potential ==================== %
%             % Smoothing of the  LFP evoked Potential here:
%             jnk = [];
%             concat_input = cat(2,(evoked(:,half_win:-1:1)),(evoked(:,:)));
%             concat_input = cat(2,concat_input, (evoked(:,end:-1:end-half_win+1)));
%             for k=1:size(evoked,1)
%                 jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
%             end
%             clear concat_input
%             smoothed_mean = jnk(:,half_win+1:end-half_win);
%             concat_input = cat(2,(evoked_std(:,half_win:-1:1)),(evoked_std(:,:)));
%             concat_input = cat(2,concat_input, (evoked_std(:,end:-1:end-half_win+1)));
%             for k=1:size(evoked,1)
%                 jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
%             end
%             clear concat_input
%             smoothed_std = jnk(:,half_win+1:end-half_win);
%             concat_input = cat(2,(percentile25(:,half_win:-1:1)),(percentile25(:,:)));
%             concat_input = cat(2,concat_input, (percentile25(:,end:-1:end-half_win+1)));
%             for k=1:size(evoked,1)
%                 jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
%             end
%             clear concat_input
%             smoothed_25 = jnk(:,half_win+1:end-half_win);
%             concat_input = cat(2,(percentile75(:,half_win:-1:1)),(percentile75(:,:)));
%             concat_input = cat(2,concat_input, (percentile75(:,end:-1:end-half_win+1)));
%             for k=1:size(evoked,1)
%                 jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
%             end
%             clear concat_input
%             smoothed_75 = jnk(:,half_win+1:end-half_win);
%             
%             sp=5;
%             sph(c,sp,t)=subplot(3,2,sp);
%             ax=sph(c,sp,t);
%             hold on;
%             plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
%             
%             %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
%             %plot(repmat(time,size(tar(t).con(c).evoked_avg,1),1)', squeeze(smoothed_mean)')
%             xlabel('Time(s)'); ylabel('sig. LFP [% of sites]');
%             title(plot_names{sp},'fontsize',10,'interpreter','none');
%             add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
%             
%             results_file{c,t} = fullfile(cfg.analyse_evoked_folder, [cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-','signed sig% of ',num2str(tar(t).nSites),' sites ',withunits]);
%             mtit(h(c,t),[ cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'- (signed) sig% of ',num2str(tar(t).nSites),' sites ' ,withunits],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
%             
%         end
%         
%     end
%     
%     for c = 1:length(cond)
%         for t=1:length(targets)
%         figure(h(c,t));
%         for sp=1:2
%             subplot(sph(c,sp,t));
%             set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
%         end
%         export_fig(h(c,t),[results_file{c,t},'.pdf']);
%         end
%     end
% end
% close all,
% 
% clear results_file sph h collim
% clc
% 
% for e=1:size(cfg.analyse_states,1)
%     E=cfg.analyse_states{e,1};
%     for t = 1:length(targets)
%         if tar(t).nSites == 0
%             continue;
%         end
%         h = figure('Name',['Max ITPCbp/POWbp in ' E ' of Target=' tar(t).target],'NumberTitle','off');
%         for c = 1:length(cond)
%             for fb = 1: length(freqb)
%                 % itpcbp
%                 subplot(2,2,2*c-1)
%                 scatter(tar(t).con(c).(E).max_itpcbp_time(fb,:),tar(t).con(c).(E).max_itpcbp(fb,:),15,colsbp(fb,:))
%                 hold on
%                 title([' max itpc-bp in ', strrep(cond{c},'_',' '),' for ',num2str(tar(t).nSites),...
%                     ' sites,',num2str(tar(t).con(c).(E).nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%                 xlabel('time (s)','Interpreter','latex'), ylabel('Max ITPC value','Interpreter','latex');
%                 
%                 % powbp
%                 subplot(2,2,2*c)
%                 scatter(tar(t).con(c).(E).max_powbp_time(fb,:),tar(t).con(c).(E).max_powbp(fb,:),15,colsbp(fb,:))
%                 hold on
%                 title([' max power-bp in ', strrep(cond{c},'_',' '),' for ',num2str(tar(t).nSites),...
%                     ' sites,',num2str(tar(t).con(c).concat.nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%                 xlabel('time (s)','Interpreter','latex'), ylabel('Max POWER value','Interpreter','latex');
%                 
%             end
%             
%             subplot(2,2,2*c-1)
%             plot([0,0],get(gca,'ylim'),'k--')
%             
%             subplot(2,2,2*c)
%             plot([0,0],get(gca,'ylim'),'k--')
%             legend(freqb,'FontSize',5)
%             legend('boxoff')
%             
%         end
%         mtit([strrep(tar(t).target,'_','-') '-' E],'xoff', 0, 'yoff', 0.05, 'fontsize', 12,'Interpreter', 'none');
%         results_file = fullfile(cfg.analyse_evoked_folder, [cfg.monkey,'-',targets{t},'-',E,'-Max_ITPCbp_POWbp_of ',num2str(tar(t).nSites),' sites ', withunits]);
%         export_fig(h,[results_file,'.pdf']);
%     end
%     close all,
%     clc    
% end
% 
% 
% %%
% % ploting the avgogram of the timing of the Max ITPC/POW:
% bins=cfg.analyse_states{1,4}:cfg.evoked.timestep*10:cfg.analyse_states{1,5};
% for t = 1: length(targets)
%     if tar(t).nSites == 0
%         continue;
%     end
%     h = figure;
%     for c = 1:length(cond)
%         for fb = 1: length(freqb)
%             % itpcbp
%             sp1=subplot(2,2,2*c-1);
%             [nelements,centers]= hist(squeeze(tar(t).con(c).concat.max_itpcbp_time(fb,:)),bins);
%             plot(centers,nelements,'-','Color',colsbp(fb,:));
%             hold on
%             title([' max itpc-bp in ', strrep(cond{c},'_',' '),' for ',num2str(size(tar(t).con(c).concat.max_itpcbp,2)),...
%                 ' sites,',num2str(tar(t).con(c).concat.nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('Max ITPC times','Interpreter','latex');
%             set(gca,'xlim',[-0.25,0.25])
%             
%             % powbp
%             sp2=subplot(2,2,2*c);
%             [nelements,centers]= hist(squeeze(tar(t).con(c).concat.max_powbp_time(fb,:)),bins);
%             plot(centers,nelements,'-','Color',colsbp(fb,:));
%             hold on
%             title([' max power-bp in ', strrep(cond{c},'_',' '),' for ',num2str(size(tar(t).con(c).concat.max_powbp,2)),...
%                 ' sites,',num2str(tar(t).con(c).concat.nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('Max POWER times','Interpreter','latex');
%             set(gca,'xlim',[-0.25,0.25])
%             
%         end
%         legend(freqb,'FontSize',5)
%         legend('boxoff')
%         
%         plot(sp1,[0,0],get(sp1,'ylim'),'k--')
%         plot(sp2,[0,0],get(sp2,'ylim'),'k--')
%     end
%     mtit(['Target=',strrep(tar(t).target,'_','-')],'xoff', 0, 'yoff', 0.05, 'Color','red', 'fontsize', 12);
%     results_file = fullfile(cfg.analyse_evoked_folder, [cfg.monkey,'-',targets{t},'-','Time_of_Max_ITPCbp_POWbp_for ',num2str(tar(t).nSites),' sites ', withunits]);
%     export_fig(h,[results_file,'.pdf']);
% %     saveas(h,[results_file]);
% end
% 
% close all,
% clc



end


function add_ticks_and_labels(ax,events,ylm,stp)

%ylm=get(ax,'ylim');
axes(ax)
event_onsets    =events.onset;
event_names     =events.name;
event_samples   =events.ticksamples;
state_ticks     =events.ticks;

for so = event_onsets
    line([so so], ylm, 'color', 'k');
    event_name = event_names{event_onsets == so};
    event_name=strrep(event_name,'_',' ');
    text(so+1, ylm(1)+(diff(ylim)/stp), event_name, 'fontsize', 8, 'fontweight', 'bold');
    %end
end

% mark state onsets
set(ax,'ylim',ylm)
set(ax,'xtick',event_samples(~isnan(state_ticks)))
set(ax,'xticklabels', state_ticks(~isnan(state_ticks)), 'fontsize', 6)
set(ax, 'xlim', [0 event_samples(end)] + 0.5);

% add 0.5 since the time value is the center of the bin
% add 0 at the beginning to make the y-axis visib
end

function out=smoothit(in,half_win,gaussian_kernel)
concat_input = cat(2,(in(:,half_win:-1:1)),in);
concat_input = cat(2,concat_input, (in(:,end:-1:end-half_win+1)));
for k=1:size(concat_input,1)
    jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');
end
out = jnk(:,half_win+1:end-half_win);
end
