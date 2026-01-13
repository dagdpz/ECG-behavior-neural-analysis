function tar = ecg_bna_compute_grand_avg_mua(cfg,withunits)

reprocess=1;
reprocess2=1;
fileName = fullfile([cfg.MUA_root_results_fldr filesep 'Grand_avg_mua_sites_',withunits,'.mat']);
fileName2 = fullfile([cfg.MUA_root_results_fldr filesep 'Grand_avg_mua_targets_',withunits,'.mat']);

if cfg.combine_hemispheres
    cfg.targets=unique(cellfun(@(x) x(1:strfind(x,'_')-1),cfg.targets,'uniformoutput',false));
end
targets = cfg.targets;%unique({sites.target}); %%cfg.monkey

if reprocess
    data_path = cfg.sites_mua_fldr;
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
                
                nTriggers = event.real.ntriggers;
                time      = event.time;
                mua_time  = event.tfr_time;
                mua       = squeeze(event.normalized.mua.mean)';
                mua_real       = squeeze(event.real.mua.mean)';
                mua_real_std       = squeeze(event.real.mua.std)';
                mua_shuffled       = squeeze(event.shuffled.mua.mean)';
                mua_shuffled_std       = squeeze(event.shuffled.mua.std)';
                p=event.prevspost.p;
                h=event.prevspost.h;
                
                mua_sig     =squeeze(event.significance.mua)'; %% weird dimension thing
                                
                sites(s).condition(c).event(e).nTriggers = nTriggers;
                sites(s).condition(c).event(e).mua = mua;               
                sites(s).condition(c).event(e).mua_sig = mua_sig;                  
                sites(s).condition(c).event(e).any_sig = any(mua_sig);   
                sites(s).condition(c).event(e).mua_real = mua_real;               
                sites(s).condition(c).event(e).mua_real_std = mua_real_std;               
                sites(s).condition(c).event(e).mua_shuffled = mua_shuffled;               
                sites(s).condition(c).event(e).mua_shuffled_std = mua_shuffled_std;        
                sites(s).condition(c).event(e).p = p;       
                sites(s).condition(c).event(e).h = h;                      
                
                sites(s).condition(c).event(e).time = time;
                sites(s).condition(c).event(e).tfr_time = mua_time;
                
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
%         concat.mua          = [];
%         concat.mua_25       = [];
%         concat.mua_75       = [];
%         concat.muasig       = [];
%         concat.mua_time     = [];
%         concat.mua_std      = [];
%         
%         for e=1:size(cfg.analyse_states,1)
%             E.nTriggers     = mean(cat(3,events(:,e).nTriggers),3);
%             E.mua           = mean(cat(3,events(:,e).mua),3);
%             E.mua_25        = prctile(cat(3,events(:,e).mua),25,3);
%             E.mua_75        = prctile(cat(3,events(:,e).mua),75,3);
%             E.mua_time      = events(1,e).time;
%             E.mua_std       = std(cat(3,events(:,e).mua),0,3);                        
%             E.mua_sig       = mean(cat(3,events(:,e).mua_sig),3);            
%             E.mua_sig_signed= mean(cat(3,events(:,e).mua_sig).*sign(cat(3,events(:,e).mua)),3);
%             
%             E.mua_sig       = E.mua_sig*100;                
%             E.mua_sig_signed= E.mua_sig_signed*100;     
%             
%             NaNseparator=100/25;
%             concat.mua          = cat(2, concat.mua,        E.mua,        nan(size(E.mua, 1),        NaNseparator));
%             concat.mua_25       = cat(2, concat.mua_25,     E.mua_25,	  nan(size(E.mua_25, 1),     NaNseparator));
%             concat.mua_75       = cat(2, concat.mua_75,     E.mua_75,	  nan(size(E.mua_75, 1),     NaNseparator));
%             
%             concat.muasig       = cat(2, concat.muasig,     E.mua_sig,    nan(size(E.mua_sig, 1),    NaNseparator));
%             concat.mua_std      = cat(2, concat.mua_std,    E.mua_std,	  nan(size(E.mua_std, 1),    NaNseparator));
%             concat.mua_time     = [concat.mua_time, E.mua_time, nan(1, NaNseparator)];
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
        cueresponse_sigbin=[taskevents(:,4).any_sig];
        cueresponse_ttest=[taskevents(:,4).p]<0.05;
        
        
        
        for c = 1:size(condition,2)
            events=cat(1,condition(:,c).event);
            
            concat.nTriggers=round(mean(reshape([events.nTriggers],size(events)),1));
            
            concat.mua          = [];
            concat.mua_abs      = [];
            concat.mua_abs_25   = [];
            concat.mua_abs_75   = [];
            concat.mua_25       = [];
            concat.mua_75       = [];
            concat.muasig       = [];
            concat.mua_time     = [];
            concat.mua_std      = [];
            concat.mua_abs_std      = [];
            
            for e=1:size(cfg.analyse_states,1)
                
                sig_sites=[taskevents(:,e).any_sig] | [restevents(:,e).any_sig];
                %S=sig_sites;
                S=true(size(sig_sites));
                %S=cueresponse_sigbin;
                %S=cueresponse_ttest;
                
                
                %             mua=zscore(cat(3,events(S,e).mua),0,2);
                %             mua_real=zscore(cat(3,events(S,e).mua_real),0,2);
                %             mua_shuffled=zscore(cat(3,events(S,e).mua_shuffled),0,2);
                mua=cat(3,events(S,e).mua);
                mua_real=cat(3,events(S,e).mua_real);
                mua_shuffled=cat(3,events(S,e).mua_shuffled);
                
                
                E.mua_p= ecg_bna_population_significant_n_bins(cat(3,events(:,e).mua_sig));
                E.nTriggers     = mean(cat(3,events(S,e).nTriggers),3);
                E.mua           = mean(mua,3);
                E.mua_real      = mean(mua_real,3);
                E.mua_shuffled  = mean(mua_shuffled,3);
                E.mua_abs       = mean(abs(mua),3);
                E.mua_abs_sterr = sterr(abs(mua),3);
                E.mua_abs_25    = prctile(abs(mua),25,3);
                E.mua_abs_75    = prctile(abs(mua),75,3);
                E.mua_25        = prctile(mua,25,3);
                E.mua_75        = prctile(mua,75,3);
                E.mua_time      = events(1,e).time;
                E.mua_std       = std(mua,0,3);
                E.mua_real_std  = std(mua_real,0,3);
                E.mua_sterr     = sterr(mua,3);
                E.mua_real_sterr= sterr(mua_real,3);
                E.mua_real_abs  = mean(abs(mua_real),3);
                E.mua_real_abs_std  = std(abs(mua_real),0,3);
                E.mua_real_abs_sterr  = sterr(abs(mua_real),3);
                E.mua_sig       = mean(cat(3,events(S,e).mua_sig),3);
                E.mua_abs_std   = std(abs(mua),0,3);
                E.mua_sig_signed= mean(cat(3,events(S,e).mua_sig).*sign(cat(3,events(S,e).mua)),3);
                [E.mua_popsig, E.mua_poppval]=ecg_bna_population_significance_normalized_parametric(mua);
                   
                
                E.mua_sig       = E.mua_sig*100;
                E.mua_sig_signed= E.mua_sig_signed*100;
                
                NaNseparator=100/25;
                concat.mua          = cat(2, concat.mua,        E.mua,        nan(size(E.mua, 1),        NaNseparator));
                concat.mua_25       = cat(2, concat.mua_25,     E.mua_25,	  nan(size(E.mua_25, 1),     NaNseparator));
                concat.mua_75       = cat(2, concat.mua_75,     E.mua_75,	  nan(size(E.mua_75, 1),     NaNseparator));
                concat.mua_abs      = cat(2, concat.mua_abs,    E.mua_abs,    nan(size(E.mua_abs, 1),    NaNseparator));
                concat.mua_abs_25   = cat(2, concat.mua_abs_25, E.mua_abs_25, nan(size(E.mua_abs_25, 1), NaNseparator));
                concat.mua_abs_75   = cat(2, concat.mua_abs_75, E.mua_abs_75, nan(size(E.mua_abs_75, 1), NaNseparator));
                
                concat.muasig       = cat(2, concat.muasig,     E.mua_sig,    nan(size(E.mua_sig, 1),    NaNseparator));
                concat.mua_std      = cat(2, concat.mua_std,    E.mua_std,	  nan(size(E.mua_std, 1),    NaNseparator));
                concat.mua_abs_std  = cat(2, concat.mua_abs_std,E.mua_abs_std,nan(size(E.mua_abs_std, 1),NaNseparator));
                concat.mua_time     = [concat.mua_time, E.mua_time, nan(1, NaNseparator)];
                
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
        mua_time=tar(t).con(c).concat.mua_time;
        tfr_events.onset        = find(mua_time == 0);
        tfr_events.name         = cfg.analyse_states(:,1);
        tfr_events.ticksamples  = sort([find(diff([0 ~isnan(mua_time)])==1) find(diff(~isnan(mua_time))==-1)]);
        tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
        tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
        tfr_events.ticks        = round(mua_time(tfr_events.ticksamples)*100)/100;
        
        mua=tar(t).con(c).concat.mua;
        mua_std=tar(t).con(c).concat.mua_std;
        %mua_se=sterr(mua,3,0);
        percentile25=tar(t).con(c).concat.mua_25;
        percentile75=tar(t).con(c).concat.mua_75;
        percentile25(isnan(percentile25))=0;
        percentile75(isnan(percentile75))=0;
        
        % create figure
        h(c,t) = figure('units','normalized','position',[0 0 1 1]);
        
        %========================== LFP evoked Potential ==================== %
        
        smoothed_mean=smoothit(mua,half_win,gaussian_kernel);
        smoothed_std=smoothit(mua_std,half_win,gaussian_kernel);
        smoothed_25=smoothit(percentile25,half_win,gaussian_kernel);
        smoothed_75=smoothit(percentile75,half_win,gaussian_kernel);
            
        sp=1;
        sph(c,sp,t)=subplot(2,1,sp);
        ax=sph(c,sp,t);
        hold on;
        lineProps={'color',[0 0 1]};
        %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,[smoothed_75-smoothed_mean;smoothed_mean-smoothed_25 ],lineProps,1);        
        shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
        %plot(repmat(time,size(tar(t).con(c).mua_avg,1),1)', squeeze(smoothed_mean)')
        
        xlabel('Time(s)'); ylabel('MUA evoked Potential');
        title(plot_names{sp},'fontsize',10,'interpreter','none');
        add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)        
        
        results_file{c} = fullfile(cfg.MUA_root_results_fldr, [cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-','avg of ',num2str(tar(t).nSites),' sites ',withunits]);
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
            
            mua_time                = tar(t).con(c).(E).mua_time;
            tfr_events.onset        = find(mua_time == 0);
            tfr_events.name         = cfg.analyse_states(e,1);
            tfr_events.ticksamples  = [1,numel(mua_time)]; 
            tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
            tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
            tfr_events.ticks        = round(mua_time(tfr_events.ticksamples)*100)/100;
            
            %========================== MUA evoked Potential ==================== %
            mua=tar(t).con(c).(E).mua;
            mua_abs=tar(t).con(c).(E).mua_abs;
            mua_std=tar(t).con(c).(E).mua_std;
            mua_sterr=tar(t).con(c).(E).mua_sterr;
            mua_abs_std=tar(t).con(c).(E).mua_abs_std;
            mua_abs_sterr=tar(t).con(c).(E).mua_abs_sterr;
            %mua_se=sterr(mua,3,0);
            percentile25=tar(t).con(c).(E).mua_25;
            percentile75=tar(t).con(c).(E).mua_75;
            percentile25(isnan(percentile25))=0;
            percentile75(isnan(percentile75))=0;
           
            % Smoothing of the  LFP evoked Potential here:
            smoothed_mean=smoothit(mua,half_win,gaussian_kernel);    
            smoothed_std=smoothit(mua_std,half_win,gaussian_kernel);
            smoothed_sterr=smoothit(mua_sterr,half_win,gaussian_kernel);
            smoothed_25=smoothit(percentile25,half_win,gaussian_kernel);
            smoothed_75=smoothit(percentile75,half_win,gaussian_kernel);
            
            %========================== mua  ==================== %
            sp=1;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            lineProps={'color',[0 0 1]};
            
            %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
            shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_sterr,lineProps,1);
            %plot(repmat(time,size(tar(t).con(c).mua_avg,1),1)', squeeze(smoothed_mean)')
            xlabel('Time(s)'); ylabel('MUA evoked Potential');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)            
            
            
            %========================== mua significance ==================== %
            mua_sig=tar(t).con(c).(E).mua_sig;
            smoothed_mean=smoothit(mua_sig,half_win,gaussian_kernel);          
            
            sp=2;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
            
            xlabel('Time(s)'); ylabel('sig. MUA [% of sites]');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)   
            
            
            %========================== mua abs  ==================== %
            smoothed_mean=smoothit(mua_abs,half_win,gaussian_kernel);    
            smoothed_std=smoothit(mua_abs_std,half_win,gaussian_kernel);
            smoothed_sterr=smoothit(mua_abs_sterr,half_win,gaussian_kernel);
            sp=3;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            lineProps={'color',[0 0 1]};
            
            %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
            shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_sterr,lineProps,1);
            %plot(repmat(time,size(tar(t).con(c).mua_avg,1),1)', squeeze(smoothed_mean)')
            xlabel('Time(s)'); ylabel('MUA Absolute');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)            
            
            
            %========================== mua significance ==================== %
            mua_sig=tar(t).con(c).(E).mua_sig_signed;
            smoothed_mean=smoothit(mua_sig,half_win,gaussian_kernel);          
            
            sp=4;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
            
            xlabel('Time(s)'); ylabel('signed sig. MUA [% of sites]');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)   
            
            
            
            results_file{c,t} = fullfile(cfg.MUA_root_results_fldr, [cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-','avg of ',num2str(tar(t).nSites),' sites ',withunits]);
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
            
            mua_time                = tar(t).con(c).(E).mua_time;
            tfr_events.onset        = find(mua_time == 0);
            tfr_events.name         = cfg.analyse_states(e,1);
            tfr_events.ticksamples  = [1,numel(mua_time)]; 
            tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
            tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
            tfr_events.ticks        = round(mua_time(tfr_events.ticksamples)*100)/100;
            
            %========================== MUA evoked Potential ==================== %
            mua=tar(t).con(c).(E).mua_real;
            mua_abs=tar(t).con(c).(E).mua_real_abs;
            mua_std=tar(t).con(c).(E).mua_real_std;
            mua_sterr=tar(t).con(c).(E).mua_real_sterr;
            mua_abs_std=tar(t).con(c).(E).mua_real_abs_std;
            mua_abs_sterr=tar(t).con(c).(E).mua_real_abs_sterr;
            %mua_se=sterr(mua,3,0);
            percentile25=tar(t).con(c).(E).mua_25;
            percentile75=tar(t).con(c).(E).mua_75;
            percentile25(isnan(percentile25))=0;
            percentile75(isnan(percentile75))=0;
           
            % Smoothing of the  LFP evoked Potential here:
            smoothed_mean=smoothit(mua,half_win,gaussian_kernel);    
            smoothed_std=smoothit(mua_std,half_win,gaussian_kernel);
            smoothed_sterr=smoothit(mua_sterr,half_win,gaussian_kernel);
            smoothed_25=smoothit(percentile25,half_win,gaussian_kernel);
            smoothed_75=smoothit(percentile75,half_win,gaussian_kernel);
            
            %========================== mua  ==================== %
            sp=1;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            lineProps={'color',[0 0 1]};
            
            %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
            shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_sterr,lineProps,1);
            %plot(repmat(time,size(tar(t).con(c).mua_avg,1),1)', squeeze(smoothed_mean)')
            xlabel('Time(s)'); ylabel('MUA evoked Potential');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)            
            
            
            %========================== mua significance ==================== %
            mua_sig=tar(t).con(c).(E).mua_sig;
            smoothed_mean=smoothit(mua_sig,half_win,gaussian_kernel);          
            
            sp=2;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
            
            xlabel('Time(s)'); ylabel('sig. MUA [% of sites]');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)   
            
            
            %========================== mua abs  ==================== %
            smoothed_mean=smoothit(mua_abs,half_win,gaussian_kernel);    
            smoothed_std=smoothit(mua_abs_std,half_win,gaussian_kernel);
            smoothed_sterr=smoothit(mua_abs_sterr,half_win,gaussian_kernel);
            sp=3;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            lineProps={'color',[0 0 1]};
            
            %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
            shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_sterr,lineProps,1);
            %plot(repmat(time,size(tar(t).con(c).mua_avg,1),1)', squeeze(smoothed_mean)')
            xlabel('Time(s)'); ylabel('MUA Absolute');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)            
            
            
            %========================== mua significance ==================== %
            mua_sig=tar(t).con(c).(E).mua_sig_signed;
            smoothed_mean=smoothit(mua_sig,half_win,gaussian_kernel);          
            
            sp=4;
            sph(c,sp,t)=subplot(2,2,sp);
            ax=sph(c,sp,t);
            hold on;
            plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
            
            xlabel('Time(s)'); ylabel('signed sig. MUA [% of sites]');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)   
            
            
            
            results_file{c,t} = fullfile(cfg.MUA_root_results_fldr, ['real_' cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-','avg of ',num2str(tar(t).nSites),' sites ',withunits]);
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
            
            mua_time                = tar(t).con(c).(E).mua_time;
            mua_p                   = tar(t).con(c).(E).mua_p;
            if mua_p<0.001
                ptoplot='<0.001';
            else
                ptoplot=num2str(round(mua_p*1000)/1000);
            end
            tfr_events.onset        = find(mua_time == 0);
            tfr_events.name         = cfg.analyse_states(e,1);
            tfr_events.ticksamples  = [1,numel(mua_time)]; 
            tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
            tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
            tfr_events.ticks        = round(mua_time(tfr_events.ticksamples)*100)/100;      
            
            
            %========================== mua significance ==================== %
            mua_sig=tar(t).con(c).(E).mua_sig;
            smoothed_mean=smoothit(mua_sig,half_win,gaussian_kernel);          
            
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
    results_file = fullfile(cfg.MUA_root_results_fldr, [cfg.monkey,'-',E,'-sigNbinstesting-',withunits]);
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
%             mua_time=tar(t).con(c).(E).tfr_time;
%             %mua_time=tar(t).con(c).concat.mua_time;
%             tfr_events.onset        = find(mua_time == 0);
%             tfr_events.name         = cfg.analyse_states(e,1);
%             tfr_events.ticksamples  = [1,numel(mua_time)]; %sort([find(diff([0 ~isnan(tfr_time)])==1) find(diff([~isnan(tfr_time)])==-1)]);
%             tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
%             tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
%             tfr_events.ticks        = round(mua_time(tfr_events.ticksamples)*100)/100;
%             
%             mua=tar(t).con(c).(E).mua_sig;
%             mua_std=tar(t).con(c).(E).mua_std;
%             %mua_se=sterr(mua,3,0);
%             percentile25=tar(t).con(c).(E).mua_25;
%             percentile75=tar(t).con(c).(E).mua_75;
%             percentile25(isnan(percentile25))=0;
%             percentile75(isnan(percentile75))=0;
%            
%             
%             
%             %========================== LFP evoked Potential ==================== %
%             % Smoothing of the  LFP evoked Potential here:
%             jnk = [];
%             concat_input = cat(2,(mua(:,half_win:-1:1)),(mua(:,:)));
%             concat_input = cat(2,concat_input, (mua(:,end:-1:end-half_win+1)));
%             for k=1:size(mua,1)
%                 jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
%             end
%             clear concat_input
%             smoothed_mean = jnk(:,half_win+1:end-half_win);
%             concat_input = cat(2,(mua_std(:,half_win:-1:1)),(mua_std(:,:)));
%             concat_input = cat(2,concat_input, (mua_std(:,end:-1:end-half_win+1)));
%             for k=1:size(mua,1)
%                 jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
%             end
%             clear concat_input
%             smoothed_std = jnk(:,half_win+1:end-half_win);
%             concat_input = cat(2,(percentile25(:,half_win:-1:1)),(percentile25(:,:)));
%             concat_input = cat(2,concat_input, (percentile25(:,end:-1:end-half_win+1)));
%             for k=1:size(mua,1)
%                 jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
%             end
%             clear concat_input
%             smoothed_25 = jnk(:,half_win+1:end-half_win);
%             concat_input = cat(2,(percentile75(:,half_win:-1:1)),(percentile75(:,:)));
%             concat_input = cat(2,concat_input, (percentile75(:,end:-1:end-half_win+1)));
%             for k=1:size(mua,1)
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
%             %plot(repmat(time,size(tar(t).con(c).mua_avg,1),1)', squeeze(smoothed_mean)')
%             xlabel('Time(s)'); ylabel('sig. LFP [% of sites]');
%             title(plot_names{sp},'fontsize',10,'interpreter','none');
%             add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
%             
%             
%             results_file{c,t} = fullfile(cfg.analyse_mua_folder, [cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-','sig% of ',num2str(tar(t).nSites),' sites ',withunits]);
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
%             mua_time=tar(t).con(c).(E).tfr_time;
%             %mua_time=tar(t).con(c).concat.mua_time;
%             tfr_events.onset        = find(mua_time == 0);
%             tfr_events.name         = cfg.analyse_states(e,1);
%             tfr_events.ticksamples  = [1,numel(mua_time)]; %sort([find(diff([0 ~isnan(tfr_time)])==1) find(diff([~isnan(tfr_time)])==-1)]);
%             tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
%             tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
%             tfr_events.ticks        = round(mua_time(tfr_events.ticksamples)*100)/100;
%             
%             powbp=tar(t).con(c).(E).powbp_sig_signed;
%             itpcbp=tar(t).con(c).(E).itpcbp_sig_signed;
%             mua=tar(t).con(c).(E).mua_sig_signed;
%             mua_std=tar(t).con(c).(E).mua_std;
%             %mua_se=sterr(mua,3,0);
%             percentile25=tar(t).con(c).(E).mua_25;
%             percentile75=tar(t).con(c).(E).mua_75;
%             percentile25(isnan(percentile25))=0;
%             percentile75(isnan(percentile75))=0;
%             
%             
%             %========================== LFP evoked Potential ==================== %
%             % Smoothing of the  LFP evoked Potential here:
%             jnk = [];
%             concat_input = cat(2,(mua(:,half_win:-1:1)),(mua(:,:)));
%             concat_input = cat(2,concat_input, (mua(:,end:-1:end-half_win+1)));
%             for k=1:size(mua,1)
%                 jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
%             end
%             clear concat_input
%             smoothed_mean = jnk(:,half_win+1:end-half_win);
%             concat_input = cat(2,(mua_std(:,half_win:-1:1)),(mua_std(:,:)));
%             concat_input = cat(2,concat_input, (mua_std(:,end:-1:end-half_win+1)));
%             for k=1:size(mua,1)
%                 jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
%             end
%             clear concat_input
%             smoothed_std = jnk(:,half_win+1:end-half_win);
%             concat_input = cat(2,(percentile25(:,half_win:-1:1)),(percentile25(:,:)));
%             concat_input = cat(2,concat_input, (percentile25(:,end:-1:end-half_win+1)));
%             for k=1:size(mua,1)
%                 jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
%             end
%             clear concat_input
%             smoothed_25 = jnk(:,half_win+1:end-half_win);
%             concat_input = cat(2,(percentile75(:,half_win:-1:1)),(percentile75(:,:)));
%             concat_input = cat(2,concat_input, (percentile75(:,end:-1:end-half_win+1)));
%             for k=1:size(mua,1)
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
%             %plot(repmat(time,size(tar(t).con(c).mua_avg,1),1)', squeeze(smoothed_mean)')
%             xlabel('Time(s)'); ylabel('sig. LFP [% of sites]');
%             title(plot_names{sp},'fontsize',10,'interpreter','none');
%             add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
%             
%             results_file{c,t} = fullfile(cfg.analyse_mua_folder, [cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-','signed sig% of ',num2str(tar(t).nSites),' sites ',withunits]);
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
%         results_file = fullfile(cfg.analyse_mua_folder, [cfg.monkey,'-',targets{t},'-',E,'-Max_ITPCbp_POWbp_of ',num2str(tar(t).nSites),' sites ', withunits]);
%         export_fig(h,[results_file,'.pdf']);
%     end
%     close all,
%     clc    
% end
% 
% 
% %%
% % ploting the avgogram of the timing of the Max ITPC/POW:
% bins=cfg.analyse_states{1,4}:cfg.mua.timestep*10:cfg.analyse_states{1,5};
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
%     results_file = fullfile(cfg.analyse_mua_folder, [cfg.monkey,'-',targets{t},'-','Time_of_Max_ITPCbp_POWbp_for ',num2str(tar(t).nSites),' sites ', withunits]);
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
