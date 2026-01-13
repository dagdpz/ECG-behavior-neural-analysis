function tar = ecg_bna_compute_grand_avg(cfg,withunits)

reprocess=1;
reprocess2=0;
plotit=0;
fileName = fullfile([cfg.analyse_lfp_folder filesep cfg.monkey,'_',cfg.analyse_states{1, 2} ,'_Triggered_target_wise_Grand_grand_avg_sessions_sites',withunits,'.mat']);

if cfg.combine_hemispheres
    cfg.targets=unique(cellfun(@(x) x(1:strfind(x,'_')-1),cfg.targets,'uniformoutput',false));
end
targets = cfg.targets;%unique({sites.target}); %%cfg.monkey

if reprocess
    data_path = cfg.sites_lfp_fldr;
    %     data_path = 'Y:\Projects\Pulv_bodysignal\LFP\ECG_Magnus_TaskRest_test_generalTrig_Reref_shamim\sample session check';
    cd(data_path)
    % read the files
    all_lfp_data = dir('*.mat');
    % storing all the Triggered parameters of the sites in the results folder:
    sites = struct;
    s=0;
    for f = 1:length(all_lfp_data)
        load(all_lfp_data(f).name)
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
            if isempty(con.event) %|| con.event(1).real.ntriggers<1000 %% very temporary condition
                sites(s).all_conditions_present = 0;
                continue;
            end
            sites(s).condition(c).condition_name = con.label;
            for e=1:size(cfg.analyse_states,1)
                event=con.event(e);
                % think about event present?
                
                nTriggers = event.real.ntriggers;
                time      = event.time;
                tfr_time  = event.tfr_time;
                %                 lfp       = squeeze(event.real.lfp.mean)';
                %                 itpc      = squeeze(event.real.itpc.mean)-squeeze(event.shuffled.itpc.mean);
                %                 power     = squeeze(event.normalized.pow.mean);
                %                 itpcbp    = squeeze(event.real.itpcbp.mean)-squeeze(event.shuffled.itpcbp.mean);
                %                 powerbp   = squeeze(event.normalized.powbp.mean);
                
                
                lfp       = squeeze(event.normalized.lfp.mean)';
                itpc      = squeeze(event.normalized.itpc.mean);
                power     = squeeze(event.normalized.pow.mean);
                itpcbp    = squeeze(event.normalized.itpcbp.mean);
                powerbp   = squeeze(event.normalized.powbp.mean);
                
                
                lfpR       = squeeze(event.real.lfp.mean)';
                itpcR      = squeeze(event.real.itpc.mean);
                powerR     = squeeze(event.real.pow.mean);
                itpcbpR    = squeeze(event.real.itpcbp.mean);
                powerbpR   = squeeze(event.real.powbp.mean);
                
                
                lfpS       = squeeze(event.shuffled.lfp.mean)';
                itpcS      = squeeze(event.shuffled.itpc.mean);
                powerS     = squeeze(event.shuffled.pow.mean);
                itpcbpS    = squeeze(event.shuffled.itpcbp.mean);
                powerbpS   = squeeze(event.shuffled.powbp.mean);
                
                
                add_index=[0:size(itpcbp,1)-1]'*size(itpcbp,2);
                [~,max_itpcbp_idx] = max(abs(itpcbp),[],2);
                max_itpcbp_time=time(max_itpcbp_idx)';
                tmp=itpcbp';
                max_itpcbp=tmp(max_itpcbp_idx+add_index);
                [~,max_powbp_idx] = max(abs(powerbp),[],2);
                tmp=powerbp';
                max_powbp=tmp(max_powbp_idx+add_index);
                max_powbp_time=time(max_powbp_idx)';
                
                lfp_sig     =squeeze(event.significance.lfp)'; %% weird dimension thing
                itpc_sig    =squeeze(event.significance.itpc);
                pow_sig     =squeeze(event.significance.pow);
                itpcbp_sig  =squeeze(event.significance.itpcbp);
                powbp_sig   =squeeze(event.significance.powbp);
                
                sites(s).condition(c).event(e).nTriggers = nTriggers;
                sites(s).condition(c).event(e).lfp = lfp;
                sites(s).condition(c).event(e).itpc = itpc;
                sites(s).condition(c).event(e).pow = power;
                sites(s).condition(c).event(e).itpcbp = itpcbp;
                sites(s).condition(c).event(e).powbp = powerbp;
                
                sites(s).condition(c).event(e).lfp_R = lfpR;
                sites(s).condition(c).event(e).itpc_R = itpcR;
                sites(s).condition(c).event(e).pow_R = powerR;
                sites(s).condition(c).event(e).itpcbp_R = itpcbpR;
                sites(s).condition(c).event(e).powbp_R = powerbpR;
                
                sites(s).condition(c).event(e).lfp_S = lfpS;
                sites(s).condition(c).event(e).itpc_S = itpcS;
                sites(s).condition(c).event(e).pow_S = powerS;
                sites(s).condition(c).event(e).itpcbp_S = itpcbpS;
                sites(s).condition(c).event(e).powbp_S = powerbpS;
                
                sites(s).condition(c).event(e).max_itpcbp_time = max_itpcbp_time;
                sites(s).condition(c).event(e).max_powbp_time = max_powbp_time;
                sites(s).condition(c).event(e).max_itpcbp = max_itpcbp;
                sites(s).condition(c).event(e).max_powbp = max_powbp;
                sites(s).condition(c).event(e).lfp_sig = lfp_sig;
                sites(s).condition(c).event(e).itpc_sig = itpc_sig;
                sites(s).condition(c).event(e).pow_sig = pow_sig;
                
                sites(s).condition(c).event(e).itpcbp_sig = itpcbp_sig;
                sites(s).condition(c).event(e).powbp_sig = powbp_sig;
                sites(s).condition(c).event(e).time = time;
                sites(s).condition(c).event(e).tfr_time = tfr_time;
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

fileName = fullfile([cfg.analyse_lfp_folder filesep cfg.monkey,'_',cfg.analyse_states{1, 2} ,'_Triggered_target_wise_Grand_grand_avg_sessions_targets',withunits,'.mat']);

if reprocess2
    for t = 1: length(targets)
        sites_for_this_target=arrayfun(@(x) any(strfind(x.target,targets{t})),sites);
        target_sites = sites(sites_for_this_target);
        condition=vertcat(target_sites.condition);
        
        for c = 1:size(condition,2)
            events=cat(1,condition(:,c).event);
            n_freqs=size(events(1).max_itpcbp,1);
            
            concat.nTriggers=round(mean(reshape([events.nTriggers],size(events)),1));
            concat.max_itpcbp=reshape([events.max_itpcbp],[n_freqs,size(events)]);
            concat.max_powbp=reshape([events.max_powbp],[n_freqs,size(events)]);
            concat.max_itpcbp_time=reshape([events.max_itpcbp_time],[n_freqs,size(events)]);
            concat.max_powbp_time=reshape([events.max_powbp_time],[n_freqs,size(events)]);
            
            concat.pow          = [];
            concat.itpc         = [];
            concat.lfp          = [];
            concat.lfp_25       = [];
            concat.lfp_75       = [];
            concat.itpcbp       = [];
            concat.powbp        = [];
            concat.powsig       = [];
            concat.itpcsig      = [];
            concat.powbpsig     = [];
            concat.itpcbpsig    = [];
            concat.lfpsig       = [];
            concat.tfr_time     = [];
            concat.lfp_time     = [];
            concat.lfp_std      = [];
            
            for e=1:size(cfg.analyse_states,1)
                to_cat={'lfp','pow','itpc','powbp','itpcbp',};
                for n_cat=1:numel(to_cat)
                    CN=to_cat{n_cat};
                    switch CN
                        case {'pow','itpc'}
                            fstart=fx;
                        case {'powbp','itpcbp'}
                            fstart=fbx;
                        case {'lfp'}
                            fstart=1;
                            
                    end
                    CAT=cat(3,events(:,e).(CN));
                    CATsig=cat(3,events(:,e).([CN '_sig']));
                    CATsigsigned=cat(3,events(:,e).([CN '_sig'])).*sign(CAT);
                    CAT_R=cat(3,events(:,e).([CN '_R']));
                    CAT_S=cat(3,events(:,e).([CN '_S']));
                    
                    E.(CN)=mean(CAT(fstart:end,:,:),3);
                    E.([CN '_sig'])=mean(CATsig(fstart:end,:,:),3)*100;
                    E.([CN '_sig_signed'])=mean(CATsigsigned(fstart:end,:,:),3)*100;
                    %if ismember(CN,{'pow','itpc'})
                    %[E.([CN '_popsig']), E.([CN '_poppval'])]=ecg_bna_population_significance_nonparametric(CAT_R(fstart:end,:,:),CAT_S(fstart:end,:,:));
                    [E.([CN '_popsig']), E.([CN '_poppval'])]=ecg_bna_population_significance_normalized_parametric(CAT(fstart:end,:,:));
                   
                    %end
                    if ismember(CN,{'powbp','itpcbp'})
                        E.([CN '_p'])= ecg_bna_population_significant_n_bins(CATsig(fstart:end,:,:));
                        tmp_max=[events(:,e).(['max_' CN])];
                        tmp_time=[events(:,e).(['max_' CN '_time'])];
                        E.(['max_' CN])=    tmp_max(fstart:end,:)  ;
                        E.(['max_' CN '_time'])=    tmp_time(fstart:end,:);
                    end
                end
                
                E.itpcbp_p= ecg_bna_population_significant_n_bins(cat(3,events(:,e).itpcbp_sig));
                E.powbp_p= ecg_bna_population_significant_n_bins(cat(3,events(:,e).powbp_sig));
                E.lfp_p= ecg_bna_population_significant_n_bins(cat(3,events(:,e).lfp_sig));
                
                E.nTriggers = mean(cat(3,events(:,e).nTriggers),3);
                E.lfp_25    = prctile(cat(3,events(:,e).lfp),25,3);
                E.lfp_75    = prctile(cat(3,events(:,e).lfp),75,3);
                E.tfr_time    = events(1,e).tfr_time;
                E.lfp_time    = events(1,e).time;
                E.lfp_std     = std(cat(3,events(:,e).lfp),0,3);
                E.lfp_sterr   = sterr(cat(3,events(:,e).lfp),3);
                
                NaNseparator=100/25;
                concat.pow          = cat(2, concat.pow,        E.pow,        nan(size(E.pow, 1),        NaNseparator));
                concat.itpc         = cat(2, concat.itpc,       E.itpc,       nan(size(E.itpc,1),        NaNseparator));
                concat.lfp          = cat(2, concat.lfp,        E.lfp,        nan(size(E.lfp, 1),        NaNseparator));
                concat.itpcbp       = cat(2, concat.itpcbp,     E.itpcbp,     nan(size(E.itpcbp, 1),     NaNseparator));
                concat.lfp_25       = cat(2, concat.lfp_25,     E.lfp_25,	  nan(size(E.lfp_25, 1),     NaNseparator));
                concat.lfp_75       = cat(2, concat.lfp_75,     E.lfp_75,	  nan(size(E.lfp_75, 1),     NaNseparator));
                
                concat.powbp        = cat(2, concat.powbp,      E.powbp,      nan(size(E.powbp, 1),      NaNseparator));
                concat.powsig       = cat(2, concat.powsig,     E.pow_sig,    nan(size(E.pow_sig, 1),    NaNseparator));
                concat.itpcsig      = cat(2, concat.itpcsig,    E.itpc_sig,   nan(size(E.itpc_sig, 1),   NaNseparator));
                concat.powbpsig     = cat(2, concat.powbpsig,   E.powbp_sig,  nan(size(E.powbp_sig, 1),  NaNseparator));
                concat.itpcbpsig    = cat(2, concat.itpcbpsig,  E.itpcbp_sig, nan(size(E.itpcbp_sig, 1), NaNseparator));
                concat.lfpsig       = cat(2, concat.lfpsig,     E.lfp_sig,    nan(size(E.lfp_sig, 1),    NaNseparator));
                concat.lfp_std      = cat(2, concat.lfp_std,    E.lfp_std,	  nan(size(E.lfp_std, 1),    NaNseparator));
                concat.tfr_time     = [concat.tfr_time, E.tfr_time, nan(1, NaNseparator)];
                concat.lfp_time     = [concat.lfp_time, E.lfp_time, nan(1, NaNseparator)];
                
                eventname=cfg.analyse_states{e,1};
                tar(t).con(c).(eventname)=E;
            end
            tar(t).con(c).cond_name = cfg.condition(c).name;
            tar(t).con(c).concat=concat;
        end
        tar(t).target = targets{t};
        tar(t).nSites = length(target_sites);
    end
    save(fileName,'tar')
else
    load(fileName);
end

if plotit
    %% plotting the results:
    cond = {cfg.condition.name};
    plot_names={'POW','ITPC','Power_BP','ITPC_BP','LFP_Evoked'};
    mtitsettings={'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none'};
    
    % Smoothing Kernel here:
    win = 1:cfg.lfp.smoothWin; win=win-(numel(win)+1)/2;
    half_win = ceil(size(win,2)/2)-1;
    gaussian_kernel=normpdf(win,0,numel(win)/6);
    gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);
    
    %% figure 1 - overview
    for t = 1: length(targets)
        if tar(t).nSites == 0
            continue;
        end
        for c = 1:length(cond)
            tfr_time=tar(t).con(c).concat.tfr_time;
            tfr_events=get_ticks_and_labels(tfr_time,cfg.analyse_states(:,1));
            
            powbp   =tar(t).con(c).concat.powbp;
            itpcbp  =tar(t).con(c).concat.itpcbp;
            lfp     =tar(t).con(c).concat.lfp;
            lfp_std =tar(t).con(c).concat.lfp_std;
            %lfp_se=sterr(lfp,3,0);
            percentile25=tar(t).con(c).concat.lfp_25;
            percentile75=tar(t).con(c).concat.lfp_75;
            percentile25(isnan(percentile25))=0;
            percentile75(isnan(percentile75))=0;
            
            % create figure
            h(c,t) = figure('units','normalized','position',[0 0 1 1]);
            toplot={squeeze(tar(t).con(c).concat.pow),squeeze(tar(t).con(c).concat.itpc)};
            for sp=1:2
                % =========================== Power & ITPC ============================= %
                sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
                image(toplot{sp},'CDataMapping','scaled');
                set(ax,'YDir','normal');
                line([0.5 0.5], ylim, 'color', 'k');
                
                % horizontal lines to separate frequency bands
                fbandstart = unique(frequency_bands(:))';
                fbandstart_idx = zeros(size(fbandstart));
                for f = fbandstart
                    f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                    line([0.5 size(toplot{sp},2)], [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                    fbandstart_idx(fbandstart == f) = f_idx;
                end
                set(ax,'TickDir','out')
                set(ax, 'ytick', fbandstart_idx);
                set(ax, 'yticklabel', fbandstart);
                set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                title(plot_names{sp},'fontsize',10,'interpreter','none');
                nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
                collim{c,sp}=[min(nonnan(:)) max(nonnan(:))];
                add_ticks_and_labels(ax,tfr_events,[0.5,numel(freq) + 0.5],8)
            end
            
            %========================== Bandpassed POWER ==================== %
            sp=3;
            sph(c,sp,t)=subplot(3,2,sp);
            ax=sph(c,sp,t);
            hold on;
            set(ax,'ColorOrder',colsbp);
            smoothed=smoothit(powbp,half_win,gaussian_kernel);
            plot(squeeze(smoothed)')
            xlabel('Time(s)'); ylabel('Power (W)');
            legend(freqb,'fontsize',3);
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
            
            %========================== Bandpassed ITPC ==================== %
            sp=4;
            sph(c,sp,t)=subplot(3,2,sp);
            ax=sph(c,sp,t);
            hold on;
            set(ax,'ColorOrder',colsbp);
            smoothed=smoothit(itpcbp,half_win,gaussian_kernel);
            plot(squeeze(smoothed)')
            xlabel('Time(s)'); ylabel('ITPC');
            legend(freqb,'fontsize',3);
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
            
            %========================== LFP evoked Potential ==================== %
            sp=5;
            sph(c,sp,t)=subplot(3,2,sp);
            ax=sph(c,sp,t);
            hold on;
            lineProps={'color',[0 0 1]};
            smoothed_mean=smoothit(lfp,half_win,gaussian_kernel);
            smoothed_std=smoothit(lfp_std,half_win,gaussian_kernel);
            shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
            % smoothed_25=smoothit(percentile25,half_win,gaussian_kernel);
            % smoothed_75=smoothit(percentile75,half_win,gaussian_kernel);
            % shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,[smoothed_75-smoothed_mean;smoothed_mean-smoothed_25 ],lineProps,1);
            
            xlabel('Time(s)'); ylabel(' LFP evoked Potential');
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
            
            %=================================================================%
            %% format spectra colors
            cbtitle = {'(P - \mu) / std','P - \mu'};
            for sp=1:2
                subplot(3,2,sp);
                cm = colormap('jet');
                cb = colorbar;
                set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
                colormap(cm);
            end
            
            results_file{c} = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-','avg of ',num2str(tar(t).nSites),' sites ',withunits]);
            mtit(h(c,t),[ cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-avg of ',num2str(tar(t).nSites),' sites ' ,withunits],mtitsettings{:})
            
        end
        
        for c = 1:length(cond)
            figure(h(c,t));
            for sp=1:2
                subplot(sph(c,sp,t));
                set(sph(c,sp,t),'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
            end
            export_fig(h(c,t),[results_file{c},'.pdf']);
        end
    end
    clearall
    
    
    %% figure 2: same as figure 1, but separately for each event
    for e=1:size(cfg.analyse_states,1)
        E=cfg.analyse_states{e,1};
        
        for t = 1: length(targets)
            for c = 1:length(cond)
                k=t+(c-1)*length(targets);
                % create figure
                h(c,t) = figure('units','normalized','position',[0 0 1 1]);
                if tar(t).nSites == 0
                    continue;
                end
                
                tfr_time=tar(t).con(c).(E).tfr_time;
                tfr_events=get_ticks_and_labels(tfr_time,cfg.analyse_states(e,1));
                
                powbp=tar(t).con(c).(E).powbp;
                itpcbp=tar(t).con(c).(E).itpcbp;
                lfp=tar(t).con(c).(E).lfp;
                lfp_std=tar(t).con(c).(E).lfp_std;
                %lfp_se=sterr(lfp,3,0);
                percentile25=tar(t).con(c).(E).lfp_25;
                percentile75=tar(t).con(c).(E).lfp_75;
                percentile25(isnan(percentile25))=0;
                percentile75(isnan(percentile75))=0;
                
                toplot={tar(t).con(c).(E).pow,tar(t).con(c).(E).itpc};
                sigplot={tar(t).con(c).(E).pow_popsig,tar(t).con(c).(E).itpc_popsig};
                for sp=1:2
                    % =========================== Power & ITPC ============================= %
                    sph(c,sp,t)=subplot(3,2,sp);
                    ax=sph(c,sp,t);
                    hold on
                    image(toplot{sp},'CDataMapping','scaled');
                    set(gca,'YDir','normal');
                    
                    significance = double(squeeze(sigplot{sp}));
                    contour(1:size(toplot{sp},2),1:size(toplot{sp},1),significance,1,'linecolor','k')
                    
                    line([0.5 0.5], ylim, 'color', 'k');
                    
                    % horizontal lines to separate frequency bands
                    fbandstart = unique(frequency_bands(:))';
                    fbandstart_idx = zeros(size(fbandstart));
                    for f = fbandstart
                        f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                        line([0.5 size(toplot{sp},2)], [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                        fbandstart_idx(fbandstart == f) = f_idx;
                    end
                    set(ax,'TickDir','out')
                    set(ax, 'ytick', fbandstart_idx);
                    set(ax, 'yticklabel', fbandstart);
                    set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                    set(ax, 'xlim', [0.5 tfr_events.ticksamples(end)] + 0.5);
                    title(plot_names{sp},'fontsize',10,'interpreter','none');
                    nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
                    collim{k,sp}=[min(nonnan(:)) max(nonnan(:))];
                    add_ticks_and_labels(ax,tfr_events,[0.5,numel(freq) + 0.5],8)
                end
                
                %========================== Bandpassed POWER ==================== %
                sp=3;
                sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
                hold on;
                set(ax,'ColorOrder',colsbp);
                smoothed=smoothit(powbp,half_win,gaussian_kernel);
                plot(squeeze(smoothed)')
                xlabel('Time(s)'); ylabel('Power (W)');
                legend(freqb,'fontsize',3);
                title(plot_names{sp},'fontsize',10,'interpreter','none');
                add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
                
                %========================== Bandpassed ITPC ==================== %
                sp=4;
                sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
                hold on;
                set(ax,'ColorOrder',colsbp);
                smoothed=smoothit(itpcbp,half_win,gaussian_kernel);
                plot(squeeze(smoothed)')
                xlabel('Time(s)'); ylabel('ITPC');
                legend(freqb,'fontsize',3);
                title(plot_names{sp},'fontsize',10,'interpreter','none');
                add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
                
                %========================== LFP evoked Potential ==================== %
                sp=5;
                sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
                hold on;
                lineProps={'color',[0 0 1]};
                smoothed_mean=smoothit(lfp,half_win,gaussian_kernel);
                smoothed_std=smoothit(lfp_std,half_win,gaussian_kernel);
                shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
                %smoothed_25=smoothit(percentile25,half_win,gaussian_kernel);
                %smoothed_75=smoothit(percentile75,half_win,gaussian_kernel);
                %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,[smoothed_75-smoothed_mean;smoothed_mean-smoothed_25 ],lineProps,1);
                %plot(repmat(time,size(tar(t).con(c).lfp_avg,1),1)', squeeze(smoothed_mean)')
                xlabel('Time(s)'); ylabel(' LFP evoked Potential');
                title(plot_names{sp},'fontsize',10,'interpreter','none');
                add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
                
                %             % format spectra colors
                %             cbtitle = {'(P - \mu) / std','P - \mu'};
                %             for sp=1:2
                %                 subplot(3,2,sp);
                %                 cm = colormap('jet');
                %                 cb = colorbar;
                %                 set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                %                 set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
                %                 colormap(cm);
                %             end
                
                results_file{c,t} = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-','avg of ',num2str(tar(t).nSites),' sites ',withunits]);
                mtit(h(c,t),[ cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-avg of ',num2str(tar(t).nSites),' sites ' ,withunits],mtitsettings{:})
                
            end
            
        end
        
        for c = 1:length(cond)
            for t=1:length(targets)
                k=[t,t+numel(targets)]; %same scaling for conditions within target
                %k=[t+(c-1)*numel(targets)]; % individual scaling
                figure(h(c,t));
                for sp=1:2
                    subplot(sph(c,sp,t));
                    % format spectra colors
                    cbtitle = {'(P - \mu) / std','P - \mu'};
                    %cm = colormap('jet');
                    if sp==1 %power
                        cm=jet(128);
                        colormap(sph(c,sp,t),cm);
                        cb = colorbar;
                        set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                        set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
                        %set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
                        M=max(abs([min([collim{k,sp}]),max([collim{k,sp}])]));
                        set(gca,'CLim',[-M M]);
                    elseif sp==2 % itpc
                        cm=jet(255);
                        cm=cm(128:end,:);
                        colormap(sph(c,sp,t),cm);
                        cb = colorbar;
                        set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                        set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
                        %set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
                        M=max(abs([min([collim{k,sp}]),max([collim{k,sp}])]));
                        set(gca,'CLim',[0 M]);
                    end
                end
                export_fig(h(c,t),[results_file{c,t},'.pdf']);
            end
        end
    end
    clearall
    
    %% figure 2.2: all power together
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
                
                tfr_time=tar(t).con(c).(E).tfr_time;
                tfr_events=get_ticks_and_labels(tfr_time,cfg.analyse_states(e,1));
                
                powbp=tar(t).con(c).(E).powbp;
                itpcbp=tar(t).con(c).(E).itpcbp;
                lfp=tar(t).con(c).(E).lfp;
                lfp_std=tar(t).con(c).(E).lfp_std;
                %lfp_se=sterr(lfp,3,0);
                percentile25=tar(t).con(c).(E).lfp_25;
                percentile75=tar(t).con(c).(E).lfp_75;
                percentile25(isnan(percentile25))=0;
                percentile75(isnan(percentile75))=0;
                
                toplot={tar(t).con(c).(E).pow,tar(t).con(c).(E).itpc};
                sigplot={tar(t).con(c).(E).pow_popsig,tar(t).con(c).(E).itpc_popsig};
                sp=1;
                % =========================== Power & ITPC ============================= %
                sph(c,sp,t)=subplot(3,2,k);
                ax=sph(c,sp,t);
                hold on
                image(toplot{sp},'CDataMapping','scaled');
                set(gca,'YDir','normal');
                
                significance = double(squeeze(sigplot{sp}));
                contour(1:size(toplot{sp},2),1:size(toplot{sp},1),significance,1,'linecolor','k')
                
                line([0.5 0.5], ylim, 'color', 'k');
                
                % horizontal lines to separate frequency bands
                fbandstart = unique(frequency_bands(:))';
                fbandstart_idx = zeros(size(fbandstart));
                for f = fbandstart
                    f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                    line([0.5 size(toplot{sp},2)], [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                    fbandstart_idx(fbandstart == f) = f_idx;
                end
                set(ax,'TickDir','out')
                set(ax, 'ytick', fbandstart_idx);
                set(ax, 'yticklabel', fbandstart);
                set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                set(ax, 'xlim', [0.5 tfr_events.ticksamples(end)] + 0.5);
                title([targets{t},'-',tar(t).con(c).cond_name '-' num2str(tar(t).nSites) 'sites'],'fontsize',10,'interpreter','none');
                nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
                collim{k,sp}=[min(nonnan(:)) max(nonnan(:))];
                add_ticks_and_labels(ax,tfr_events,[0.5,numel(freq) + 0.5],8)
                
            end
            
        end
        
        k=0;
        for t=1:length(targets)
            for c = 1:length(cond)
                k=(t-1)*length(cond)+[1 2]; %same scaling for conditions within target
                %k=[t+(c-1)*numel(targets)]; % individual scaling
                sp=1;
                subplot(sph(c,sp,t));
                % format spectra colors
                cbtitle = {'(P - \mu) / std','P - \mu'};
                %cm = colormap('jet');
                if sp==1 %power
                    cm=jet(128);
                    colormap(sph(c,sp,t),cm);
                    cb = colorbar;
                    set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                    set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
                    %set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
                    M=max(abs([min([collim{k,sp}]),max([collim{k,sp}])]));
                    set(gca,'CLim',[-M M]);
                elseif sp==2 % itpc
                    cm=jet(255);
                    cm=cm(128:end,:);
                    colormap(sph(c,sp,t),cm);
                    cb = colorbar;
                    set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                    set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
                    %set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
                    M=max(abs([min([collim{k,sp}]),max([collim{k,sp}])]));
                    set(gca,'CLim',[0 M]);
                end
            end
        end
        
        results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',E,'-power,' withunits]);
        mtit(h(e),[ cfg.monkey,'-',E,'-power,' withunits],mtitsettings{:})
        export_fig(h(e),[results_file,'.pdf']);
    end
    
    %% figure 2.3: all itpc together
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
                
                tfr_time=tar(t).con(c).(E).tfr_time;
                tfr_events=get_ticks_and_labels(tfr_time,cfg.analyse_states(e,1));
                
                powbp=tar(t).con(c).(E).powbp;
                itpcbp=tar(t).con(c).(E).itpcbp;
                lfp=tar(t).con(c).(E).lfp;
                lfp_std=tar(t).con(c).(E).lfp_std;
                %lfp_se=sterr(lfp,3,0);
                percentile25=tar(t).con(c).(E).lfp_25;
                percentile75=tar(t).con(c).(E).lfp_75;
                percentile25(isnan(percentile25))=0;
                percentile75(isnan(percentile75))=0;
                
                toplot={tar(t).con(c).(E).pow,tar(t).con(c).(E).itpc};
                sigplot={tar(t).con(c).(E).pow_popsig,tar(t).con(c).(E).itpc_popsig};
                sp=2;
                % =========================== Power & ITPC ============================= %
                sph(c,sp,t)=subplot(3,2,k);
                ax=sph(c,sp,t);
                hold on
                image(toplot{sp},'CDataMapping','scaled');
                set(gca,'YDir','normal');
                
                significance = double(squeeze(sigplot{sp}));
                contour(1:size(toplot{sp},2),1:size(toplot{sp},1),significance,1,'linecolor','k')
                
                line([0.5 0.5], ylim, 'color', 'k');
                
                % horizontal lines to separate frequency bands
                fbandstart = unique(frequency_bands(:))';
                fbandstart_idx = zeros(size(fbandstart));
                for f = fbandstart
                    f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                    line([0.5 size(toplot{sp},2)], [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                    fbandstart_idx(fbandstart == f) = f_idx;
                end
                set(ax,'TickDir','out')
                set(ax, 'ytick', fbandstart_idx);
                set(ax, 'yticklabel', fbandstart);
                set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                set(ax, 'xlim', [0.5 tfr_events.ticksamples(end)] + 0.5);
                title([targets{t},'-',tar(t).con(c).cond_name '-' num2str(tar(t).nSites) 'sites'],'fontsize',10,'interpreter','none');
                nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
                collim{k,sp}=[min(nonnan(:)) max(nonnan(:))];
                add_ticks_and_labels(ax,tfr_events,[0.5,numel(freq) + 0.5],8)
                
            end
            
        end
        
        k=0;
        
        for t=1:length(targets)
            for c = 1:length(cond)
                k=(t-1)*length(cond)+[1 2]; %same scaling for conditions within target
                %k=[t+(c-1)*numel(targets)]; % individual scaling
                sp=2;
                subplot(sph(c,sp,t));
                % format spectra colors
                cbtitle = {'(P - \mu) / std','P - \mu'};
                %cm = colormap('jet');
                if sp==1 %power
                    cm=jet(128);
                    colormap(sph(c,sp,t),cm);
                    cb = colorbar;
                    set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                    set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
                    %set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
                    M=max(abs([min([collim{k,sp}]),max([collim{k,sp}])]));
                    set(gca,'CLim',[-M M]);
                elseif sp==2 % itpc
                    cm=jet(255);
                    cm=cm(128:end,:);
                    colormap(sph(c,sp,t),cm);
                    cb = colorbar;
                    set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                    set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
                    %set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
                    M=max(abs([min([collim{k,sp}]),max([collim{k,sp}])]));
                    set(gca,'CLim',[0 M]);
                end
            end
        end
        
        results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',E,'-itpc,' withunits]);
        mtit(h(e),[ cfg.monkey,'-',E,'-itpc,' withunits],mtitsettings{:})
        export_fig(h(e),[results_file,'.pdf']);
    end
    
    
    clearall
    clear results_file
    
    %% same as figure 2, but now for significances
    for e=1:size(cfg.analyse_states,1)
        E=cfg.analyse_states{e,1};
        
        for t = 1: length(targets)
            for c = 1:length(cond)
                k=t+(c-1)*length(targets);
                
                % create figure
                h(c,t) = figure('units','normalized','position',[0 0 1 1]);
                if tar(t).nSites == 0
                    continue;
                end
                
                tfr_time=tar(t).con(c).(E).tfr_time;
                tfr_events=get_ticks_and_labels(tfr_time,cfg.analyse_states(:,1));
                powbp=tar(t).con(c).(E).powbp_sig;
                itpcbp=tar(t).con(c).(E).itpcbp_sig;
                lfp=tar(t).con(c).(E).lfp_sig;
                
                % =========================== Power & ITPC ============================= %
                toplot={tar(t).con(c).(E).pow_sig,tar(t).con(c).(E).itpc_sig};
                for sp=1:2
                    sph(c,sp,t)=subplot(3,2,sp);
                    ax=sph(c,sp,t);
                    hold on
                    image(toplot{sp},'CDataMapping','scaled');
                    set(gca,'YDir','normal');
                    line([0.5 0.5], ylim, 'color', 'k');
                    
                    % horizontal lines to separate frequency bands
                    fbandstart = unique(frequency_bands(:))';
                    fbandstart_idx = zeros(size(fbandstart));
                    for f = fbandstart
                        f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                        line([0.5 size(toplot{sp},2)], [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                        fbandstart_idx(fbandstart == f) = f_idx;
                    end
                    set(ax,'TickDir','out')
                    set(ax, 'ytick', fbandstart_idx);
                    set(ax, 'yticklabel', fbandstart);
                    set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                    set(ax, 'xlim', [0.5 tfr_events.ticksamples(end)] + 0.5);
                    title(plot_names{sp},'fontsize',10,'interpreter','none');
                    nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
                    collim{k,sp}=[min(nonnan(:)) max(nonnan(:))];
                    add_ticks_and_labels(ax,tfr_events,[0.5,numel(freq) + 0.5],8)
                end
                
                %========================== Bandpassed POWER ==================== %
                sp=3;
                sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
                hold on;
                set(ax,'ColorOrder',colsbp);
                smoothed=smoothit(powbp,half_win,gaussian_kernel);
                plot(squeeze(smoothed)')
                xlabel('Time(s)'); ylabel('sig. Power [% of sites]');
                legend(freqb,'fontsize',3);
                title(plot_names{sp},'fontsize',10,'interpreter','none');
                add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
                
                %========================== Bandpassed ITPC ==================== %
                sp=4;
                sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
                hold on;
                set(ax,'ColorOrder',colsbp);
                smoothed=smoothit(itpcbp,half_win,gaussian_kernel);
                plot(squeeze(smoothed)')
                xlabel('Time(s)'); ylabel('sig. ITPC [% of sites]');
                legend(freqb,'fontsize',3);
                title(plot_names{sp},'fontsize',10,'interpreter','none');
                add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
                
                %========================== LFP evoked Potential ==================== %
                sp=5;
                sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
                hold on;
                smoothed_mean=smoothit(lfp,half_win,gaussian_kernel);
                plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
                xlabel('Time(s)'); ylabel('sig. LFP [% of sites]');
                title(plot_names{sp},'fontsize',10,'interpreter','none');
                add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
                
                % format spectra colors
                cbtitle = {'% sig','% sig'};
                for sp=1:2
                    subplot(3,2,sp);
                    cm = colormap('jet');
                    cb = colorbar;
                    set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                    set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
                    colormap(cm);
                end
                
                results_file{c,t} = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-','sig% of ',num2str(tar(t).nSites),' sites ',withunits]);
                mtit(h(c,t),[cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-sig% of ',num2str(tar(t).nSites),' sites ' ,withunits],mtitsettings{:})
            end
        end
        
        for c = 1:length(cond)
            for t=1:length(targets)
                figure(h(c,t));
                for sp=1:2
                    subplot(sph(c,sp,t));
                    set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
                end
                export_fig(h(c,t),[results_file{c,t},'.pdf']);
            end
        end
    end
    clearall
    
    for e=1:size(cfg.analyse_states,1)
        E=cfg.analyse_states{e,1};
        k=0;
        h(e) = figure('units','normalized','position',[0 0 1 1]);
        
        tfr_time=tar(t).con(c).(E).tfr_time;
        tfr_events=get_ticks_and_labels(tfr_time,cfg.analyse_states(e,1));
        for t = 1: length(targets)
            for c = 1:length(cond)
                k=k+1;
                lfp=tar(t).con(c).(E).lfp_sig;
                pval=tar(t).con(c).(E).lfp_p;
                if pval<0.001
                    ptoplot='<0.001';
                else
                    ptoplot=num2str(round(pval*1000)/1000);
                end
                %========================== LFP evoked Potential ==================== %
                sp=5;
                sph(c,sp,t)=subplot(3,2,k);
                ax=sph(c,sp,t);
                hold on;
                smoothed_mean=smoothit(lfp,half_win,gaussian_kernel);
                plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
                xlabel('Time(s)'); ylabel('sig. LFP [% of sites]');
                title([targets{t} '-' cond{c} '-' num2str(tar(t).nSites) 'sites-p=' ptoplot],'fontsize',10,'interpreter','none');
                
                add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
                
            end
        end       %mtit(h(c,t),[cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-sig% of ',num2str(tar(t).nSites),' sites ' ,withunits],mtitsettings{:})
        
        results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',E,'-','sig% of lfp, ',withunits]);
        
        export_fig(h(e),[results_file,'.pdf']);
        %     for c = 1:length(cond)
        %         for t=1:length(targets)
        %             figure(h(c,t));
        %             for sp=1:2
        %                 subplot(sph(c,sp,t));
        %                 set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
        %             end
        %
        %         end
        %     end
    end
    clearall
    
    
    
    for e=1:size(cfg.analyse_states,1)
        E=cfg.analyse_states{e,1};
        
        for t = 1: length(targets)
            for c = 1:length(cond)
                k=t+(c-1)*length(targets);
                
                % create figure
                h(c,t) = figure('units','normalized','position',[0 0 1 1]);
                if tar(t).nSites == 0
                    continue;
                end
                
                tfr_time=tar(t).con(c).(E).tfr_time;
                tfr_events=get_ticks_and_labels(tfr_time,cfg.analyse_states(:,1));
                powbp=tar(t).con(c).(E).powbp_sig;
                itpcbp=tar(t).con(c).(E).itpcbp_sig;
                lfp=tar(t).con(c).(E).lfp_sig;
                
                % =========================== Power & ITPC ============================= %
                toplot={tar(t).con(c).(E).pow_sig,tar(t).con(c).(E).itpc_sig};
                for sp=1:2
                    sph(c,sp,t)=subplot(3,2,sp);
                    ax=sph(c,sp,t);
                    hold on
                    image(toplot{sp},'CDataMapping','scaled');
                    set(gca,'YDir','normal');
                    line([0.5 0.5], ylim, 'color', 'k');
                    
                    % horizontal lines to separate frequency bands
                    fbandstart = unique(frequency_bands(:))';
                    fbandstart_idx = zeros(size(fbandstart));
                    for f = fbandstart
                        f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                        line([0.5 size(toplot{sp},2)], [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                        fbandstart_idx(fbandstart == f) = f_idx;
                    end
                    set(ax,'TickDir','out')
                    set(ax, 'ytick', fbandstart_idx);
                    set(ax, 'yticklabel', fbandstart);
                    set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                    set(ax, 'xlim', [0.5 tfr_events.ticksamples(end)] + 0.5);
                    title(plot_names{sp},'fontsize',10,'interpreter','none');
                    nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
                    collim{k,sp}=[min(nonnan(:)) max(nonnan(:))];
                    add_ticks_and_labels(ax,tfr_events,[0.5,numel(freq) + 0.5],8)
                end
                
                %========================== Bandpassed POWER ==================== %
                sp=3;
                sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
                hold on;
                set(ax,'ColorOrder',colsbp);
                smoothed=smoothit(powbp,half_win,gaussian_kernel);
                plot(squeeze(smoothed)')
                xlabel('Time(s)'); ylabel('sig. Power [% of sites]');
                legend(freqb,'fontsize',3);
                title(plot_names{sp},'fontsize',10,'interpreter','none');
                add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
                
                %========================== Bandpassed ITPC ==================== %
                sp=4;
                sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
                hold on;
                set(ax,'ColorOrder',colsbp);
                smoothed=smoothit(itpcbp,half_win,gaussian_kernel);
                plot(squeeze(smoothed)')
                xlabel('Time(s)'); ylabel('sig. ITPC [% of sites]');
                legend(freqb,'fontsize',3);
                title(plot_names{sp},'fontsize',10,'interpreter','none');
                add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
                
                %========================== LFP evoked Potential ==================== %
                sp=5;
                sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
                hold on;
                smoothed_mean=smoothit(lfp,half_win,gaussian_kernel);
                plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
                xlabel('Time(s)'); ylabel('sig. LFP [% of sites]');
                title(plot_names{sp},'fontsize',10,'interpreter','none');
                add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
                
                % format spectra colors
                cbtitle = {'% sig','% sig'};
                for sp=1:2
                    subplot(3,2,sp);
                    cm = colormap('jet');
                    cb = colorbar;
                    set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                    set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
                    colormap(cm);
                end
                
                results_file{c,t} = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-','sig% of ',num2str(tar(t).nSites),' sites ',withunits]);
                mtit(h(c,t),[cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-sig% of ',num2str(tar(t).nSites),' sites ' ,withunits],mtitsettings{:})
            end
        end
        
        for c = 1:length(cond)
            for t=1:length(targets)
                figure(h(c,t));
                for sp=1:2
                    subplot(sph(c,sp,t));
                    set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
                end
                export_fig(h(c,t),[results_file{c,t},'.pdf']);
            end
        end
    end
    clearall
    
    %% same as figure 2, but now for SIGNED significances
    for e=1:size(cfg.analyse_states,1)
        E=cfg.analyse_states{e,1};
        for t = 1: length(targets)
            for c = 1:length(cond)
                k=t+(c-1)*length(targets);
                nsites=[' of ',num2str(size(tar(t).con(c).concat.max_itpcbp,2)),' sites,'];
                ntriggers=[num2str(tar(t).con(c).concat.nTriggers),' avg nTriggers'];
                ncond=['-' strrep(cond{c},'_',' ') '-'];
                
                % create figure
                h(c,t) = figure('units','normalized','position',[0 0 1 1]);
                if tar(t).nSites == 0
                    continue;
                end
                
                tfr_time=tar(t).con(c).(E).tfr_time;
                tfr_events=get_ticks_and_labels(tfr_time,cfg.analyse_states(:,1));
                powbp=tar(t).con(c).(E).powbp_sig_signed;
                itpcbp=tar(t).con(c).(E).itpcbp_sig_signed;
                lfp=tar(t).con(c).(E).lfp_sig_signed;
                
                % =========================== Power & ITPC ============================= %
                toplot={tar(t).con(c).(E).pow_sig_signed,tar(t).con(c).(E).itpc_sig_signed};
                for sp=1:2
                    sph(c,sp,t)=subplot(3,2,sp);
                    ax=sph(c,sp,t);
                    hold on
                    image(toplot{sp},'CDataMapping','scaled');
                    set(gca,'YDir','normal');
                    
                    % horizontal lines to separate frequency bands
                    fbandstart = unique(frequency_bands(:))';
                    fbandstart_idx = zeros(size(fbandstart));
                    for f = fbandstart
                        f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                        line([0.5 size(toplot{sp},2)], [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                        fbandstart_idx(fbandstart == f) = f_idx;
                    end
                    set(ax,'TickDir','out')
                    set(ax, 'ytick', fbandstart_idx);
                    set(ax, 'yticklabel', fbandstart);
                    set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                    set(ax, 'xlim', [0.5 tfr_events.ticksamples(end)] + 0.5);
                    title(plot_names{sp},'fontsize',10,'interpreter','none');
                    nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
                    collim{k,sp}=[min(nonnan(:)) max(nonnan(:))];
                    add_ticks_and_labels(ax,tfr_events,[0.5,numel(freq) + 0.5],8)
                end
                
                %========================== Bandpassed POWER ==================== %
                sp=3;
                sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
                hold on;
                set(ax,'ColorOrder',colsbp);
                smoothed=smoothit(powbp,half_win,gaussian_kernel);
                plot(squeeze(smoothed)')
                xlabel('Time(s)'); ylabel('sig. Power [% of sites]');
                legend(freqb,'fontsize',3);
                title(plot_names{sp},'fontsize',10,'interpreter','none');
                add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
                
                %========================== Bandpassed ITPC ==================== %
                sp=4;
                sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
                hold on;
                set(ax,'ColorOrder',colsbp);
                smoothed=smoothit(itpcbp,half_win,gaussian_kernel);
                plot(squeeze(smoothed)')
                xlabel('Time(s)'); ylabel('sig. ITPC [% of sites]');
                legend(freqb,'fontsize',3);
                title(plot_names{sp},'fontsize',10,'interpreter','none');
                add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
                
                %========================== LFP evoked Potential ==================== %
                sp=5;
                sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
                hold on;
                smoothed_mean=smoothit(lfp,half_win,gaussian_kernel);
                plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
                xlabel('Time(s)'); ylabel('sig. LFP [% of sites]');
                title(plot_names{sp},'fontsize',10,'interpreter','none');
                add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
                
                %% format spectra colors
                cbtitle = {'% sig','% sig'};
                for sp=1:2
                    subplot(3,2,sp);
                    cm = colormap('jet');
                    cb = colorbar;
                    set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                    set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
                    colormap(cm);
                end
                
                results_file{c,t} = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{t},ncond,E,'-','signed sig%',nsites,withunits]);
                mtit(h(c,t),[ cfg.monkey,'-',targets{t},ncond,E,'- (signed) sig%',nsites,withunits],mtitsettings{:})
            end
        end
        
        for c = 1:length(cond)
            for t=1:length(targets)
                figure(h(c,t));
                for sp=1:2
                    subplot(sph(c,sp,t));
                    set(gca,'CLim',[min([collim{:,sp}]) max([collim{:,sp}])]);
                end
                export_fig(h(c,t),[results_file{c,t},'.pdf']);
            end
        end
    end
    clearall
    
    %% plotting max values versus time of max for ITPC/POW:
    for e=1:size(cfg.analyse_states,1)
        E=cfg.analyse_states{e,1};
        for t = 1:length(targets)
            if tar(t).nSites == 0
                continue;
            end
            h = figure('Name',['Max ITPCbp/POWbp in ' E ' of Target=' tar(t).target],'NumberTitle','off');
            for c = 1:length(cond)
                for fb = 1: length(freqb)
                    % itpcbp
                    subplot(2,2,2*c-1)
                    scatter(tar(t).con(c).(E).max_itpcbp_time(fb,:),tar(t).con(c).(E).max_itpcbp(fb,:),15,colsbp(fb,:))
                    hold on
                    title([' max itpc-bp in ', strrep(cond{c},'_',' '),' for ',num2str(tar(t).nSites),...
                        ' sites,',num2str(tar(t).con(c).(E).nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
                    xlabel('time (s)','Interpreter','latex'), ylabel('Max ITPC value','Interpreter','latex');
                    
                    % powbp
                    subplot(2,2,2*c)
                    scatter(tar(t).con(c).(E).max_powbp_time(fb,:),tar(t).con(c).(E).max_powbp(fb,:),15,colsbp(fb,:))
                    hold on
                    title([' max power-bp in ', strrep(cond{c},'_',' '),' for ',num2str(tar(t).nSites),...
                        ' sites,',num2str(tar(t).con(c).concat.nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
                    xlabel('time (s)','Interpreter','latex'), ylabel('Max POWER value','Interpreter','latex');
                    
                end
                
                subplot(2,2,2*c-1)
                plot([0,0],get(gca,'ylim'),'k--')
                
                subplot(2,2,2*c)
                plot([0,0],get(gca,'ylim'),'k--')
                legend(freqb,'FontSize',5)
                legend('boxoff')
                
            end
            mtit([strrep(tar(t).target,'_','-') '-' E],mtitsettings{:});
            results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{t},'-',E,'-Max_ITPCbp_POWbp_of ',num2str(tar(t).nSites),' sites ', withunits]);
            export_fig(h,[results_file,'.pdf']);
        end
        close all,
        clc
    end
    clearall
    
    %% plotting the avgogram of the timing of the Max ITPC/POW:
    for e=1:size(cfg.analyse_states,1)
        E=cfg.analyse_states{e,1};
        bins=cfg.analyse_states{e,4}:cfg.lfp.timestep*10:cfg.analyse_states{e,5};
        for t = 1:length(targets)
            if tar(t).nSites == 0
                continue;
            end
            h = figure('Name',['Max ITPCbp/POWbp in ' E ' of Target=' tar(t).target],'NumberTitle','off');
            for c = 1:length(cond)
                ncond=strrep(cond{c},'_',' ');
                ntriggers=[num2str(tar(t).con(c).(E).nTriggers),' avg nTriggers'];
                ns=size(tar(t).con(c).(E).max_itpcbp_time,2);
                nsites=[' of ',num2str(ns),' sites,'];
                for fb = 1: length(freqb)
                    % itpcbp
                    sp1(c)=subplot(2,2,2*c-1);
                    
                    [nelements,centers]= hist(squeeze(tar(t).con(c).(E).max_itpcbp_time(fb,:)),bins);
                    nelements=nelements/ns;
                    plot(centers,nelements,'-','Color',colsbp(fb,:));
                    hold on
                    title([' max itpc-bp in ',E,' ',ncond,nsites,ntriggers],'FontSize',6,'Interpreter','latex');
                    xlabel('Max ITPC times','Interpreter','latex');
                    set(gca,'xlim',[bins(1) bins(end)])
                    
                    % powbp
                    sp2(c)=subplot(2,2,2*c);
                    [nelements,centers]= hist(squeeze(tar(t).con(c).(E).max_powbp_time(fb,:)),bins);
                    nelements=nelements/ns;
                    plot(centers,nelements,'-','Color',colsbp(fb,:));
                    hold on
                    title([' max power-bp in ',E,' ',ncond,nsites,ntriggers],'FontSize',6,'Interpreter','latex');
                    xlabel('Max POWER times','Interpreter','latex');
                    set(gca,'xlim',[bins(1) bins(end)])
                end
            end
            legend(freqb,'FontSize',5)
            legend('boxoff')
            ylims1=get(sp1,'ylim');
            set(sp1,'ylim',[0 max([ylims1{:}])]);
            ylims2=get(sp2,'ylim');
            set(sp2,'ylim',[0 max([ylims2{:}])]);
            
            
            for c = 1:length(cond)
                plot(sp1(c),[0,0],get(sp1(c),'ylim'),'k--')
                plot(sp2(c),[0,0],get(sp2(c),'ylim'),'k--')
            end
            
            mtit([strrep(tar(t).target,'_','-') '-' E],mtitsettings{:});
            results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{t},'-',E,'-Time_of_Max_ITPCbp_POWbp ',num2str(tar(t).nSites),' sites ', withunits]);
            export_fig(h,[results_file,'.pdf']);
        end
        close all,
        clc
    end
    clearall
    
    % %% plotting the avgogram of the timing of the Max ITPC/POW:
    % bins=cfg.analyse_states{1,4}:cfg.lfp.timestep*10:cfg.analyse_states{1,5};
    % for t = 1: length(targets)
    %     if tar(t).nSites == 0
    %         continue;
    %     end
    %     h = figure;
    %     for c = 1:length(cond)
    %         nsites=[' for ',num2str(size(tar(t).con(c).concat.max_itpcbp,2)),' sites,'];
    %         ntriggers=[num2str(tar(t).con(c).concat.nTriggers),' avg nTriggers'];
    %         ncond=strrep(cond{c},'_',' ');
    %         for fb = 1: length(freqb)
    %
    %             % itpcbp
    %             sp1=subplot(2,2,2*c-1);
    %             [nelements,centers]= hist(squeeze(tar(t).con(c).concat.max_itpcbp_time(fb,:)),bins);
    %             plot(centers,nelements,'-','Color',colsbp(fb,:));
    %             hold on
    %             title([' max itpc-bp in ',ncond,nsites,ntriggers],'FontSize',6,'Interpreter','latex');
    %             xlabel('Max ITPC times','Interpreter','latex');
    %             set(gca,'xlim',[-0.25,0.25])
    %
    %             % powbp
    %             sp2=subplot(2,2,2*c);
    %             [nelements,centers]= hist(squeeze(tar(t).con(c).concat.max_powbp_time(fb,:)),bins);
    %             plot(centers,nelements,'-','Color',colsbp(fb,:));
    %             hold on
    %             title([' max power-bp in ',ncond,nsites,ntriggers],'FontSize',6,'Interpreter','latex');
    %             xlabel('Max POWER times','Interpreter','latex');
    %             set(gca,'xlim',[-0.25,0.25])
    %         end
    %         legend(freqb,'FontSize',5)
    %         legend('boxoff')
    %         plot(sp1,[0,0],get(sp1,'ylim'),'k--')
    %         plot(sp2,[0,0],get(sp2,'ylim'),'k--')
    %     end
    %     mtit(['Target=',strrep(tar(t).target,'_','-')],mtitsettings{:});
    %     results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{t},'-','Time_of_Max_ITPCbp_POWbp',nsites,withunits]);
    %     export_fig(h,[results_file,'.pdf']);
    % end
    % clearall
    %
end

end

function clearall
close all,
clear results_file sph h collim
clc
end

function out=smoothit(in,half_win,gaussian_kernel)
jnk=[];
concat_input = cat(2,(in(:,half_win:-1:1)),in);
concat_input = cat(2,concat_input, (in(:,end:-1:end-half_win+1)));
for k=1:size(concat_input,1)
    jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');
end
out = jnk(:,half_win+1:end-half_win);
end

function tfr_events=get_ticks_and_labels(time,eventnames)
tfr_events.onset        = find(time == 0);
tfr_events.name         = eventnames;
tfr_events.ticksamples  = [1,numel(time)];
%sort([find(diff([0 ~isnan(tfr_time)])==1) find(diff([~isnan(tfr_time)])==-1)]);
tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
tfr_events.ticks        = round(time(tfr_events.ticksamples)*100)/100;
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
