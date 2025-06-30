function tar = ecg_bna_compute_grand_avg(cfg,withunits)

reprocess=0;
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
                tfr_time  = event.tfr_time;
                lfp       = squeeze(event.real.lfp.mean)';
%                 lfp       = squeeze(event.real.lfp.mean)' - squeeze(event.shuffled.lfp.mean)';
%                 lfp       = squeeze(event.normalized.lfp.mean)';
                itpc      = squeeze(event.real.itpc.mean)-squeeze(event.shuffled.itpc.mean);
                power     = squeeze(event.normalized.pow.mean);
                itpcbp    = squeeze(event.real.itpcbp.mean)-squeeze(event.shuffled.itpcbp.mean);
                powerbp   = squeeze(event.normalized.powbp.mean);
                
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

startfreq=3.5;
colsbp=jet(length(cfg.lfp.freqb));
fbx=find(cfg.lfp.frequency_bands(:,1)>=startfreq,1,'first');
frequency_bands=cfg.lfp.frequency_bands(fbx:end,:);
freqb = cfg.lfp.freqb(fbx:end);
%freqName = cfg.lfp.freqName(fbx:end);
colsbp=colsbp(fbx:end,:);

fx=find(cfg.lfp.foi>=startfreq,1,'first');
freq = cfg.lfp.foi(fx:end);


%% Computing the Target-wise averaging of total available sites
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
            E.nTriggers = mean(cat(3,events(:,e).nTriggers),3);
            E.pow       = mean(cat(3,events(:,e).pow),3);
            E.itpc      = mean(cat(3,events(:,e).itpc),3);
            E.lfp       = mean(cat(3,events(:,e).lfp),3);
            E.lfp_25    = prctile(cat(3,events(:,e).lfp),25,3);
            E.lfp_75    = prctile(cat(3,events(:,e).lfp),75,3);
            E.itpcbp    = mean(cat(3,events(:,e).itpcbp),3);
            E.powbp       = mean(cat(3,events(:,e).powbp),3);
            E.tfr_time    = events(1,e).tfr_time;
            E.lfp_time    = events(1,e).time;
            E.lfp_std     = std(cat(3,events(:,e).lfp),0,3);
                        
            E.pow_sig     = mean(cat(3,events(:,e).pow_sig),3);
            E.itpc_sig    = mean(cat(3,events(:,e).itpc_sig),3);
            E.powbp_sig   = mean(cat(3,events(:,e).powbp_sig),3);
            E.itpcbp_sig  = mean(cat(3,events(:,e).itpcbp_sig),3);
            E.lfp_sig     = mean(cat(3,events(:,e).lfp_sig),3);
            
            E.pow_sig_signed     = mean(cat(3,events(:,e).pow_sig).*sign(cat(3,events(:,e).pow)),3);
            E.itpc_sig_signed    = mean(cat(3,events(:,e).itpc_sig).*sign(cat(3,events(:,e).itpc)),3);
            E.powbp_sig_signed   = mean(cat(3,events(:,e).powbp_sig).*sign(cat(3,events(:,e).powbp)),3);
            E.itpcbp_sig_signed  = mean(cat(3,events(:,e).itpcbp_sig).*sign(cat(3,events(:,e).itpcbp)),3);
            E.lfp_sig_signed     = mean(cat(3,events(:,e).lfp_sig).*sign(cat(3,events(:,e).lfp)),3);
                        
            
            E.pow       = E.pow(fx:end,:);
            E.itpc       = E.itpc(fx:end,:);
            E.itpcbp =E.itpcbp(fbx:end,:);
            E.powbp =E.powbp(fbx:end,:);
            
            E.pow_sig    = E.pow_sig(fx:end,:)*100;
            E.itpc_sig   = E.itpc_sig(fx:end,:)*100;  
            E.powbp_sig  = E.powbp_sig(fbx:end,:)*100;
            E.itpcbp_sig = E.itpcbp_sig(fbx:end,:)*100;            
            E.lfp_sig    = E.lfp_sig*100;
            
            E.pow_sig_signed    = E.pow_sig_signed(fx:end,:)*100;
            E.itpc_sig_signed   = E.itpc_sig_signed(fx:end,:)*100;  
            E.powbp_sig_signed  = E.powbp_sig_signed(fbx:end,:)*100;
            E.itpcbp_sig_signed = E.itpcbp_sig_signed(fbx:end,:)*100;            
            E.lfp_sig_signed    = E.lfp_sig_signed*100;            
            
            
            E.max_itpcbp=[events(:,e).max_itpcbp];
            E.max_powbp=[events(:,e).max_powbp];
            E.max_itpcbp_time=[events(:,e).max_itpcbp_time];
            E.max_powbp_time=[events(:,e).max_powbp_time];
            
            
            %(fbx:end,:)
            
            
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


%% plotting the results:
cond = {cfg.condition.name};
plot_names={'POW','ITPC','Power_BP','ITPC_BP','LFP_Evoked'};

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
        tfr_time=tar(t).con(c).concat.tfr_time;
        tfr_events.onset        = find(tfr_time == 0);
        tfr_events.name         = cfg.analyse_states(:,1);
        tfr_events.ticksamples  = sort([find(diff([0 ~isnan(tfr_time)])==1) find(diff(~isnan(tfr_time))==-1)]);
        tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
        tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
        tfr_events.ticks        = round(tfr_time(tfr_events.ticksamples)*100)/100;
        
        powbp=tar(t).con(c).concat.powbp;
        itpcbp=tar(t).con(c).concat.itpcbp;
        lfp=tar(t).con(c).concat.lfp;
        lfp_std=tar(t).con(c).concat.lfp_std;
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
        % Smoothing of the POWbp here:
        jnk = [];
        concat_input = cat(2,(powbp(:,half_win:-1:1)),(powbp(:,:)));
        concat_input = cat(2,concat_input, (powbp(:,end:-1:end-half_win+1)));
        for k=1:size(powbp,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed = jnk(:,half_win+1:end-half_win);
        
        sp=3;
        sph(c,sp,t)=subplot(3,2,sp);
        ax=sph(c,sp,t);
        hold on;
        set(ax,'ColorOrder',colsbp);
        plot(squeeze(smoothed)')
        xlabel('Time(s)'); ylabel('Power (W)');
        legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
        title(plot_names{sp},'fontsize',10,'interpreter','none');
        add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
        
        %========================== Bandpassed ITPC ==================== %
        % Smoothing of the ITPCbp here:
        jnk = [];
        concat_input = cat(2,(itpcbp(:,half_win:-1:1)),(itpcbp(:,:)));
        concat_input = cat(2,concat_input, (itpcbp(:,end:-1:end-half_win+1)));
        for k=1:size(itpcbp,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed = jnk(:,half_win+1:end-half_win);
        
        sp=4;
        sph(c,sp,t)=subplot(3,2,sp);
        ax=sph(c,sp,t);
        hold on;
        set(ax,'ColorOrder',colsbp);
        plot(squeeze(smoothed)')
        xlabel('Time(s)'); ylabel('ITPC');
        legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
        title(plot_names{sp},'fontsize',10,'interpreter','none');
        add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
        
        %========================== LFP evoked Potential ==================== %
        % Smoothing of the  LFP evoked Potential here:
        jnk = [];
        concat_input = cat(2,(lfp(:,half_win:-1:1)),(lfp(:,:)));
        concat_input = cat(2,concat_input, (lfp(:,end:-1:end-half_win+1)));
        for k=1:size(lfp,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed_mean = jnk(:,half_win+1:end-half_win);
        concat_input = cat(2,(lfp_std(:,half_win:-1:1)),(lfp_std(:,:)));
        concat_input = cat(2,concat_input, (lfp_std(:,end:-1:end-half_win+1)));
        for k=1:size(lfp,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed_std = jnk(:,half_win+1:end-half_win);
        concat_input = cat(2,(percentile25(:,half_win:-1:1)),(percentile25(:,:)));
        concat_input = cat(2,concat_input, (percentile25(:,end:-1:end-half_win+1)));
        for k=1:size(lfp,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed_25 = jnk(:,half_win+1:end-half_win);
        concat_input = cat(2,(percentile75(:,half_win:-1:1)),(percentile75(:,:)));
        concat_input = cat(2,concat_input, (percentile75(:,end:-1:end-half_win+1)));
        for k=1:size(lfp,1)
            jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
        end
        clear concat_input
        smoothed_75 = jnk(:,half_win+1:end-half_win);
        
        sp=5;
        sph(c,sp,t)=subplot(3,2,sp);
        ax=sph(c,sp,t);
        hold on;
        %set(ax,'ColorOrder',jet(size(lfp,1)));
        lineProps={'color',[0 0 1]};
        %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,[smoothed_75-smoothed_mean;smoothed_mean-smoothed_25 ],lineProps,1);        
        shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
        %plot(repmat(time,size(tar(t).con(c).lfp_avg,1),1)', squeeze(smoothed_mean)')
        
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
        mtit(h(c,t),[ cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-avg of ',num2str(tar(t).nSites),' sites ' ,withunits],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
        
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
close all,
clear results_file sph h collim


% same as figure 1, but separately for each event
for e=1:size(cfg.analyse_states,1)
    E=cfg.analyse_states{e,1};
    
    for t = 1: length(targets)
        for c = 1:length(cond)
            k=t+(c-1)*length(targets);
            bins=tar(t).con(c).(E).tfr_time;
            % create figure
            h(c,t) = figure('units','normalized','position',[0 0 1 1]);
            if tar(t).nSites == 0
                continue;
            end
            
            tfr_time=tar(t).con(c).(E).tfr_time;
            tfr_events.onset        = find(tfr_time == 0);
            tfr_events.name         = cfg.analyse_states(e,1);
            tfr_events.ticksamples  = [1,numel(tfr_time)]; 
            tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
            tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
            tfr_events.ticks        = round(tfr_time(tfr_events.ticksamples)*100)/100;
            
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
            for sp=1:2
                % =========================== Power & ITPC ============================= %
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
            % Smoothing of the POWbp here:
            jnk = [];
            concat_input = cat(2,(powbp(:,half_win:-1:1)),(powbp(:,:)));
            concat_input = cat(2,concat_input, (powbp(:,end:-1:end-half_win+1)));
            for k=1:size(powbp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed = jnk(:,half_win+1:end-half_win);
            
            sp=3;
            sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
            hold on;
            set(ax,'ColorOrder',colsbp);
            plot(squeeze(smoothed)')
            xlabel('Time(s)'); ylabel('Power (W)');
            legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
            
            %========================== Bandpassed ITPC ==================== %
            % Smoothing of the ITPCbp here:
            jnk = [];
            concat_input = cat(2,(itpcbp(:,half_win:-1:1)),(itpcbp(:,:)));
            concat_input = cat(2,concat_input, (itpcbp(:,end:-1:end-half_win+1)));
            for k=1:size(itpcbp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed = jnk(:,half_win+1:end-half_win);
            
            sp=4;
            sph(c,sp,t)=subplot(3,2,sp);
            ax=sph(c,sp,t);
            hold on;
            set(ax,'ColorOrder',colsbp);
            plot(squeeze(smoothed)')
            xlabel('Time(s)'); ylabel('ITPC');
            legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
            
            %========================== LFP evoked Potential ==================== %
            % Smoothing of the  LFP evoked Potential here:
            jnk = [];
            concat_input = cat(2,(lfp(:,half_win:-1:1)),(lfp(:,:)));
            concat_input = cat(2,concat_input, (lfp(:,end:-1:end-half_win+1)));
            for k=1:size(lfp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_mean = jnk(:,half_win+1:end-half_win);
            concat_input = cat(2,(lfp_std(:,half_win:-1:1)),(lfp_std(:,:)));
            concat_input = cat(2,concat_input, (lfp_std(:,end:-1:end-half_win+1)));
            for k=1:size(lfp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_std = jnk(:,half_win+1:end-half_win);
            concat_input = cat(2,(percentile25(:,half_win:-1:1)),(percentile25(:,:)));
            concat_input = cat(2,concat_input, (percentile25(:,end:-1:end-half_win+1)));
            for k=1:size(lfp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_25 = jnk(:,half_win+1:end-half_win);
            concat_input = cat(2,(percentile75(:,half_win:-1:1)),(percentile75(:,:)));
            concat_input = cat(2,concat_input, (percentile75(:,end:-1:end-half_win+1)));
            for k=1:size(lfp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_75 = jnk(:,half_win+1:end-half_win);
            
            sp=5;
            sph(c,sp,t)=subplot(3,2,sp);
            ax=sph(c,sp,t);
            hold on;
            %set(ax,'ColorOrder',jet(size(lfp,1)));
            lineProps={'color',[0 0 1]};
            %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,[smoothed_75-smoothed_mean;smoothed_mean-smoothed_25 ],lineProps,1);
            
            shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
            %plot(repmat(time,size(tar(t).con(c).lfp_avg,1),1)', squeeze(smoothed_mean)')
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
            
            results_file{c,t} = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-','avg of ',num2str(tar(t).nSites),' sites ',withunits]);
            mtit(h(c,t),[ cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-avg of ',num2str(tar(t).nSites),' sites ' ,withunits],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
            
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
close all,

clear results_file sph h collim
clc

%<<<<<<< HEAD
% same as figure 1, but now for significances
for e=1:size(cfg.analyse_states,1)
    E=cfg.analyse_states{e,1};
    
    for t = 1: length(targets)
        for c = 1:length(cond)
            k=t+(c-1)*length(targets);
            bins=tar(t).con(c).(E).tfr_time;
            % create figure
            h(c,t) = figure('units','normalized','position',[0 0 1 1]);
            if tar(t).nSites == 0
                continue;
            end
% =======
% %%
% freqb = {'delta 2-4 Hz','theta 4-8 Hz','alpha 8-14 Hz','beta 14-30 Hz','lowgamma 30-50 Hz','highgamma 70-150 Hz'};
% freqName = {'delta','theta','alpha','beta','lowGamma','highGamma'};
% % freqb = {'theta 4-8 Hz','alpha 8-14 Hz','beta 14-30 Hz','lowgamma 30-50 Hz','highgamma 70-150 Hz'};
% % freqName = {'theta','alpha','beta','lowGamma','highGamma'};
% color = jet(length(frequency_bands));
% 
% for tr = 1:length(targets)
%     if grand_avg(tr).nSites == 0
%         continue;
%     end
%     h = figure('Name',['Max ITPCbp/POWbp of Target=',grand_avg(tr).target],'NumberTitle','off');
%     for c = 1:length(cond)
%         for fb = 1: length(freqb)
%             % itpcbp
%             subplot(2,2,2*c-1)
%             scatter(grand_avg(tr).avg(c).max_itpcbp_time(fb,:),grand_avg(tr).avg(c).max_itpcbp(fb,:),15,color(fb,:))
%             hold on
%             title([' max itpc-bp in ', strrep(cond{c},'_',' '),' for ',num2str(grand_avg(tr).nSites),...
%                 ' sites,',num2str(grand_avg(tr).avg(c).nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('time (s)','Interpreter','latex'), ylabel('Max ITPC value','Interpreter','latex');
% %             set(gca,'xlim',[-0.5,0.5],'ylim',[0,1])
%             set(gca,'xlim',[-0.25,0.25],'ylim',[0,1])
% >>>>>>> 86fc46a25bdc94ee94288cf71ea552db2eba85fa
            
            % tfr_time not defined well!
            tfr_time=tar(t).con(c).(E).tfr_time;
            %lfp_time=tar(t).con(c).concat.lfp_time;
            tfr_events.onset        = find(tfr_time == 0);
            tfr_events.name         = cfg.analyse_states(e,1);
            tfr_events.ticksamples  = [1,numel(tfr_time)]; %sort([find(diff([0 ~isnan(tfr_time)])==1) find(diff([~isnan(tfr_time)])==-1)]);
            tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
            tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
            tfr_events.ticks        = round(tfr_time(tfr_events.ticksamples)*100)/100;
            
            powbp=tar(t).con(c).(E).powbp_sig;
            itpcbp=tar(t).con(c).(E).itpcbp_sig;
            lfp=tar(t).con(c).(E).lfp_sig;
            lfp_std=tar(t).con(c).(E).lfp_std;
            %lfp_se=sterr(lfp,3,0);
            percentile25=tar(t).con(c).(E).lfp_25;
            percentile75=tar(t).con(c).(E).lfp_75;
            percentile25(isnan(percentile25))=0;
            percentile75(isnan(percentile75))=0;
            
            toplot={tar(t).con(c).(E).pow_sig,tar(t).con(c).(E).itpc_sig};
            for sp=1:2
                % =========================== Power & ITPC ============================= %
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
            % Smoothing of the POWbp here:
            jnk = [];
            concat_input = cat(2,(powbp(:,half_win:-1:1)),(powbp(:,:)));
            concat_input = cat(2,concat_input, (powbp(:,end:-1:end-half_win+1)));
            for k=1:size(powbp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed = jnk(:,half_win+1:end-half_win);
            
            sp=3;
            sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
            hold on;
            set(ax,'ColorOrder',colsbp);
            plot(squeeze(smoothed)')
            xlabel('Time(s)'); ylabel('sig. Power [% of sites]');
            legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
            
            %========================== Bandpassed ITPC ==================== %
            % Smoothing of the ITPCbp here:
            jnk = [];
            concat_input = cat(2,(itpcbp(:,half_win:-1:1)),(itpcbp(:,:)));
            concat_input = cat(2,concat_input, (itpcbp(:,end:-1:end-half_win+1)));
            for k=1:size(itpcbp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed = jnk(:,half_win+1:end-half_win);
            
            sp=4;
            sph(c,sp,t)=subplot(3,2,sp);
            ax=sph(c,sp,t);
            hold on;
            set(ax,'ColorOrder',colsbp);
            plot(squeeze(smoothed)')
            xlabel('Time(s)'); ylabel('sig. ITPC [% of sites]');
            legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
            
            %========================== LFP evoked Potential ==================== %
            % Smoothing of the  LFP evoked Potential here:
            jnk = [];
            concat_input = cat(2,(lfp(:,half_win:-1:1)),(lfp(:,:)));
            concat_input = cat(2,concat_input, (lfp(:,end:-1:end-half_win+1)));
            for k=1:size(lfp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_mean = jnk(:,half_win+1:end-half_win);
            concat_input = cat(2,(lfp_std(:,half_win:-1:1)),(lfp_std(:,:)));
            concat_input = cat(2,concat_input, (lfp_std(:,end:-1:end-half_win+1)));
            for k=1:size(lfp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_std = jnk(:,half_win+1:end-half_win);
            concat_input = cat(2,(percentile25(:,half_win:-1:1)),(percentile25(:,:)));
            concat_input = cat(2,concat_input, (percentile25(:,end:-1:end-half_win+1)));
            for k=1:size(lfp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_25 = jnk(:,half_win+1:end-half_win);
            concat_input = cat(2,(percentile75(:,half_win:-1:1)),(percentile75(:,:)));
            concat_input = cat(2,concat_input, (percentile75(:,end:-1:end-half_win+1)));
            for k=1:size(lfp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_75 = jnk(:,half_win+1:end-half_win);
            
            sp=5;
            sph(c,sp,t)=subplot(3,2,sp);
            ax=sph(c,sp,t);
            hold on;
            plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
            
            %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
            %plot(repmat(time,size(tar(t).con(c).lfp_avg,1),1)', squeeze(smoothed_mean)')
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
            
            results_file{c,t} = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-','sig% of ',num2str(tar(t).nSites),' sites ',withunits]);
            mtit(h(c,t),[cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-sig% of ',num2str(tar(t).nSites),' sites ' ,withunits],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
            
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
close all,

clear results_file sph h collim
clc


% same as figure 1, but now for SIGNED significances
for e=1:size(cfg.analyse_states,1)
    E=cfg.analyse_states{e,1};
    
    for t = 1: length(targets)
        for c = 1:length(cond)
            k=t+(c-1)*length(targets);
            bins=tar(t).con(c).(E).tfr_time;
            % create figure
            h(c,t) = figure('units','normalized','position',[0 0 1 1]);
            if tar(t).nSites == 0
                continue;
            end
            
            % tfr_time not defined well!
            tfr_time=tar(t).con(c).(E).tfr_time;
            %lfp_time=tar(t).con(c).concat.lfp_time;
            tfr_events.onset        = find(tfr_time == 0);
            tfr_events.name         = cfg.analyse_states(e,1);
            tfr_events.ticksamples  = [1,numel(tfr_time)]; %sort([find(diff([0 ~isnan(tfr_time)])==1) find(diff([~isnan(tfr_time)])==-1)]);
            tfr_events.startsamples = tfr_events.ticksamples(1:3:end);
            tfr_events.endsamples   = tfr_events.ticksamples(3:3:end);
            tfr_events.ticks        = round(tfr_time(tfr_events.ticksamples)*100)/100;
            
            powbp=tar(t).con(c).(E).powbp_sig_signed;
            itpcbp=tar(t).con(c).(E).itpcbp_sig_signed;
            lfp=tar(t).con(c).(E).lfp_sig_signed;
            lfp_std=tar(t).con(c).(E).lfp_std;
            %lfp_se=sterr(lfp,3,0);
            percentile25=tar(t).con(c).(E).lfp_25;
            percentile75=tar(t).con(c).(E).lfp_75;
            percentile25(isnan(percentile25))=0;
            percentile75(isnan(percentile75))=0;
            
            toplot={tar(t).con(c).(E).pow_sig_signed,tar(t).con(c).(E).itpc_sig_signed};
            for sp=1:2
                % =========================== Power & ITPC ============================= %
                sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
                hold on
                image(toplot{sp},'CDataMapping','scaled');
                set(gca,'YDir','normal');
                %line([0 0], ylim, 'color', 'k');
                
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
            % Smoothing of the POWbp here:
            jnk = [];
            concat_input = cat(2,(powbp(:,half_win:-1:1)),(powbp(:,:)));
            concat_input = cat(2,concat_input, (powbp(:,end:-1:end-half_win+1)));
            for k=1:size(powbp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed = jnk(:,half_win+1:end-half_win);
            
            sp=3;
            sph(c,sp,t)=subplot(3,2,sp);
                ax=sph(c,sp,t);
            hold on;
            set(ax,'ColorOrder',colsbp);
            plot(squeeze(smoothed)')
            xlabel('Time(s)'); ylabel('sig. Power [% of sites]');
            legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
            
            %========================== Bandpassed ITPC ==================== %
            % Smoothing of the ITPCbp here:
            jnk = [];
            concat_input = cat(2,(itpcbp(:,half_win:-1:1)),(itpcbp(:,:)));
            concat_input = cat(2,concat_input, (itpcbp(:,end:-1:end-half_win+1)));
            for k=1:size(itpcbp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed = jnk(:,half_win+1:end-half_win);
            
            sp=4;
            sph(c,sp,t)=subplot(3,2,sp);
            ax=sph(c,sp,t);
            hold on;
            set(ax,'ColorOrder',colsbp);
            plot(squeeze(smoothed)')
            xlabel('Time(s)'); ylabel('sig. ITPC [% of sites]');
            legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
            title(plot_names{sp},'fontsize',10,'interpreter','none');
            add_ticks_and_labels(ax,tfr_events,get(ax,'ylim'),8)
            
            %========================== LFP evoked Potential ==================== %
            % Smoothing of the  LFP evoked Potential here:
            jnk = [];
            concat_input = cat(2,(lfp(:,half_win:-1:1)),(lfp(:,:)));
            concat_input = cat(2,concat_input, (lfp(:,end:-1:end-half_win+1)));
            for k=1:size(lfp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_mean = jnk(:,half_win+1:end-half_win);
            concat_input = cat(2,(lfp_std(:,half_win:-1:1)),(lfp_std(:,:)));
            concat_input = cat(2,concat_input, (lfp_std(:,end:-1:end-half_win+1)));
            for k=1:size(lfp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_std = jnk(:,half_win+1:end-half_win);
            concat_input = cat(2,(percentile25(:,half_win:-1:1)),(percentile25(:,:)));
            concat_input = cat(2,concat_input, (percentile25(:,end:-1:end-half_win+1)));
            for k=1:size(lfp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_25 = jnk(:,half_win+1:end-half_win);
            concat_input = cat(2,(percentile75(:,half_win:-1:1)),(percentile75(:,:)));
            concat_input = cat(2,concat_input, (percentile75(:,end:-1:end-half_win+1)));
            for k=1:size(lfp,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
            clear concat_input
            smoothed_75 = jnk(:,half_win+1:end-half_win);
            
            sp=5;
            sph(c,sp,t)=subplot(3,2,sp);
            ax=sph(c,sp,t);
            hold on;
            plot(1:numel(smoothed_mean),smoothed_mean,'color',[0 0 1]);
            
            %shadedErrorBar(1:numel(smoothed_mean),smoothed_mean,smoothed_std,lineProps,1);
            %plot(repmat(time,size(tar(t).con(c).lfp_avg,1),1)', squeeze(smoothed_mean)')
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
            
            results_file{c,t} = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'-','signed sig% of ',num2str(tar(t).nSites),' sites ',withunits]);
            mtit(h(c,t),[ cfg.monkey,'-',targets{t},'-',tar(t).con(c).cond_name,'-',E,'- (signed) sig% of ',num2str(tar(t).nSites),' sites ' ,withunits],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
            
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
close all,

clear results_file sph h collim
clc

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
        mtit([strrep(tar(t).target,'_','-') '-' E],'xoff', 0, 'yoff', 0.05, 'fontsize', 12,'Interpreter', 'none');
        results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{t},'-',E,'-Max_ITPCbp_POWbp_of ',num2str(tar(t).nSites),' sites ', withunits]);
        export_fig(h,[results_file,'.pdf']);
    end
    close all,
    clc    
end

% for t = 1:length(targets)
%     if tar(t).nSites == 0
%         continue;
%     end
%     h = figure('Name',['Max ITPCbp/POWbp of Target=',tar(t).target],'NumberTitle','off');
%     for c = 1:length(cond)
%         for fb = 1: length(freqb)
%             % itpcbp
%             subplot(2,2,2*c-1)
%             scatter(tar(t).con(c).concat.max_itpcbp_time(fb,:),tar(t).con(c).concat.max_itpcbp(fb,:),15,colsbp(fb,:))
%             hold on
%             title([' max itpc-bp in ', strrep(cond{c},'_',' '),' for ',num2str(tar(t).nSites),...
%                 ' sites,',num2str(tar(t).con(c).concat.nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('time (s)','Interpreter','latex'), ylabel('Max ITPC value','Interpreter','latex');
% %             set(gca,'xlim',[-0.5,0.5],'ylim',[0,1])
%             set(gca,'xlim',[-0.25,0.25],'ylim',[0,1])
%             
%             % powbp
%             subplot(2,2,2*c)
%             scatter(tar(t).con(c).concat.max_powbp_time(fb,:),tar(t).con(c).concat.max_powbp(fb,:),15,colsbp(fb,:))
%             hold on
%             title([' max power-bp in ', strrep(cond{c},'_',' '),' for ',num2str(tar(t).nSites),...
%                 ' sites,',num2str(tar(t).con(c).concat.nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('time (s)','Interpreter','latex'), ylabel('Max POWER value','Interpreter','latex');
% %             set(gca,'xlim',[-0.5,0.5])
%             set(gca,'xlim',[-0.25,0.25])
%             
%         end
%         
%         subplot(2,2,2*c-1)
%         plot([0,0],get(gca,'ylim'),'k--')
%         %
%         subplot(2,2,2*c)
%         plot([0,0],get(gca,'ylim'),'k--')
%         legend(freqb,'FontSize',5)
%         legend('boxoff')
%         
%     end
%     mtit(['Target=',strrep(tar(t).target,'_','-')],'xoff', 0, 'yoff', 0.05, 'Color','red', 'fontsize', 12);
%     results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{t},' - ','Max_ITPCbp_POWbp_of ',num2str(tar(t).nSites),' sites ', withunits]);
%     %     mtit([ ecg_bna_cfg.monkey,'-',targets{t},'-avg of ',num2str(tar(t).nSites),' sites in all sessions - ',tar(t).con(c).cond_name],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
%     export_fig(h,[results_file,'.pdf']);
% %     saveas(h,[results_file]);
% end
% close all,
% clc

%%
% ploting the avgogram of the timing of the Max ITPC/POW:
bins=cfg.analyse_states{1,4}:cfg.lfp.timestep*10:cfg.analyse_states{1,5};
for t = 1: length(targets)
    if tar(t).nSites == 0
        continue;
    end
    h = figure;
    for c = 1:length(cond)
        for fb = 1: length(freqb)
            % itpcbp
            sp1=subplot(2,2,2*c-1);
            [nelements,centers]= hist(squeeze(tar(t).con(c).concat.max_itpcbp_time(fb,:)),bins);
            plot(centers,nelements,'-','Color',colsbp(fb,:));
            hold on
            title([' max itpc-bp in ', strrep(cond{c},'_',' '),' for ',num2str(size(tar(t).con(c).concat.max_itpcbp,2)),...
                ' sites,',num2str(tar(t).con(c).concat.nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
            xlabel('Max ITPC times','Interpreter','latex');
            set(gca,'xlim',[-0.25,0.25])
            
            % powbp
            sp2=subplot(2,2,2*c);
            [nelements,centers]= hist(squeeze(tar(t).con(c).concat.max_powbp_time(fb,:)),bins);
            plot(centers,nelements,'-','Color',colsbp(fb,:));
            hold on
            title([' max power-bp in ', strrep(cond{c},'_',' '),' for ',num2str(size(tar(t).con(c).concat.max_powbp,2)),...
                ' sites,',num2str(tar(t).con(c).concat.nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
            xlabel('Max POWER times','Interpreter','latex');
            set(gca,'xlim',[-0.25,0.25])
            
        end
        legend(freqb,'FontSize',5)
        legend('boxoff')
        
        plot(sp1,[0,0],get(sp1,'ylim'),'k--')
        plot(sp2,[0,0],get(sp2,'ylim'),'k--')
    end
    mtit(['Target=',strrep(tar(t).target,'_','-')],'xoff', 0, 'yoff', 0.05, 'Color','red', 'fontsize', 12);
    results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{t},'-','Time_of_Max_ITPCbp_POWbp_for ',num2str(tar(t).nSites),' sites ', withunits]);
    export_fig(h,[results_file,'.pdf']);
%     saveas(h,[results_file]);
end

close all,
clc


%<<<<<<< HEAD
% % ploting number of significant sites in each bin for ITPC/POW/evoked lfp:
% for e=1:size(cfg.analyse_states,1)
%     h = figure;
%     E=cfg.analyse_states{e,1};
%     for t = 1: length(targets)
%         if tar(t).nSites == 0
%             continue;
%         end
%         
%         for c = 1:length(cond)   
%             bins=tar(t).con(c).(E).tfr_time;
%             nsites=num2str(tar(t).nSites);
%             ntrigs=num2str(round(tar(t).con(c).(E).nTriggers));
%             
%             % itpcbp
%             sp1=subplot(3,2,c);
%             hold on
%             title([' fraction significant itpc-bp in ', strrep(cond{c},'_',' '),' for ',nsites,...
%                 ' sites,',ntrigs,' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('Significant ITPC times','Interpreter','latex');
%             xlim([bins(1) bins(end)]);
%             
%             % powbp
%             sp2=subplot(3,2,c+2);
%             hold on
%             title([' fraction significant power-bp in ', strrep(cond{c},'_',' '),' for ',nsites,...
%                 ' sites,',ntrigs,' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('Significant POWER times','Interpreter','latex');
%             xlim([bins(1) bins(end)]);
%             
%             for fb = 1: length(freqb)
%                 plot(sp1,bins,squeeze(tar(t).con(c).(E).itpcbp_sig(fb,:)),'-','Color',colsbp(fb,:));
%                 plot(sp2,bins,squeeze(tar(t).con(c).(E).powbp_sig(fb,:)),'-','Color',colsbp(fb,:));
%             end   
%             plot(sp1,[0,0],get(sp1,'ylim'),'k--')
%             plot(sp2,[0,0],get(sp2,'ylim'),'k--')         
%             legend(freqb,'FontSize',5)
%             legend('boxoff')
%             
%             % evoked lfp
%             sp3=subplot(3,2,c+4);
%             plot(bins,tar(t).con(c).(E).lfp_sig,'-','Color',[0 0 1]);
%             hold on
%             title([' fraction significant evoked lfp in ', strrep(cond{c},'_',' '),' for ',nsites,...
%                 ' sites,',ntrigs,' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('Significant evoked LFP times','Interpreter','latex');
%             xlim([bins(1) bins(end)]);
%         end
%         
%         %% adjust y scales so that they are the same across conditions        
%         mtit(['Target=',strrep(tar(t).target,'_','-'), ', ' E],'Interpreter','latex','Color','red', 'fontsize', 12);
%         %sgtitle(['Target=',strrep(tar(t).target,'_','-'), ' ,' E],'interpreter','none','Color','red', 'fontsize', 12);
%         results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{t},' ',E,'-','Significant_bins_for ',nsites,' sites ', withunits]);
%         export_fig(h,[results_file,'.pdf']);
%     end
% end
% close all,
% clc
% =======
% % ploting number of significant sites in each bin for ITPC/POW/evoked lfp:
% bins=tfr_time;
% for tr = 1: length(targets)
%     if grand_avg(tr).nSites == 0
%         continue;
%     end
%     h = figure;
%     for c = 1:length(cond)
%         for fb = 1: length(freqb)
%             % itpcbp
%             sp1=subplot(3,2,c);
%             plot(bins,squeeze(grand_avg(tr).avg(c).itpcbpsig(fb,:)),'-','Color',color(fb,:));
%             hold on
%             title([' fraction significant itpc-bp in ', strrep(cond{c},'_',' '),' for ',num2str(size(grand_avg(tr).avg(c).max_itpcbp,2)),...
%                 ' sites,',num2str(grand_avg(tr).avg(c).nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('Significant ITPC times','Interpreter','latex');
%             set(gca,'xlim',[-0.25,0.25])
% %             xlim([bins(1) bins(end)]);
%             
%             % powbp
%             sp2=subplot(3,2,c+2);
%             plot(bins,squeeze(grand_avg(tr).avg(c).powbpsig(fb,:)),'-','Color',color(fb,:));
%             hold on
%             title([' fraction significant power-bp in ', strrep(cond{c},'_',' '),' for ',num2str(size(grand_avg(tr).avg(c).max_powbp,2)),...
%                 ' sites,',num2str(grand_avg(tr).avg(c).nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%             xlabel('Significant POWER times','Interpreter','latex');
%             set(gca,'xlim',[-0.25,0.25])
% %             xlim([bins(1) bins(end)]);            
%         end
%         
%         legend(freqb,'FontSize',5)
%         legend('boxoff')
%         
%         % evoked lfp
%         sp3=subplot(3,2,c+4);
%         plot(bins,grand_avg(tr).avg(c).lfpsig,'-','Color',[0 0 1]);
%         hold on
%         title([' fraction significant evoked lfp in ', strrep(cond{c},'_',' '),' for ',num2str(size(grand_avg(tr).avg(c).lfpsig,2)),...
%             ' sites,',num2str(grand_avg(tr).avg(c).nTriggers),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
%         xlabel('Significant evoked LFP times','Interpreter','latex');
%         set(gca,'xlim',[-0.25,0.25])
% %         xlim([bins(1) bins(end)]);
%         
%         
%         plot(sp1,[0,0],get(sp1,'ylim'),'k--')
%         plot(sp2,[0,0],get(sp2,'ylim'),'k--')
%     end
%     mtit(['Target=',strrep(grand_avg(tr).target,'_','-')],'xoff', 0, 'yoff', 0.05, 'Color','red', 'fontsize', 12);
%     results_file = fullfile(cfg.analyse_lfp_folder, [cfg.monkey,'-',targets{tr},' - ','Significant_bins_for ',num2str(grand_avg(tr).nSites),' sites ', withunits]);
%     export_fig(h,[results_file,'.pdf']);
% %     saveas(h,[results_file]);
% end
% 
% close all,
% clc
% >>>>>>> 86fc46a25bdc94ee94288cf71ea552db2eba85fa


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
