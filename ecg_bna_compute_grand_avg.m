function grand_avg = ecg_bna_compute_grand_avg(ecg_bna_cfg)


data_path = ecg_bna_cfg.sites_lfp_fldr;
cd(data_path)
% read the files
all_lfp_data = dir('*.mat');
%%
% storing all the Triggered parameters of the sites in the results folder:
out = struct;
for i = 1:length(all_lfp_data)
    load(all_lfp_data(i).name)
    session = triggered_site_data.session;
    site_ID = triggered_site_data.site_ID;
    target = triggered_site_data.target;
    
    for cn = 1:length(triggered_site_data.condition)
        if isempty(triggered_site_data.condition(cn).event)
            continue;
        end
        condition = triggered_site_data.condition(cn).label;
        nTriggers = triggered_site_data.condition(cn).event.real.ntriggers;
        lfp = squeeze(triggered_site_data.condition(cn).event.normalized.lfp.mean)';
        itpc = squeeze(triggered_site_data.condition(cn).event.real.itpc.mean)-squeeze(triggered_site_data.condition(cn).event.shuffled.itpc.mean);
        power = squeeze(triggered_site_data.condition(cn).event.normalized.pow.mean);
        itpcbp = squeeze(triggered_site_data.condition(cn).event.real.itpcbp.mean)-squeeze(triggered_site_data.condition(cn).event.shuffled.itpcbp.mean);
        powerbp = squeeze(triggered_site_data.condition(cn).event.normalized.powbp.mean);
        time = triggered_site_data.condition(cn).event.time;
        tfr_time = triggered_site_data.condition(cn).event.tfr_time;
        
        max_itpcbp = [];
        max_powbp = [];
        for fb = 1:6
            % itpcbp
            [~,max_itpcbp_idx] = max(abs(itpcbp(fb,:)));
            max_itpcbp_time = time(max_itpcbp_idx);
            max_itpcbp = [max_itpcbp;[itpcbp(fb,max_itpcbp_idx),max_itpcbp_time]];
            % powbp
            [~,max_powbp_idx] = max(abs(powerbp(fb,:)));
            max_powbp_time = time(max_powbp_idx);
            max_powbp = [max_powbp;[powerbp(fb,max_powbp_idx),max_powbp_time]];
        end
        
        

        out(i).site_ID = site_ID;
        out(i).target = target;
        out(i).condition(cn).condition_name = condition;
        out(i).condition(cn).nTriggers = nTriggers;
        out(i).condition(cn).lfp = lfp;
        out(i).condition(cn).itpc = itpc;
        out(i).condition(cn).pow = power;
        out(i).condition(cn).itpcbp = itpcbp;
        out(i).condition(cn).powbp = powerbp;
        out(i).condition(cn).max_itpcbp = max_itpcbp;
        out(i).condition(cn).max_powbp = max_powbp;
        out(i).condition(cn).time = time;
        out(i).condition(cn).tfr_time = tfr_time;
        
    end
end
out_mask = cellfun(@isempty, {out.site_ID});
out = out(~out_mask);
%%
% Computing the Target-wise averaging of total available sites
avg = struct;
grand_avg = struct;
targets = unique({out.target});
%
for tr = 1: length(targets)
    %
    avg.lfp = [];
    avg.itpc = [];
    avg.pow = [];
    avg.itpcbp = [];
    avg.powbp = [];
    avg.max_itpcbp = [];
    avg.max_powbp = [];
    avg.nTriggers = [];
    avg.cond_name = [];
    %
    sites_for_this_target = ismember({out.target},targets{tr});
    target_sites = out(sites_for_this_target);
    cond = {'Task','Rest'};
    for cn = 1:length(cond)
        conditions = [target_sites.condition];
        mask = cellfun(@isempty, {conditions.condition_name});
        conditions = conditions(~mask);
        %         all(arrayfun(@(x) x.condition(cn).nTriggers==0, target_sites)
        for i = 1:length(conditions)
            if ~strcmp(conditions(i).condition_name,cond{cn})
                continue;
            else
                avg(cn).cond_name = cond{cn};
                avg(cn).itpc = cat(3,avg(cn).itpc,conditions(i).itpc);
                avg(cn).pow = cat(3,avg(cn).pow,conditions(i).pow);
                avg(cn).lfp = cat(3,avg(cn).lfp,conditions(i).lfp);
                avg(cn).itpcbp = cat(3,avg(cn).itpcbp,conditions(i).itpcbp);
                avg(cn).powbp = cat(3,avg(cn).powbp,conditions(i).powbp);
                avg(cn).max_itpcbp = cat(3,avg(cn).max_itpcbp,conditions(i).max_itpcbp);
                avg(cn).max_powbp = cat(3,avg(cn).max_powbp,conditions(i).max_powbp);
                avg(cn).nTriggers = [avg(cn).nTriggers,conditions(i).nTriggers];
                %
                avg(cn).itpc_avg = mean(avg(cn).itpc,3);
                avg(cn).pow_avg = mean(avg(cn).pow,3);
                avg(cn).lfp_avg = mean(avg(cn).lfp,3);
                avg(cn).itpcbp_avg = mean(avg(cn).itpcbp,3);
                avg(cn).powbp_avg = mean(avg(cn).powbp,3);
                avg(cn).nTriggers_avg = round(mean(avg(cn).nTriggers,2));
                avg(cn).time = conditions(i).time;
                avg(cn).tfr_time = conditions(i).tfr_time;
                
                
                max_itpcbp = [];
                max_powbp = [];
                for fb = 1:6
                    % itpcbp
                    [~,max_itpcbp_idx] = max(abs(avg(cn).itpcbp_avg(fb,:)));
                    max_itpcbp_time = time(max_itpcbp_idx);
                    max_itpcbp = [max_itpcbp;[avg(cn).itpcbp_avg(fb,max_itpcbp_idx),max_itpcbp_time]];
                    % powbp
                    [~,max_powbp_idx] = max(abs(avg(cn).powbp_avg(fb,:)));
                    max_powbp_time = time(max_powbp_idx);
                    max_powbp = [max_powbp;[avg(cn).powbp_avg(fb,max_powbp_idx),max_powbp_time]];
                end
                
                avg(cn).max_itpcbp_avg = max_itpcbp;
                avg(cn).max_powbp_avg = max_powbp;
            end
        end
    end
    grand_avg(tr).target = targets{tr};
    grand_avg(tr).nSites = length(target_sites);
    grand_avg(tr).avg = avg;
    clear avg

end
clear out
clear conditions
fileName = fullfile([ecg_bna_cfg.analyse_lfp_folder filesep ecg_bna_cfg.monkey,'_',ecg_bna_cfg.analyse_states{1, 2} ,'_Triggered_target_wise_Grand_grand_avg_sessions_sites','.mat']);
save(fileName,'grand_avg')

%%
% plotting the results:
cond = ecg_bna_cfg.conditionname;
freq = logspace(log10(2),log10(120),60);
frequency_bands=[2 4; 4 8; 8 14; 14 30; 30 50; 70 150];
plot_names={'POW','ITPC','Power_BP','ITPC_BP','LFP_Evoked'};

% Smoothing Kernel here:
win = 1:ecg_bna_cfg.smoothWin; win=win-(numel(win)+1)/2;
half_win = ceil(size(win,2)/2)-1;
gaussian_kernel=normpdf(win,0,numel(win)/6);
gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);

%%

for tr = 1: length(targets)
    %
    for cn = 1:length(cond)
        % create figure
        h = figure('units','normalized','position',[0 0 1 1]);
        toplot={grand_avg(tr).avg(cn).pow_avg,grand_avg(tr).avg(cn).itpc_avg};
        collim{1}=[];
        collim{2}=[];
        % =========================== Power ============================= %
        sp=1;
        subplot(3,2,sp)
        image(grand_avg(tr).avg(cn).tfr_time, 1:numel(freq), squeeze(grand_avg(tr).avg(cn).pow_avg),'CDataMapping','scaled');
        set(gca,'YDir','normal');
        line([0 0], ylim, 'color', 'k');
        
        % horizontal lines to separate frequency bands
        fbandstart = [2, 4, 8, 12, 18, 32, 80];
        fbandstart_idx = zeros(size(fbandstart));
        for f = fbandstart
            f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
            line([grand_avg(tr).avg(cn).tfr_time(1) grand_avg(tr).avg(cn).tfr_time(end)], [f_idx f_idx], 'color', 'k', 'linestyle', '--');
            fbandstart_idx(fbandstart == f) = f_idx;
        end
        set(gca,'TickDir','out')
        set(gca, 'ytick', fbandstart_idx);
        set(gca, 'yticklabel', fbandstart);
        set(gca, 'ylim', [0.5,numel(freq) + 0.5]);
        title(plot_names{sp},'fontsize',10,'interpreter','none');
        nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
        collim{sp}=[min([collim{sp}(:); nonnan(:)]) max([collim{sp}(:); nonnan(:)])];
        
        % =========================== ITPC ============================= %
        sp=2;
        subplot(3,2,sp)
        image(grand_avg(tr).avg(cn).tfr_time, 1:numel(freq), squeeze(grand_avg(tr).avg(cn).itpc_avg),'CDataMapping','scaled');
        set(gca,'YDir','normal');
        line([0 0], ylim, 'color', 'k');
        
        % horizontal lines to separate frequency bands
        fbandstart = [2, 4, 8, 12, 18, 32, 80];
        fbandstart_idx = zeros(size(fbandstart));
        for f = fbandstart
            f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
            line([grand_avg(tr).avg(cn).tfr_time(1) grand_avg(tr).avg(cn).tfr_time(end)], [f_idx f_idx], 'color', 'k', 'linestyle', '--');
            fbandstart_idx(fbandstart == f) = f_idx;
        end
        set(gca,'TickDir','out')
        set(gca, 'ytick', fbandstart_idx);
        set(gca, 'yticklabel', fbandstart);
        set(gca, 'ylim', [0.5,numel(freq) + 0.5]);
        title(plot_names{sp},'fontsize',10,'interpreter','none');
        nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
        collim{sp}=[min([collim{sp}(:); nonnan(:)]) max([collim{sp}(:); nonnan(:)])];
        
        %========================== Bandpassed POWER ==================== %
        % Smoothing of the POWbp here:
        jnk = [];
        concat_input = cat(2,(grand_avg(tr).avg(cn).powbp_avg(:,half_win:-1:1)),(grand_avg(tr).avg(cn).powbp_avg(:,:)));
        concat_input = cat(2,concat_input, (grand_avg(tr).avg(cn).powbp_avg(:,end:-1:end-half_win+1)));
            for k=1:size(grand_avg(tr).avg(cn).powbp_avg,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
        clear concat_input
        smoothed = jnk(:,half_win+1:end-half_win);
        
        sp=3;
        subplot(3,2,sp)
        hold on;
        set(gca,'ColorOrder',jet(size(grand_avg(tr).avg(cn).powbp_avg,1)));
        plot(repmat(grand_avg(tr).avg(cn).time,size(grand_avg(tr).avg(cn).powbp_avg,1),1)', squeeze(smoothed)')
        line([0 0], ylim, 'color', 'k');
        xlabel('Time(s)'); ylabel('Power (W)');
        legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
        title(plot_names{sp},'fontsize',10,'interpreter','none');
        
        %========================== Bandpassed ITPC ==================== %
        % Smoothing of the ITPCbp here:
        jnk = [];
        concat_input = cat(2,(grand_avg(tr).avg(cn).itpcbp_avg(:,half_win:-1:1)),(grand_avg(tr).avg(cn).itpcbp_avg(:,:)));
        concat_input = cat(2,concat_input, (grand_avg(tr).avg(cn).itpcbp_avg(:,end:-1:end-half_win+1)));
            for k=1:size(grand_avg(tr).avg(cn).itpcbp_avg,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
        clear concat_input
        smoothed = jnk(:,half_win+1:end-half_win);
        
        sp=4;
        subplot(3,2,sp)
        hold on;
        set(gca,'ColorOrder',jet(size(grand_avg(tr).avg(cn).itpcbp_avg,1)));
        plot(repmat(grand_avg(tr).avg(cn).time,size(grand_avg(tr).avg(cn).itpcbp_avg,1),1)', squeeze(smoothed)')
        line([0 0], ylim, 'color', 'k');
        xlabel('Time(s)'); ylabel('ITPC');
        legend({strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz')},'fontsize',3);
        title(plot_names{sp},'fontsize',10,'interpreter','none');
        
        %========================== LFP evoked Potential ==================== %
        % Smoothing of the  LFP evoked Potential here:
        jnk = [];
        concat_input = cat(2,(grand_avg(tr).avg(cn).lfp_avg(:,half_win:-1:1)),(grand_avg(tr).avg(cn).lfp_avg(:,:)));
        concat_input = cat(2,concat_input, (grand_avg(tr).avg(cn).lfp_avg(:,end:-1:end-half_win+1)));
            for k=1:size(grand_avg(tr).avg(cn).lfp_avg,1)
                jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
            end
        clear concat_input
        smoothed = jnk(:,half_win+1:end-half_win);
        
        sp=5;
        subplot(3,2,sp)
        hold on;
        set(gca,'ColorOrder',jet(size(grand_avg(tr).avg(cn).lfp_avg,1)));
        plot(repmat(grand_avg(tr).avg(cn).time,size(grand_avg(tr).avg(cn).lfp_avg,1),1)', squeeze(smoothed)')
        line([0 0], ylim, 'color', 'k');
        xlabel('Time(s)'); ylabel(' LFP evoked Potential');
        title(plot_names{sp},'fontsize',10,'interpreter','none');

        %=================================================================%
        %% format spectra colors
        cbtitle = {'(P - \mu) / std','P - \mu'};
        for sp=1:2
            subplot(3,2,sp);
            set(gca,'CLim',collim{sp})
        cm = colormap('jet');
        cb = colorbar;
        set(get(cb,'title'),'string', cbtitle{sp}, 'fontsize',8);
        colormap(cm);
        end
    
    results_file = fullfile(ecg_bna_cfg.analyse_lfp_folder, [ecg_bna_cfg.monkey,'-',targets{tr},' - ','avg of ',num2str(grand_avg(tr).nSites),' sites in all sessions - ',grand_avg(tr).avg(cn).cond_name]);
    mtit([ ecg_bna_cfg.monkey,'-',targets{tr},'-avg of ',num2str(grand_avg(tr).nSites),' sites in all sessions - ',grand_avg(tr).avg(cn).cond_name],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
    export_fig(h,[results_file,'.pdf']);
    end
end
close all,
clc

%%
freqb = {'delta 2-4 Hz','theta 4-8 Hz','alpha 8-14 Hz','beta 14-30 Hz','lowgamma 30-50 Hz','highgamma 70-150 Hz'};
freqName = {'delta','theta','alpha','beta','lowGamma','highGamma'};
color = {[0 0 1],[0 1 1],[0 1 0],[1 1 0],[1 0 1],[1 0 0]};

for tr = 1:length(targets)
    h = figure('Name',['Max ITPCbp/POWbp of Target=',grand_avg(tr).target],'NumberTitle','off');
    for cn = 1:length(cond)
        for fb = 1: length(freqb)
            % itpcbp
            subplot(2,2,2*cn-1)
            scatter(squeeze(grand_avg(tr).avg(cn).max_itpcbp(fb,2,:)),squeeze(grand_avg(tr).avg(cn).max_itpcbp(fb,1,:)),15,color{fb},'fill')
            hold on           
            title([' max itpc-bp in ', strrep(cond{cn},'_',' '),' for ',num2str(grand_avg(tr).nSites),...
                ' sites,',num2str(grand_avg(tr).avg(cn).nTriggers_avg),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
            xlabel('time (s)','Interpreter','latex'), ylabel('Max ITPC value','Interpreter','latex');
            set(gca,'xlim',[-0.5,0.5],'ylim',[0,1])
            
            % powbp
            subplot(2,2,2*cn)
            scatter(squeeze(grand_avg(tr).avg(cn).max_powbp(fb,2,:)),squeeze(grand_avg(tr).avg(cn).max_powbp(fb,1,:)),15,color{fb},'fill')
            hold on            
            title([' max power-bp in ', strrep(cond{cn},'_',' '),' for ',num2str(grand_avg(tr).nSites),...
                ' sites,',num2str(grand_avg(tr).avg(cn).nTriggers_avg),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
            xlabel('time (s)','Interpreter','latex'), ylabel('Max POWER value','Interpreter','latex');
            set(gca,'xlim',[-0.5,0.5])
           
        end
        
        for fb = 1:6
            hold on
            subplot(2,2,2*cn-1)
            plot([0,0],get(gca,'ylim'),'k--')
            %
            subplot(2,2,2*cn)
            plot([0,0],get(gca,'ylim'),'k--')
        end
        sgtitle(['Target=',strrep(grand_avg(tr).target,'_','-')],'Color','red')
        legend(freqb,'FontSize',5)
        legend('boxoff')
        
    end
    results_file = fullfile(ecg_bna_cfg.analyse_lfp_folder, [ecg_bna_cfg.monkey,'-',targets{tr},' - ','Max_ITPCbp_POWbp_of ',num2str(grand_avg(tr).nSites),' sites in all sessions']);
    %     mtit([ ecg_bna_cfg.monkey,'-',targets{tr},'-avg of ',num2str(grand_avg(tr).nSites),' sites in all sessions - ',grand_avg(tr).avg(cn).cond_name],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
    export_fig(h,[results_file,'.pdf']);
end
close all,
clc

%%
% ploting the avgogram of the timing of the Max ITPC/POW:
color = [0 0 1;0 1 1;0 1 0;1 1 0;1 0 1;1 0 0];
for tr = 1: length(targets)
    
    for cn = 1:length(cond)
        figure(tr+5)
        for fb = 1: length(freqb)
            % itpcbp
            subplot(2,2,2*cn-1)
            [nelements,centers]= hist(squeeze(grand_avg(tr).avg(cn).max_itpcbp(fb,2,:)));
            plot(centers,nelements,'o-','Color',color(fb,:)); 
            hold on           
            title([' max itpc-bp in ', strrep(cond{cn},'_',' '),' for ',num2str(size(grand_avg(tr).avg(cn).max_itpcbp,3)),...
                ' sites,',num2str(grand_avg(tr).avg(cn).nTriggers_avg),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
            xlabel('Max ITPC times','Interpreter','latex');
            plot([0,0],get(gca,'ylim'),'k--')
         
%             % powbp
            subplot(2,2,2*cn)
             [nelements,centers]= hist(squeeze(grand_avg(tr).avg(cn).max_powbp(fb,2,:)));
             plot(centers,nelements,'o-','Color',color(fb,:)); 
            hold on
            title([' max power-bp in ', strrep(cond{cn},'_',' '),' for ',num2str(size(grand_avg(tr).avg(cn).max_powbp,3)),...
                ' sites,',num2str(grand_avg(tr).avg(cn).nTriggers_avg),' avg nTriggers'],'FontSize',6,'Interpreter','latex');
            xlabel('Max POWER times','Interpreter','latex');
            plot([0,0],get(gca,'ylim'),'k--')

        end
        sgtitle(['Target = ',strrep(grand_avg(tr).target,'_','-')],'Color','red')
        legend(freqb,'FontSize',5)
        legend('boxoff')

    end
    results_file = fullfile(ecg_bna_cfg.analyse_lfp_folder, [ecg_bna_cfg.monkey,'-',targets{tr},' - ','Time_of_Max_ITPCbp_POWbp_for ',num2str(grand_avg(tr).nSites),' sites in all sessions']);
    export_fig(h,[results_file,'.pdf']);
end

close all,
clc

end