function ecg_bna_compute_ITPC_vs_POWER_Corr
close all,
clear all
clc

%%
withunits = 'all_sites';
monkey = 'Bacchus';%'Magnus'; 
targets = {'VPL', 'dPul', 'MD'};% unique({sites.target});
cond = {'Rest','Task'};
path_to_save = ['Y:\Projects\Pulv_bodysignal\LFP\ECG_',monkey,'_TaskRest_generalTrig_finalRun\not filtered\grand_average_wo combine_hemispheres_LFP real_New GrandAvg version\'];

%%
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
% ====================================================================
% read Per site data real itpc and real pow, then zscore them based on
% their own mean and std to have normalized itpc/ pow in same readable
% range, then calculate the target weise struct
% ====================================================================
data_path = ['Y:\Projects\Pulv_bodysignal\LFP\ECG_',monkey,'_TaskRest_generalTrig_finalRun\not filtered\Per_Site'];
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
            % think about event present?
            
            nTriggers = event.real.ntriggers;
            time      = event.time;
            tfr_time  = event.tfr_time;
            lfp       = squeeze(event.real.lfp.mean)';
            % lfp       = squeeze(event.real.lfp.mean)' - squeeze(event.shuffled.lfp.mean)';
            % lfp       = squeeze(event.normalized.lfp.mean)';
            %             itpc      = squeeze(event.real.itpc.mean)-squeeze(event.shuffled.itpc.mean);
            itpc      = zscore(squeeze(event.real.itpc.mean),1,'all');
            %             power     = squeeze(event.normalized.pow.mean);
            power     = zscore(squeeze(event.real.pow.mean),1,'all');
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
            E              = [];
            E.pow          = [];
            E.itpc         = [];
            E.lfp          = [];
            E.itpcbp       = [];
            E.powbp        = [];
            E.tfr_time     = [];
            E.lfp_time     = [];
            E.lfp_std      = [];
            
            E.nTriggers = mean([events(:,e).nTriggers]);
            E.pow       = (cat(3,events(:,e).pow));
            E.itpc      = (cat(3,events(:,e).itpc));
            E.lfp       = (cat(3,events(:,e).lfp));
            
            E.itpcbp    = (cat(3,events(:,e).itpcbp));
            E.powbp     = (cat(3,events(:,e).powbp));
            E.tfr_time  = events(1,e).tfr_time;
            E.lfp_time  = events(1,e).time;
            E.lfp_std   = std(cat(3,events(:,e).lfp),0,3);
            
            E.pow       = E.pow(fx:end,:,:);
            E.itpc      = E.itpc(fx:end,:,:);
            E.itpcbp    = E.itpcbp(fbx:end,:,:);
            E.powbp     = E.powbp(fbx:end,:,:);
            
            
            eventname = cfg.analyse_states{e,1};
            tar(t).con(c).(eventname)=E;
        end
        tar(t).con(c).cond_name = cond{c};
    end
    tar(t).target = targets{t};
    tar(t).nSites = length(target_sites);
end

% clear sites,
%%
tfr_freq  = freq;
tfr_time  = tar(1).con(1).R.tfr_time;
tmwnd = [0.15,0.20];% sec
frwnd = [4,8]; % freq

[~, idx1] = min(abs(tfr_time - tmwnd(1)));
[~, idx2] = min(abs(tfr_time - tmwnd(2)));
tmidx = [idx1,idx2];

[~, idx1] = min(abs(tfr_freq - frwnd(1)));
[~, idx2] = min(abs(tfr_freq - frwnd(2)));
fridx = [idx1,idx2];
%%
for e = 1%: length(cfg.analyse_states)
    eventname = cfg.analyse_states{e,1};
    for c = 1:2
        for t = 1: length(targets)
            data_itpc = tar(t).con(c).(eventname).itpc; % grand adverage of MD target, Rest Condition , across grand_avg(2).nSites;
            avg_data_itpc = squeeze(mean(mean(data_itpc(fridx(1):fridx(2),tmidx(1):tmidx(2),:),1),2));
            data_pow  = tar(t).con(c).(eventname).pow; % grand adverage of MD target, Task Condition ,across grand_avg(2).nSites;
            avg_data_pow = squeeze(mean(mean(data_pow(fridx(1):fridx(2),tmidx(1):tmidx(2),:),1),2));
            
            x = avg_data_pow;
            y = avg_data_itpc;
            [rho,pval] = corr(x,y);
            h(t)= figure;
            scatter(x, y,'k', 'filled','MarkerFaceAlpha',0.5); % Scatter plot with blue markers
            hold on;
            coeffs = polyfit(x, y, 1); % First-degree polynomial (linear fit)
            y_fit = polyval(coeffs, x); % Predicted y values based on the fit
            
            % Plot the correlation line
            plot(x, y_fit, 'r-', 'LineWidth', 2); % Red line for the fit
            xlabel("Mean Power (Z-score)")
            ylabel("Mean ITPC (Z-score)","FontSize",6)
            title([monkey, ' - ',targets{t},' - nSite: ', num2str(size(data_itpc,3)),'- ITPC vs POW correlation in ',cond{c},' - r = ', num2str(rho) , ' p = ', num2str(pval)],'FontSize',10);
            results_file = fullfile([path_to_save,filesep,monkey,'-',targets{t},'- ITPC vs POW correlation in ',cond{c}]);
            export_fig(h(t),[results_file,'.pdf']);
        end
    end
end
close all
end