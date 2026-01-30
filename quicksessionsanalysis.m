
monkeys={'Bacchus','Magnus'};
for m=1:numel(monkeys)
monkey=monkeys{m};

% reading all the electrode depth and store it perSite:
if exist(['Y:\Projects\Pulv_bodysignal\MUA\ECG_',monkey,'_TaskRest_generalTrig_MUA',filesep,'allSites_info.mat'],'file')
    load(['Y:\Projects\Pulv_bodysignal\MUA\ECG_',monkey,'_TaskRest_generalTrig_MUA',filesep,'allSites_info.mat']);
else
    % this part not so clear to me
    epbys_path=['Y:\Projects\Pulv_bodysignal\ephys\ECG_TaskRest_' monkey '_MUA'];
    all_perSite_raw_data = dir([epbys_path filesep 'sites_',monkey,'_*.mat']);
    for si = 1:length(all_perSite_raw_data)
        load([epbys_path filesep all_perSite_raw_data(si).name])
        siteFields = fieldnames(sites);
        site_info(si) = rmfield(sites,{siteFields{8:end-2}}); %% remove everuything exepth required info
    end
    clear sites;
    save(fullfile(['Y:\Projects\Pulv_bodysignal\MUA\ECG_',monkey,'_TaskRest_generalTrig_MUA',filesep,'allSites_info.mat']),'site_info');
end


% prepare sites, remove wrong targets

load(['Y:\Projects\Pulv_bodysignal\LFP\ECG_' monkey '_unfiltered\' monkey '_Rpeak_Triggered_target_wise_Grand_grand_avg_sessions_sitesall.mat'])

val_targets={'dPul_L','dPul_R','MD_L','MD_R','VPL_L','VPL_R'};
targets={sites.target};
ti=ismember(targets,val_targets);
sites=sites(ti);

for s =1:numel(sites)
    sites(s).session=sites(s).site_ID(5:12);
end


%% add x/y/z to sites:
for s=1:numel(sites)
    site_ID=sites(s).site_ID;
    ix=ismember({site_info.site_ID},site_ID);
    sites(s).x=site_info(ix).grid_x; %*0.8;
    sites(s).y=site_info(ix).grid_y; %*0.8;
    sites(s).z=round(site_info(ix).electrode_depth*10)/10;
end

siteinfo=rmfield(sites,'condition');
save(['Y:\Projects\Pulv_bodysignal\LFP\ECG_' monkey '_unfiltered\' monkey '_siteinfo.mat'],'siteinfo')
sessions=unique({sites.session});
end












load('Y:\Projects\Pulv_bodysignal\LFP\ECG_Magnus_unfiltered\Magnus_Rpeak_Triggered_target_wise_Grand_grand_avg_sessions_sitesall.mat');

site_IDs.M={sites.site_ID};
e=1;
for s=1:numel(sites)
 con1.M(s)=max(max(sites(s).condition(1).event(e).itpc));
 con2.M(s)=max(max(sites(s).condition(2).event(e).itpc));
 ntrigs1.M(s)=max(max(sites(s).condition(1).event(e).nTriggers));
 ntrigs2.M(s)=max(max(sites(s).condition(2).event(e).nTriggers));
end


load('Y:\Projects\Pulv_bodysignal\LFP\ECG_Bacchus_unfiltered\Bacchus_Rpeak_Triggered_target_wise_Grand_grand_avg_sessions_sitesall.mat');

site_IDs.B={sites.site_ID};
e=1;
for s=1:numel(sites)
 con1.B(s)=max(max(sites(s).condition(1).event(e).itpc));
 con2.B(s)=max(max(sites(s).condition(2).event(e).itpc));
 ntrigs1.B(s)=max(max(sites(s).condition(1).event(e).nTriggers));
 ntrigs2.B(s)=max(max(sites(s).condition(2).event(e).nTriggers));
end

bins=1:1000:20000;
bins=1:1000:20000;
figure
subplot(2,2,1)
hist(ntrigs1.M,bins)
title('Magnus rest')
xlabel('N triggers')
ylabel('N sites')
xlim([bins(1) bins(end)])
subplot(2,2,2)
hist(ntrigs2.M,bins)
title('Magnus task')
xlabel('N triggers')
ylabel('N sites')
xlim([bins(1) bins(end)])
subplot(2,2,3)
hist(ntrigs1.B,bins)
title('Bacchus rest')
xlabel('N triggers')
ylabel('N sites')
xlim([bins(1) bins(end)])
subplot(2,2,4)
hist(ntrigs2.B,bins)
title('Bacchus task')
xlabel('N triggers')
ylabel('N sites')
xlim([bins(1) bins(end)])