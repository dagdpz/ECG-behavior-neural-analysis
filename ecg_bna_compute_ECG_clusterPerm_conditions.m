function ecg_bna_compute_ECG_clusterPerm_conditions
close all,
clear all,
clc
%%
reprosess = 0;

project = 'Pulv_bodysignal';
version = 'ECG_Magnus_TaskRest_generalTrig_ECG'; % 'ECG_Magnus_TaskRest_generalTrig_ECG' , 'ECG_Bacchus_TaskRest_generalTrig_ECG'


driver_path = 'Y:';%'/home/shamim/fileserver';

cfg = [];
cfg.ecg.timestep = 1;

cfg.process_ecg=0;
cfg.project = project;
cfg.results_folder = [driver_path,filesep,'Projects',filesep,cfg.project];
ecg_bna_location     =which('ecg_bna_define_folders');
github_folder        =ecg_bna_location(1:strfind(ecg_bna_location,['ECG-behavior-neural-analysis' filesep 'ecg_bna_define_folders'])-1);
cfg.version = version;
run([github_folder filesep 'Settings' filesep cfg.project filesep 'ECG_bna' filesep cfg.version '.m']);
cfg = ecg_bna_define_folders(cfg);

%% Get info about sessions to be analysed
sessions_info = cfg.session_info;


monkeys = unique({cfg.session_info.Monkey});
cfg.monkey = [monkeys{:}];
cfg.session_ecg_fldr = fullfile(cfg.ECG_root_results_fldr, 'Per_Session');
cfg.sites_ecg_fldr   = fullfile(cfg.ECG_root_results_fldr, 'Per_Site');

if reprosess
    
    for i = 1:length(sessions_info)
        % First make same seed for each session -> that way even if
        % the data changes in a subsequent re-running, shuffles in
        % unaffected sessions are unaffected...
        seed_filename=[cfg.ECG_root_results_fldr filesep 'seed.mat']; %% not quite sure yet where to put seed
        if exist(seed_filename,'file')
            load(seed_filename);
            rng(seed);
        else
            seed=rng;
            save(seed_filename,'seed');
        end
        
        java.lang.System.gc() % added by Luba to control the number of graphical device interface handles (hope to handle the problem with freezing plots while creating figures)
        
        load(sessions_info(i).Input_trials);
        
        % reading in actual TDT clock block starts (in seconds - inprecise, but that is irrelevant)
        blocks=unique([trials.block]);
        blockstart=ecg_bna_get_anchor_times(monkey,sessions_info(i).Date,blocks,driver_path);
        
        cfg.event_types=cfg.analyse_states(:,2);
        cfg.events=cfg.analyse_states(:,1);
        for e=1:numel(cfg.events)
            current_event=cfg.analyse_states(e,:);
            Event=ecg_bna_get_events(sessions_info(i),trials,current_event,blockstart);
            Triggers.(cfg.events{e})=ecg_bna_jitter(Event,cfg.spk);
        end
        
        fprintf('Analysing for session %s\n', [sessions_info(i).Monkey '_' sessions_info(i).Date]);
        cfg.session_ecg_fldr = fullfile(cfg.ECG_root_results_fldr, 'Per_Session');
        cfg.sites_ecg_fldr   = fullfile(cfg.ECG_root_results_fldr, 'Per_Site');
        cfg.sites_fldr       = cfg.sites_ecg_fldr;
        cfg.to_trigger   = 'ECG';
        
        %% this is new
        sitesdir=fileparts(sessions_info(i).Input_LFP{:});
        [sitefiles]=dir(sessions_info(i).Input_LFP{:});
        sr=unique([trials.TDT_ECG1_SR]);
        ts_original=1/sr;
        
        for s = 1:length(sitefiles) %% loop only through valid sites
            load([sitesdir filesep sitefiles(s).name], 'sites');
            site_ECG = ecg_bna_process_ECG(sites, cfg, ts_original);
            n_ECG_samples_per_block=site_ECG.tfs.n_samples_per_block;
            
            site_triggers     = ecg_bna_resample_triggers2(Triggers,[blocks;blockstart(blocks)],n_ECG_samples_per_block,site_ECG.tfs.sr);
            site_data=ecg_bna_compute_triggered_variables_ECG(site_ECG,site_triggers,trials,cfg );
            
            triggered_session_data.sites(s) = site_data;
            triggered_session_data.session = site_data.session;
        end
        
        
        
        
        
    end
else
    
    withunits = 'all_sites';
    path_to_save = ['Y:\Projects\Pulv_bodysignal\ECG\ECG_',monkey,'_TaskRest_generalTrig_ECG\grand_avg_all\'];
    perSitedata_path = ['Y:\Projects\Pulv_bodysignal\ECG\ECG_',monkey,'_TaskRest_generalTrig_ECG\Per_Site'];
    %%
    
    targets = {'VPL', 'dPul', 'MD'};% unique({sites.target});
    cond = {'Rest','Task'};

    % Smoothing Kernel here:
    win = 1:cfg.lfp.smoothWin; win=win-(numel(win)+1)/2;
    half_win = ceil(size(win,2)/2)-1;
    gaussian_kernel=normpdf(win,0,numel(win)/6);
    gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);

    cfg.analyse_states = {'R',    'Rpeak',1,-0.25, 0.25};
    %%
    
    cd(perSitedata_path)
    % read the files
    all_ecg_data = dir('*.mat');
    % storing all the Triggered parameters of the sites in the results folder:
    sites = struct;
    s=0;
    for f = 1:length(all_ecg_data)
        load(all_ecg_data(f).name)
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
                ecg_time  = event.tfr_time;
                %                     ecg       = squeeze(event.real.ecg.mean)'; % not normalized ecg
                ecg       = squeeze(event.normalized.ecg.mean)';
                
                
                ecg_sig   =squeeze(event.significance.ecg)'; %% weird dimension thing
                
                sites(s).condition(c).event(e).nTriggers = nTriggers;
                sites(s).condition(c).event(e).ecg = ecg;
                sites(s).condition(c).event(e).ecg_sig = ecg_sig;
                
                sites(s).condition(c).event(e).time = time;
                sites(s).condition(c).event(e).tfr_time = ecg_time;
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
                E.ecg           = [];
                E.ecg_time      = [];
                E.ecg_std       = [];
                E.ecg_25        = [];
                E.ecg_75        = [];
                E.ecg_sig       = [];
                E.ecg_sig_signed= [];
                
                E.nTriggers     = mean([events(:,e).nTriggers]);
                E.ecg           = (cat(3,events(:,e).ecg));
                E.ecg_25        = prctile(cat(3,events(:,e).ecg),25,3);
                E.ecg_75        = prctile(cat(3,events(:,e).ecg),75,3);
                E.ecg_time      = events(1,e).time;
                E.ecg_std       = std(cat(3,events(:,e).ecg),0,3);
                E.ecg_sig       = mean(cat(3,events(:,e).ecg_sig),3);
                E.ecg_sig_signed= mean(cat(3,events(:,e).ecg_sig).*sign(cat(3,events(:,e).ecg)),3);
                
                E.ecg_sig       = E.ecg_sig*100;
                E.ecg_sig_signed= E.ecg_sig_signed*100;
                
                
                eventname = cfg.analyse_states{e,1};
                tar(t).con(c).(eventname)=E;
            end
            tar(t).con(c).cond_name = cond{c};
        end
        tar(t).target = targets{t};
        tar(t).nSites = length(target_sites);
    end
    %%
    % Permutation Based on Lukas function
    numPermutation = 1000;
    pthreshold = 0.05;
    clusterthreshold = 0.05;
    % convert p threshold into a t threshold
    tThreshold = (1-pthreshold)/2;
    Fin = {'ecg'};
    for t = 1: length(targets)
        for fin = 1:length(Fin)
            for e=1:size(cfg.analyse_states,1)
                
                data_rest = squeeze(tar(t).con(1).(eventname).(Fin{fin})); % grand adverage of MD target, Rest Condition , across grand_avg(2).nSites;
                data_task = squeeze(tar(t).con(2).(eventname).(Fin{fin})); % grand adverage of MD target, Task Condition ,across grand_avg(2).nSites;
                ecg_time = tar(t).con(1).(eventname).ecg_time;
                
                
                rest_R = permute(data_rest, [2, 1]);  % [trials x  time ]
                task_R = permute(data_task, [2, 1]);  % [trials x  time ]
                
                
                % generating the shuffled files :
                tic
                task_shf = zeros(numPermutation, size(task_R,2));
                rest_shf = zeros(numPermutation, size(rest_R,2));
                for s = 1:numPermutation
                    inverseindexes = (randi(2, size(task_R,1), 1) - 1.5) * 2;  % [trials x 1]
                    task_shf(s,:) = mean(task_R .* inverseindexes, 1);
                    rest_shf(s,:) = mean(rest_R .* inverseindexes, 1);
                end
                toc
                
                tic
                shf=task_shf-rest_shf;
                s_ix=1:numPermutation;
                for s = 1:numPermutation
                    T = sum(shf(s_ix~=s,:) < repmat(shf(s,:), [numPermutation-1,1]), 1) / (numPermutation-1) - 0.5;
                    
                    sigpos = T > tThreshold;
                    signeg = T < -tThreshold;
                    
                    CCpos = bwconncomp(sigpos);  % 1D cluster
                    CCneg = bwconncomp(signeg);
                    
                    sumTpos = cellfun(@(x) sum(T(x)), CCpos.PixelIdxList);
                    sumTneg = cellfun(@(x) sum(T(x)), CCneg.PixelIdxList);
                    
                    T_pos(s) = max([sumTpos 0]);
                    T_neg(s) = min([sumTneg 0]);
                end
                toc
                
                
                %% now find significant clusters in the real data
                T_max_thr=prctile(T_pos,100*(1-clusterthreshold));   % sum of
                T_min_thr=prctile(T_neg,100*clusterthreshold);    % sum of
                
                T=squeeze(sum(shf<repmat(mean((task_R-rest_R),1),[size(shf,1),1,1]),1));
                T=T/(size(shf,1)-1)-0.5;
                
                sigpos=T>tThreshold;
                signeg=T<-tThreshold;
                
                CCpos = bwconncomp(sigpos);
                CCneg = bwconncomp(signeg);
                sigCpos = cellfun(@(x) sum(T(x))>=T_max_thr,CCpos.PixelIdxList);
                sigCneg = cellfun(@(x) sum(T(x))<=T_min_thr,CCneg.PixelIdxList);
                sigpixpos=vertcat(CCpos.PixelIdxList{sigCpos});
                sigpixneg=vertcat(CCneg.PixelIdxList{sigCneg});
                significance_pos=zeros(size(T));
                significance_neg=zeros(size(T));
                significance_pos(sigpixpos)=true;
                significance_neg(sigpixneg)=true;
                
                
                % ==============================================================
                % plotting the significance
                h(e) = figure;
                h(e).WindowState = 'maximized';
                
                % Plotting the Task and Rest data
                linePropsT ={'color',[1 0 0]};
                linePropsR ={'color',[0 0 1]};
                
                task_R_std = squeeze(tar(t).con(2).(eventname).([Fin{fin},'_std']));
                rest_R_std = squeeze(tar(t).con(1).(eventname).([Fin{fin},'_std']));
                
                smoothed_meanT = smoothit(mean(task_R,1),half_win,gaussian_kernel);
                smoothed_stdT = smoothit(task_R_std,half_win,gaussian_kernel);
                
                smoothed_meanR = smoothit(mean(rest_R,1),half_win,gaussian_kernel);
                smoothed_stdR = smoothit(rest_R_std,half_win,gaussian_kernel);
                
                shadedErrorBar(ecg_time,smoothed_meanT,smoothed_stdT,linePropsT,1); hold on;
                shadedErrorBar(ecg_time,smoothed_meanR,smoothed_stdR,linePropsR,1);
                
% %                 if strcmp(Fin{fin},'ecg') && strcmp(monkey,'Magnus')
% %                     y1 = [-12 8];
% %                 elseif strcmp(Fin{fin},'ecg') && strcmp(monkey,'Bacchus')
% %                     y1 = [-5 10];
% %                 end
                y1 = ylim;
                for i = 1:CCpos.NumObjects
                    idx = CCpos.PixelIdxList{i};
                    t_start = ecg_time(idx(1));
                    t_end = ecg_time(idx(end));
                    
                    % Draw a semi-transparent rectangle over the significant cluster
                    patch([t_start t_end t_end t_start], ...
                        [y1(1) y1(1) y1(2) y1(2)], ...
                        [0.1 0 0.1], ...       % [0.8 0.8 0.8], ...       % gray color
                        'FaceAlpha', 0.3, ...    % transparency
                        'EdgeColor', 'none');
                end
                
                for i = 1:CCneg.NumObjects
                    idx = CCneg.PixelIdxList{i};
                    t_start = ecg_time(idx(1));
                    t_end = ecg_time(idx(end));
                    
                    % Draw a semi-transparent rectangle over the significant cluster
                    patch([t_start t_end t_end t_start], ...
                        [y1(1) y1(1) y1(2) y1(2)], ...
                        [0 0.1 0.1], ...       % [0.8 0.8 0.8], ...       % gray color
                        'FaceAlpha', 0.3, ...    % transparency
                        'EdgeColor', 'none');
                end
                axis square;
                set(gca,'ylim',y1);
                xlabel('Time(s)');
                line([0 0], y1, 'color', 'k');
                title([monkey,'-',targets{t},'-',Fin{fin},' - nSite: ', num2str(size(data_task,2)),'-',cond{2},'-',cond{1},'- Significant difference - numPerm =',...
                    num2str(numPermutation)],'fontsize', 12,'Interpreter', 'none');
                
                ylabel('ecg evoked Potential (A.U.)');
                results_file = fullfile([path_to_save,filesep,monkey,'-',targets{t},'-',Fin{fin},'-',cond{2},' vs. ',cond{1},'_Significant difference']);
                
                
                export_fig(h(e),[results_file,'.pdf']);
            end
        end
    end
%     close all
end

end


function out=smoothit(in,half_win,gaussian_kernel)
concat_input = cat(2,(in(:,half_win:-1:1)),in);
concat_input = cat(2,concat_input, (in(:,end:-1:end-half_win+1)));
for k=1:size(concat_input,1)
    jnk(k,:)= conv(squeeze(concat_input(k,:)),gaussian_kernel,'same');
end
out = jnk(:,half_win+1:end-half_win);
end
