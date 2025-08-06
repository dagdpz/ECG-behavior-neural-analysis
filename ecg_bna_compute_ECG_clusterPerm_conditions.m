function ecg_bna_compute_ECG_clusterPerm_conditions
close all,
clear all,
clc
%%

project = 'Pulv_bodysignal';
version = 'ECG_Bacchus_TaskRest_generalTrig_ECG'; % 'ECG_Magnus_TaskRest_generalTrig_ECG' , 'ECG_Bacchus_TaskRest_generalTrig_ECG'


driver_path = 'Y:';%'/home/shamim/fileserver';

cfg = [];
cfg.ecg.timestep = 1;
 
cfg.process_MUA=0;
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
end
