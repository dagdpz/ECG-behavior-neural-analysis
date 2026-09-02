function ecg_bna_analysis_main(project,versions)
%ecg_bna_analysis_main('Pulv_bodysignal',{'ver_LS2'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for LFP time frequency analysis
% Runs functions for reading processed LFP data, rejection of noise trials
% and task specific analysis using TFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% INITIALIZATION
% loop through settings file
% addpath(genpath_exclude('/home/shamim/fileserver/Projects/_Shamim/rework_LS_SS','.git'));
root_drive = 'Y:';%'/home/shamim/fileserver';

cfg = [];
cfg.process_MUA=0;
cfg.project = project;
cfg.fldr.root     = [root_drive,filesep,'Projects',filesep,cfg.project];
ecg_bna_location     =which('ecg_bna_define_folders');
github_folder        =ecg_bna_location(1:strfind(ecg_bna_location,['ECG-behavior-neural-analysis' filesep 'ecg_bna_define_folders'])-1);

sanitycheck=0;
for v = 1:length(versions)
    cfg.version = versions{v};
    run([github_folder filesep 'Settings' filesep cfg.project filesep 'ECG_bna' filesep cfg.version '.m']);
    cfg = ecg_bna_define_folders(cfg);
    
    %% Get info about sessions to be analysed
    sessions_info = cfg.session_info;
    
    %% per session processing..
    if cfg.process_per_session
        for i = 1:length(sessions_info)
            % First make a seed for each session -> that way even if the
            % new sessios are added in a subsequent re-running, already
            % processed sessions are unaffected
            seed_filename=[cfg.fldr.seeds filesep 'seed_' sessions_info(i).Date '.mat']; %% make seeds subfolder and add one per session
            if exist(seed_filename,'file')
                load(seed_filename);
                rng(seed);
            else
                seed=rng;
                save(seed_filename,'seed');
            end
            
            %java.lang.System.gc() % added by Luba to control the number of graphical device interface handles (hope to handle the problem with freezing plots while creating figures)
            
            load(sessions_info(i).Input_trials);
            
            % reading in actual TDT clock block starts (in seconds - inprecise, but that is irrelevant)
            blocks=unique([trials.block]);
            blockstart=ecg_bna_get_anchor_times(sessions_info(i).Monkey,sessions_info(i).Date,blocks,root_drive);
            
            cfg.event_types=cfg.analyse_states(:,2);
            cfg.events=cfg.analyse_states(:,1);
            for e=1:numel(cfg.events)
                current_event=cfg.analyse_states(e,:);
                Event=ecg_bna_get_events(sessions_info(i),trials,current_event,blockstart);
                Triggers.(cfg.events{e})=ecg_bna_jitter(Event,cfg.spk);
            end
            
            
            if cfg.process_spikes
                cfg.Input_WC=sessions_info(i).Input_WC;
                
                %% apply exclusion criteria and save lists of units - we do it once
                if cfg.spk.compute_unit_subsets
                    ecg_bna_get_unit_list(cfg,1);
                end
                %% copy selected units separately - we do it once
                if cfg.spk.move_files
                    ecg_bna_copy_selected_units(cfg)
                end
                
                %% do ECG spike analysis and computations related to cardioballistic effect
                if cfg.spk.compute_spike_histograms || cfg.spk.compute_spike_phase
                    if sanitycheck
                        load(strrep(sessions_info(i).Input_spikes,'population','by_block'));
                        ecg_bna_sanity_check(by_block,trials,Triggers,blockstart,cfg)
                    end
                    load(sessions_info(i).Input_spikes);
                    
                    ecg_bna_compute_session_phase_response_analysis(trials,population,Rpeaks,cfg)
                    
                    if cfg.spk.compute_spike_histograms
                        ecg_bna_PSTH(trials,population,Triggers,blockstart,cfg)
                    end
                    
                    if cfg.spk.compute_spike_phase
                        %ecg_bna_compute_session_ECG_related_spikePhase(trials,population,Rpeaks,sessions_info(i),cfg)
                        ecg_bna_compute_PEPH(trials,population,Triggers,blockstart,cfg)
                    end
                    
                    if cfg.spk.compute_correlation
                        ecg_bna_compute_session_correlation_analysis(trials,population,Rpeaks,cfg)
                    end
                end
                
                if cfg.spk.plot_spike_histograms
                    EPO=ecg_bna_get_plotoptions(cfg,'spk');
                    ecg_bna_plot_PSTH(sessions_info(i),EPO,cfg);
                    %ecg_bna_plot_session_spike_histogram(sessions_info(i),cfg);
                end
                
                
                if cfg.spk.plot_spike_phase
                    ecg_bna_plot_session_ECG_related_spikePhase(sessions_info(i),cfg)
                end
                if cfg.spk.plot_correlation
                    ecg_bna_plot_session_correlation(sessions_info(i),cfg)
                end
                %                 if cfg.spk.plot_spike_phase
                %                     ecg_bna_plot_session_ECG_related_spikePhase(sessions_info(i),cfg)
                %                 end
            end
            if cfg.process_MUA || cfg.process_LFP
                fprintf('Analysing for session %s\n', [sessions_info(i).Monkey '_' sessions_info(i).Date]);
                
                sitesdir=fileparts(sessions_info(i).Input_LFP{:});
                [sitefiles]=dir(sessions_info(i).Input_LFP{:});
                
                % excluding sites and ICA correction
                if cfg.process_LFP && (cfg.lfp.Reref==1 || cfg.lfp.runICA==1)
                    sr=unique([trials.TDT_LFPx_SR]);
                    cfg.LFP_ts_original=1/sr;
                    allSitesData = ecg_bna_remove_rawLFP_outliers(sitesdir,sitefiles,cfg,cfg.LFP_ts_original,sr,Triggers,blocks,blockstart);
                    removed_sites = find(~ismember({sitefiles.name},{allSitesData.name}));
                    if ~isempty(removed_sites)
                        removed_sites_name = sitefiles(removed_sites).name;
                        fprintf('\n\n removed_sites: %s\n',removed_sites_name);
                    else
                        fprintf('\n\n*****  NO site was removed ! *****\n\n');
                    end
                    sitefiles = allSitesData;
                    clear allSitesData;
                    %                 elseif (isfield(cfg.lfp, 'runICA') && cfg.lfp.runICA==0)
                    
                    %                     smplTblfilename = fullfile([cfg.analyse_lfp_folder,filesep,...
                    %                         cfg.session_info(i).Monkey(1:3),'_', cfg.session_info(i).Date , '_sampleTble_allSitesdata']);
                    %                     load([smplTblfilename,'.mat']);
                    %
                    %                     sitefiles = allSitesData;
                    %                     clear allSitesData;
                end
                
                for s = 1:length(sitefiles) %% loop only through valid sites
                    if cfg.lfp.Reref==1 || cfg.lfp.runICA==1
                        sites = sitefiles(s).site;
                    else
                        load([sitesdir filesep sitefiles(s).name], 'sites');
                    end
                    
                    if cfg.process_MUA
                        cfg.to_trigger='MUA';
                        sr=unique([trials.TDT_MUAx_SR]);
                        cfg.ts_original=1/sr;
                        site_prep    = ecg_bna_resampled_spectra(sites,cfg);
                        n_samples_per_block=site_prep.mua.n_samples_per_block;
                        site_triggers     = ecg_bna_resample_triggers(Triggers,[blocks;blockstart(blocks)],n_samples_per_block,site_prep.mua.sr);
                        site_triggered                 = ecg_bna_compute_triggered_variables(site_prep,site_triggers,trials,cfg);
                        triggered_session_MUA.sites(s) = site_triggered;
                        triggered_session_MUA.session  = site_triggered.session;
                    end
                    
                    if cfg.process_LFP
                        cfg.to_trigger='LFP';
                        sr=unique([trials.TDT_LFPx_SR]);
                        cfg.ts_original=1/sr;
                        site_prep    = ecg_bna_resampled_spectra(sites,cfg);
                        n_samples_per_block=site_prep.lfp.n_samples_per_block;
                        site_triggers     = ecg_bna_resample_triggers(Triggers,[blocks;blockstart(blocks)],n_samples_per_block,site_prep.lfp.sr);
                        site_triggered                 = ecg_bna_compute_triggered_variables(site_prep,site_triggers,trials,cfg);
                        triggered_session_LFP.sites(s) = site_triggered;
                        triggered_session_LFP.session  = site_triggered.session;
                    end
                    
                end
                if cfg.process_MUA
                    triggered_session_data=triggered_session_MUA;
                    save(fullfile(cfg.fldr.LFP_sessions, ['Triggered_session_' triggered_session_MUA.session '.mat']), 'triggered_session_data');
                    clear triggered_session_MUA triggered_session_data;
                end
                
                if cfg.process_LFP
                    triggered_session_data=triggered_session_LFP;
                    save(fullfile(cfg.fldr.LFP_sessions, ['Triggered_session_' triggered_session_LFP.session '.mat']), 'triggered_session_data');
                    clear triggered_session_LFP triggered_session_data;
                end
                clear site_prep site_triggers site_triggered;
            end
        end
    end
    
    %% average across sessions
    if cfg.process_population
        monkeys = unique({cfg.session_info.Monkey});
        cfg.monkey = [monkeys{:}];
        
%        keys=ecg_bna_get_unit_list(cfg,0);
 %       cfg.site_IDS=keys.tuning_table(2:end,find_column_index(keys.tuning_table,'site_ID'));
        
        if cfg.process_MUA
            cfg.site_IDS = ecg_bna_get_R_peak_units(cfg);
            %ecg_bna_compute_grand_avg_mua(cfg,'all');
            ecg_bna_compute_grand_avg_mua(cfg,'w_units');
            ecg_bna_compute_grand_avg_mua(cfg,'wo_units');
            
            %             load('Y:\Projects\Pulv_bodysignal\LFP\Figures\cueresponsivesites.mat');
            %             cfg.site_IDS=cueresponsivesites;
            %             grand_avg = ecg_bna_compute_grand_avg_mua(cfg,'w_units');
        end
        
        if cfg.process_LFP
            %cfg.site_IDS = ecg_bna_get_R_peak_units(cfg);
            %ecg_bna_compute_grand_avg(cfg,'all');
            cfg.combine_monkeys=0;
            ecg_bna_population_comparisons(cfg);
            cfg.combine_monkeys=1;
            ecg_bna_population_comparisons(cfg);
%             ecg_bna_compute_grand_avg(cfg,'w_units');
%             ecg_bna_compute_grand_avg(cfg,'wo_units');
            
            %             load('Y:\Projects\Pulv_bodysignal\LFP\Figures\cueresponsivesites.mat');
            %             cfg.site_IDS=cueresponsivesites;
            %             grand_avg = ecg_bna_compute_grand_avg(cfg,'w_units');
        end
        
        if cfg.process_spikes
            SPK_PSTH=load_stuff(sessions_info,cfg,'SPK_root_fldr','','per_unit','Output');
            %ecg_bna_avg_spike_histogram(SPK_PSTH,sessions_info, cfg);
            ecg_bna_avg_spike_histogram_clean(SPK_PSTH,cfg);
        end
    end
end
end

function Out = load_stuff(sessions_info,cfg,subfolder,namepart,per,varname)
out_idx=0;
for i = 1:length(sessions_info)
    monkey=sessions_info(i).Monkey(1:3);
    
    results_folder = fullfile(cfg.(subfolder),per);
    file=dir(fullfile(results_folder, [namepart monkey '_' sessions_info(i).Date '*.mat']));
    for f=1:numel(file)
        out_idx=out_idx+1;
        to_load = load([results_folder filesep file(f).name],varname);
        Out(out_idx)=to_load.(varname);
    end
end
end

