function ecg_bna_analysis_main(project,versions)
%ecg_bna_analysis_main('Pulv_bodysignal',{'ver_LS2'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for LFP time frequency analysis
% Runs functions for reading processed LFP data, rejection of noise trials
% and task specific analysis using TFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% whether the LFP should be processed (true) or not (false)
% if the LFP for the sessions to analyse has already been processed, and no
% settings need to be changed, this flag can be set to false to skip LFP
% being processed again
% If LFP was not previously processed and the flag is set to false,
% analysis won't happen
% TODO: check if LFP is processed, of not, process LFP even if flag is set
% to false

%% INITIALIZATION
% loop through settings file

cfg = [];
cfg.project = project;
cfg.results_folder = ['Y:\Projects\' cfg.project];
ecg_bna_location     =which('ecg_bna_define_folders');
github_folder        =ecg_bna_location(1:strfind(ecg_bna_location,['ECG-behavior-neural-analysis' filesep 'ecg_bna_define_folders'])-1);

for v = 1:length(versions)
    cfg.version = versions{v};
    run([github_folder filesep 'Settings' filesep cfg.project filesep 'ECG_bna' filesep cfg.version '.m']);
    cfg = ecg_bna_define_folders(cfg);
    
    %% Get info about sessions to be analysed
    % Read the info about sessions to analyse
    sessions_info = cfg.session_info;
    
    %% per session processing..
    if cfg.process_per_session
        for i = 1:length(sessions_info)
            java.lang.System.gc() % added by Luba to control the number of graphical device interface handles (hope to handle the problem with freezing plots while creating figures)
            
            session_name = [sessions_info(i).Monkey '_' sessions_info(i).Date];
            load(sessions_info(i).Input_trials);
            
            % First make seed and ecg shuffles, then use those shuffles for all subfunctions
            seed_filename=[cfg.ECG_root_results_fldr filesep 'seed.mat']; %% not quite sure yet where to put seed
            if exist(seed_filename,'file')
                load(seed_filename);
                rng(seed);
            else
                seed=rng;
                save(seed_filename,'seed');
            end
            
            
            if cfg.process_spikes
                
                
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
                    Rpeaks=ecg_bna_compute_session_shuffled_Rpeaks(sessions_info(i),cfg.spk);
                    %Rpeaks=ecg_bna_jitter(sessions_info(i),cfg.spk);
                    cfg.Input_WC=sessions_info(i).Input_WC;
                    load(sessions_info(i).Input_spikes);
                    if cfg.spk.compute_spike_histograms
                        ecg_bna_compute_session_spike_histogram(trials,population,Rpeaks,cfg);
                    end
                    
                    if cfg.spk.compute_spike_phase
                        ecg_bna_compute_session_ECG_related_spikePhase(trials,population,Rpeaks,cfg)
                    end
                end
                
                if cfg.spk.plot_spike_histograms
                    ecg_bna_plot_session_spike_histogram(sessions_info(i),cfg);
                end
                if cfg.spk.plot_spike_phase
                    ecg_bna_plot_session_ECG_related_spikePhase(sessions_info(i),cfg)
                end
                aa=1;
            end
            if cfg.process_LFP
                
                Rpeaks=ecg_bna_compute_session_shuffled_Rpeaks(sessions_info(i),cfg.lfp);
                %Rpeaks=ecg_bna_jitter(sessions_info(i),cfg.lfp);
                
                fprintf('Analysing for session %s\n', session_name);
                
                % Read LFP data
                cfg.session_lfp_fldr = fullfile(cfg.analyse_lfp_folder, 'Per_Session');
                cfg.sites_lfp_fldr   = fullfile(cfg.analyse_lfp_folder, 'Per_Site');
                
                %% this is new
                sitesdir=fileparts(sessions_info(i).Input_LFP{:});
                [sitefiles]=dir(sessions_info(i).Input_LFP{:});
                sr=unique([trials.TDT_LFPx_SR]);
                ts_original=1/sr;
                
                for s = 1:length(sitefiles)
                    load([sitesdir filesep sitefiles(s).name], 'sites');
                    site_LFP= ecg_bna_process_LFP(sites, cfg, ts_original);
                    n_LFP_samples_per_block=site_LFP.tfs.n_samples_per_block;
                    
                    
                    %% CHECK BLOCKS OVERLAP - remove LFP for which there is no trigger block - we did the opposite already in resample_triggers!
                    blocks_not_present_in_triggers=n_LFP_samples_per_block(1, ~ismember(n_LFP_samples_per_block(1,:),[Rpeaks.block]));
                    if numel(blocks_not_present_in_triggers)==size(n_LFP_samples_per_block,2)
                        continue;
                    end
                    
                    tfr_samples_to_remove=[];
                    for b=1:numel(blocks_not_present_in_triggers)
                        B=blocks_not_present_in_triggers(b);
                        start_to_remove=sum(n_LFP_samples_per_block(2,n_LFP_samples_per_block(1,:)<B))+1;
                        end_to_remove=n_LFP_samples_per_block(2,n_LFP_samples_per_block(1,:)==B)+start_to_remove-1;
                        b_idx=[site_LFP.block]==B;
                        FN_tr={'LFP_samples','dataset','block','run','n','LFP_tStart','LFP_t0'};
                        for f=1:numel(FN_tr)
                            FN=FN_tr{f};
                            site_LFP.(FN)(b_idx)=[];
                        end
                        tfr_samples_to_remove=[tfr_samples_to_remove start_to_remove:end_to_remove];
                    end
                    FN_tfs={'phabp','powbp','pha','pow','lfp'};
                    for f=1:numel(FN_tfs)
                        FN=FN_tfs{f};
                        site_LFP.tfs.(FN)(:,tfr_samples_to_remove)=[];
                    end
                    
                    blockstarts=[1, find(diff([trials.block]))+1];
                    t_offset_per_block=[[trials(blockstarts).block];[trials(blockstarts).TDT_LFPx_tStart]];
                    triggers.ecg.shuffled = ecg_bna_resample_triggers(Rpeaks,'shuffled_ts',t_offset_per_block,n_LFP_samples_per_block,1/site_LFP.tfs.sr);
                    triggers.ecg.real     = ecg_bna_resample_triggers(Rpeaks,'RPEAK_ts',t_offset_per_block,n_LFP_samples_per_block,1/site_LFP.tfs.sr);
                    
                    site_data=ecg_bna_compute_session_Rpeak_triggered_variables( site_LFP,triggers,trials,cfg );
                    triggered_session_data.sites(s) = site_data;
                    triggered_session_data.session = site_data.session;
                end
                
                % make a folder to save figures
                session_result_folder = fullfile(cfg.session_lfp_fldr);
                if ~exist(session_result_folder, 'dir')
                    mkdir(session_result_folder);
                end
                save(fullfile(session_result_folder, ['Triggered_session_' triggered_session_data.session '.mat']), 'triggered_session_data');
                clear session_proc_lfp site_data triggered_session_data;
                
            end
        end
    end
    
    %% average across sessions
    if cfg.process_population
        
        keys=ecg_bna_get_unit_list(cfg,0);
        cfg.site_IDS=keys.tuning_table(2:end,find_column_index(keys.tuning_table,'site_ID'));
        
        if cfg.process_LFP
            monkeys = unique({cfg.session_info.Monkey});
            cfg.monkey = [monkeys{:}];
            cfg.session_lfp_fldr = fullfile(cfg.analyse_lfp_folder, 'Per_Session');
            cfg.sites_lfp_fldr   = fullfile(cfg.analyse_lfp_folder, 'Per_Site');
            
            grand_avg = ecg_bna_compute_grand_avg(cfg,'w_units');
            grand_avg = ecg_bna_compute_grand_avg(cfg,'wo_units');
            grand_avg = ecg_bna_compute_grand_avg(cfg,'all');
        end
        
        if cfg.process_spikes
            SPK_PSTH=load_stuff(sessions_info,cfg,'SPK_root_results_fldr','','per_unit','Output');
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

