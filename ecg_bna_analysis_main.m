function ecg_bna_analysis_main(project,versions)
% project='Pulv_distractor_spatial_choice';
% versions={'ver_LS'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for LFP time frequency analysis
% Runs functions for reading processed LFP data, rejection of noise trials
% and task specific analysis using TFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 

% whether the LFP should be processed (true) or not (false)
% if the LFP for the sessions to analyse has already been processed, and no
% settings need to be changed, this flag can be set to false to skip LFP
% being processed again
% If LFP was not previously processed and the flag is set to false,
% analysis won't happen
% TODO: check if LFP is processed, of not, process LFP even if flag is set
% to false

%% INITIALIZATION
close all;

% loop through settings file
for v = 1:length(versions)
    version = versions{v};

    ecg_bna_cfg = ecg_bna_define_settings(project,version);
    
    %% Get info about sessions to be analysed

    % Read the info about sessions to analyse
    sessions_info = ecg_bna_cfg.session_info;

    %% initialize structs to store intermediate results
    % struct to store average ECG evoked response for different conditions
    Rpeak_evoked_ECG = struct();
    Rpeak_evoked_event_prob = struct();
    ECG_R2Rt_evoked = struct();
    % struct to store average LFP TFR for different conditions
    Rpeak_evoked_lfp_tfr = struct();
    % struct to store average LFP evoked response for different conditions
    Rpeak_evoked_lfp_raw = struct();
    
    % loop through each session to process
            SPK_PSTH = []; 
% 
%     for i = 1:length(sessions_info)
%         
%     
%         % name of session = [Monkey name '_' Recording date]
%         session_name = [sessions_info(i).Monkey '_' sessions_info(i).Date];
%         colormap gray;
%         cmp = colormap; 
%                
%         ecg_histogram(sessions_info(i),ecg_bna_cfg, cmp(i,:));
% 
%    end
        %% ECG spike histogram per session..
        for i = 1:length(sessions_info)

        if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_spike_histogram')) %% not quite sure yet where to put seed
            seed_filename=[ecg_bna_cfg.ECG_root_results_fldr filesep 'seed.mat'];
            if exist(seed_filename,'file');
                load(seed_filename);
                rng(seed);
            else
                seed=rng;
                save(seed_filename,'seed');
            end
            SPK_PSTH{i}=ecg_bna_compute_session_spike_histogram(sessions_info(i),ecg_bna_cfg);
        end
        
        
        if ecg_bna_cfg.process_ECG
            fprintf('Reading ECG for session %s\n', session_name);
            ecg_bna_cfg.session = session_name;
            % read ECG data for each session
            
%             if ~isfield(session_info, 'Input_ECG') %% create ECG Rpeaks file!
%                 bsa_ecg_analyze_one_session(session_path,pathExcel,settings_filename,varargin);
%             end
            
            if isfield(sessions_info(i), 'Input_ECG_combined') && ~isempty(sessions_info(i).Input_ECG_combined)
                session_ecg = ecg_bna_read_combined_ECG(sessions_info(i), ecg_bna_cfg.plottrials); %this one isnt fixed at all (?)
            elseif isfield(sessions_info(i), 'Input_ECG_preproc') && ~isempty(sessions_info(i).Input_ECG_preproc)
                session_ecg = ecg_bna_read_preproc_ECG(sessions_info(i), ecg_bna_cfg.plottrials); %this one is somewhat fixed
            end
            
            % Read LFP data
            if isfield(sessions_info(i), 'Input_LFP') && ~isempty(sessions_info(i).Input_LFP) && ecg_bna_cfg.process_LFP
                sessions_info(i) = ecg_bna_process_combined_LFP_ECG(sessions_info(i), ecg_bna_cfg); %this one should be fixed
            end
            
            if isempty(fieldnames(session_ecg))
                continue;
            end
        else
            disp('!!comment out - bug- KK')
            % load session ecg for the session
%             session_ecg_filename = fullfile(sessions_info(i).proc_ecg_fldr, ['session_ecg_' sessions_info(i).session_name '.mat']);
%             %session_ecg_filename = fullfile(sessions_info(i).proc_ecg_fldr, ['session_ecg_' sessions_info(i).Date '.mat']);
% 
%             if exist(session_ecg_filename, 'file')
%                 load(session_ecg_filename, 'session_ecg');
%             else
%                 continue;
%             end
        end
         if ecg_bna_cfg.process_LFP
        
        fprintf('Analysing for session %s\n', session_name);
        
        % folder to which results of analysis of this session should be
        % stored
        ecg_bna_cfg.session_ecg_fldr = fullfile(ecg_bna_cfg.analyse_ecg_folder, session_name);
        
        % Calculate and plot the session average ECG,
        % evoked response for different conditions
        
        Rpeak_evoked_ECG.session(i) = ecg_bna_compute_session_evoked_ECG( session_ecg, sessions_info(i), ecg_bna_cfg.analyse_states, ecg_bna_cfg );
        Rpeak_evoked_event_prob.session(i) = ecg_bna_compute_session_Rpeak_evoked_state_onsets( session_ecg, sessions_info(i), ecg_bna_cfg.analyse_Rpeak_states, ecg_bna_cfg );
        ECG_R2Rt_evoked.session(i) = ecg_bna_compute_session_evoked_ECG_R2Rt( session_ecg, sessions_info(i), ecg_bna_cfg.event_triggers, ecg_bna_cfg );
        
        %clear session_ecg;
        
        if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_LFP')) || any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_TFS'))
            ecg_bna_cfg.session_lfp_fldr = fullfile(ecg_bna_cfg.analyse_lfp_folder, 'Per_Session');
            ecg_bna_cfg.sites_lfp_fldr   = fullfile(ecg_bna_cfg.analyse_lfp_folder, 'Per_Site');
            % read the processed lfp mat files for all sites of this session
            sites_lfp_files = dir(fullfile(sessions_info(i).proc_lfp_fldr, ['site_lfp_pow_' sessions_info(i).Monkey(1:3) '_' sessions_info(i).Date '*.mat']));
            session_proc_lfp = [];
            for file = {sites_lfp_files.name}
                %fprintf('Reading processed LFP for site %s\n', file{1});
                if ~isempty(strfind(file{1}, 'site_lfp_'))
                    fprintf('Reading processed LFP for site %s\n', file{:});
                    load(fullfile(sessions_info(i).proc_lfp_fldr, file{1}))
                    session_proc_lfp = [session_proc_lfp site_lfp];
                end
            end
        end
        
        if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_LFP'))
            Rpeak_evoked_lfp_raw.session(i) = ecg_bna_compute_session_Rpeak_evoked_LFP( session_proc_lfp, ecg_bna_cfg.analyse_states, ecg_bna_cfg ); % this one is pretty fixed
        end
        if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_TFS'))
            Rpeak_evoked_lfp_tfr.session(i) = ecg_bna_compute_session_Rpeak_evoked_TFS( session_proc_lfp, ecg_bna_cfg.analyse_states, ecg_bna_cfg );
        end
        
        clear session_proc_lfp;
        
        %         tfs_ecg.session(i) = lfp_tfa_plot_session_tfs_ECG( session_ecg, ...
        %             sessions_info(i), lfp_tfa_cfg.event_triggers, lfp_tfa_cfg );
         end

         
         path_population=[sessions_info(1).SPK_fldr filesep 'Population' filesep 'SPK_PSTH'];
         if ~exist(path_population,'dir')
             mkdir(path_population);
         end
         save([sessions_info(1).SPK_fldr filesep 'Population' filesep 'SPK_PSTH'],'SPK_PSTH')
         
    end


    %% average across sessions
    if length(sessions_info) > 1
        % load all sessions
        if isempty( SPK_PSTH)
            ephys_folder=['Y:\Projects\' project '\ECG_triggered_spikes\' versions{1} filesep 'per_unit' filesep];
            Files =  dir([ephys_folder '*.mat']);
            fprintf(['Reading processed ECG-triggered spike data: ', Files.name]);

            for i = 1: length(Files)
                load([ephys_folder, Files(i).name]);
                SPK_PSTH{i} =Output; 
            end
        end
         
        if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_spike_histogram'))
            ecg_bna_avg_spike_histogram(SPK_PSTH,sessions_info);
        end
        
        % Average task evoked ECG
        Rpeak_evoked_ECG.sessions_avg = ecg_bna_avg_sessions_Rpeak_evoked_ECG(Rpeak_evoked_ECG, ecg_bna_cfg);
        % Average task evoked ECG b2bt
        

    % %     % Average task evoked ECG time frequency spectrogram
    %     tfs_ecg.sessions_avg = lfp_tfa_avg_sessions_ECG_tfs(tfs_ecg, lfp_tfa_cfg);
        % Average Rpeak evoked state onset probability
        
        %% these two functions here still have big issues!!
%         ECG_R2Rt_evoked.sessions_avg = ecg_bna_avg_sessions_ECGb2bt_evoked(ECG_R2Rt_evoked, ecg_bna_cfg);
%         Rpeak_evoked_event_prob.sessions_avg = ecg_bna_avg_sessions_Rpeak_evoked_state_onsets( Rpeak_evoked_event_prob, ecg_bna_cfg);
        
        if any(strcmp(ecg_bna_cfg.compute_avg_across, 'sessions'))
            if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_TFS'))
                ecg_bna_cfg.root_results_fldr=ecg_bna_cfg.results_folder;
                Rpeak_evoked_lfp_tfr.sessions_avg = lfp_tfa_avg_tfr_across_sessions(Rpeak_evoked_lfp_tfr.session, ecg_bna_cfg, fullfile(ecg_bna_cfg.root_results_fldr, 'ECG_triggered_avg_across_sessions'));
            end
            if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_LFP'))
                Rpeak_evoked_lfp_raw.sessions_avg = ecg_bna_avg_sessions_Rpeak_evoked_LFP(Rpeak_evoked_lfp_raw, ecg_bna_cfg);
            end
        end
        
        if any(strcmp(ecg_bna_cfg.compute_avg_across, 'sites'))
            if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_TFS'))
                ecg_bna_cfg.root_results_fldr=ecg_bna_cfg.results_folder;
                Rpeak_evoked_lfp_tfr.sessions_avg = lfp_tfa_avg_tfr_across_sites(Rpeak_evoked_lfp_tfr.session, ecg_bna_cfg, fullfile(ecg_bna_cfg.root_results_fldr, 'ECG_triggered_avg_across_sites'));
            end
            if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_LFP'))
                Rpeak_evoked_lfp_raw.sessions_avg = ecg_bna_avg_sites_Rpeak_evoked_LFP(Rpeak_evoked_lfp_raw, ecg_bna_cfg);
            end
        end
    end
end
end