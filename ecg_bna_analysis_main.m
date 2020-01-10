%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for LFP time frequency analysis
% Runs functions for reading processed LFP data, rejection of noise trials
% and task specific analysis using TFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; 

% file containing settings for LFP analysis

settings_filepaths = ...%{'C:\Users\snair\Documents\GitHub\ECG-behavior-neural-analysis\settings\Magnus\ecg_bna_settings_Magnus_dPul_ECG_LFP.m'};
   {'C:\Users\snair\Documents\GitHub\ECG-behavior-neural-analysis\settings\Curius\ecg_bna_settings_Curius_inactivation.m'};%, ...
%     'C:\Users\snair\Documents\GitHub\ECG-behavior-neural-analysis\settings\Curius\ecg_bna_settings_Curius_baseline1.m', ...
%     'C:\Users\snair\Documents\GitHub\ECG-behavior-neural-analysis\settings\Cornelius\ecg_bna_settings_Cornelius_inactivation_ecg_by_block.m', ...
%     'C:\Users\snair\Documents\GitHub\ECG-behavior-neural-analysis\settings\Cornelius\ecg_bna_settings_Cornelius_baseline_ecg_by_block.m'};


% whether the LFP should be processed (true) or not (false)
% if the LFP for the sessions to analyse has already been processed, and no
% settings need to be changed, this flag can be set to false to skip LFP
% being processed again
% If LFP was not previously processed and the flag is set to false,
% analysis won't happen
% TODO: check if LFP is processed, of not, process LFP even if flag is set
% to false
process_ECG = 1; process_LFP = 0;

%% INITIALIZATION
close all;

% loop through settings file
for s = 1:length(settings_filepaths)
    settings_filepath = settings_filepaths{s};

    ecg_bna_cfg = ecg_bna_define_settings(settings_filepath);

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

    %try
        % loop through each session to process
        for i = 1:length(sessions_info)
            % name of session = [Monkey name '_' Recording date]
            session_name = [sessions_info(i).Monkey '_' sessions_info(i).Date];

            if process_ECG
                fprintf('Reading ECG for session %s\n', session_name);
                ecg_bna_cfg.session = session_name;
                % read ECG data for each session
                if isfield(sessions_info(i), 'Input_ECG_combined') && ...
                        ~isempty(sessions_info(i).Input_ECG_combined)
                    session_ecg = ...
                        ecg_bna_read_combined_ECG(sessions_info(i), ecg_bna_cfg.plottrials);
                elseif isfield(sessions_info(i), 'Input_ECG_preproc') && ...
                        ~isempty(sessions_info(i).Input_ECG_preproc)
                    session_ecg = ...
                        ecg_bna_read_preproc_ECG(sessions_info(i), ecg_bna_cfg.plottrials);
                end
                
                % Read LFP data
                if isfield(sessions_info(i), 'Input_LFP') && ...
                        ~isempty(sessions_info(i).Input_LFP) && ecg_bna_cfg.process_LFP
                    sessions_info(i) = ...
                        ecg_bna_process_combined_LFP_ECG(sessions_info(i), ecg_bna_cfg);
                end

                if isempty(fieldnames(session_ecg))
                    continue;
                end
            else
                % load session ecg for the session
                session_ecg_filename = fullfile(sessions_info(i).proc_ecg_fldr, ...
                    ['session_ecg_' session_name '.mat']);
                if exist(session_ecg_filename, 'file')
                    load(session_ecg_filename, 'session_ecg');
                else
                    continue;
                end
            end

            fprintf('Analysing for session %s\n', session_name);

            % folder to which results of analysis of this session should be
            % stored
            ecg_bna_cfg.session_ecg_fldr = ...
                fullfile(ecg_bna_cfg.analyse_ecg_folder, session_name);  
        
            % Calculate and plot the session average ECG, 
            % evoked response for different conditions 
            
            Rpeak_evoked_ECG.session(i) = ecg_bna_compute_session_evoked_ECG( session_ecg, ...
                sessions_info(i), ecg_bna_cfg.analyse_states, ecg_bna_cfg );
            
            Rpeak_evoked_event_prob.session(i) = ecg_bna_compute_session_Rpeak_evoked_state_onsets...
                ( session_ecg, sessions_info(i), ecg_bna_cfg.analyse_Rpeak_states, ...
                ecg_bna_cfg );
            
            ECG_R2Rt_evoked.session(i) = ecg_bna_compute_session_evoked_ECG_R2Rt( session_ecg, ...
                    sessions_info(i), ecg_bna_cfg.event_triggers, ecg_bna_cfg );
            
            %clear session_ecg;
            
            if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_LFP')) || ...
                    any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_TFS'))
                ecg_bna_cfg.session_lfp_fldr = ...
                    fullfile(ecg_bna_cfg.analyse_lfp_folder, session_name);
                % read the processed lfp mat files for all sites of this session
                sites_lfp_files = dir(fullfile(sessions_info(i).proc_lfp_fldr, 'site_lfp_*.mat'));
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
                Rpeak_evoked_lfp_raw.session(i) = ...
                    ecg_bna_compute_session_Rpeak_evoked_LFP( session_proc_lfp, ...
                    ecg_bna_cfg.analyse_states, ecg_bna_cfg );
            end
            if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_TFS'))
                Rpeak_evoked_lfp_tfr.session(i) = ...
                    ecg_bna_compute_session_Rpeak_evoked_TFS( session_proc_lfp, ...
                    ecg_bna_cfg.analyse_states, ecg_bna_cfg );
            end

            clear session_proc_lfp;
            
    %         tfs_ecg.session(i) = lfp_tfa_plot_session_tfs_ECG( session_ecg, ...
    %             sessions_info(i), lfp_tfa_cfg.event_triggers, lfp_tfa_cfg );       

        end
    %catch e
    %    error(e.message());
    %end

    %% average across sessions
    if length(sessions_info) > 1
        % Average task evoked ECG
        Rpeak_evoked_ECG.sessions_avg = ecg_bna_avg_sessions_Rpeak_evoked_ECG(Rpeak_evoked_ECG, ...
            ecg_bna_cfg);
        % Average task evoked ECG b2bt
        ECG_R2Rt_evoked.sessions_avg = ecg_bna_avg_sessions_ECGb2bt_evoked(ECG_R2Rt_evoked, ...
            ecg_bna_cfg);
    % %     % Average task evoked ECG time frequency spectrogram
    %     tfs_ecg.sessions_avg = lfp_tfa_avg_sessions_ECG_tfs(tfs_ecg, ...
    %         lfp_tfa_cfg);
        % Average Rpeak evoked state onset probability
        Rpeak_evoked_event_prob.sessions_avg = ecg_bna_avg_sessions_Rpeak_evoked_state_onsets( ...
            Rpeak_evoked_event_prob, ecg_bna_cfg);
        
        if any(strcmp(ecg_bna_cfg.compute_avg_across, 'sessions'))
            if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_TFS'))
                Rpeak_evoked_lfp_tfr.sessions_avg = ...
                    lfp_tfa_avg_tfr_across_sessions(Rpeak_evoked_lfp_tfr, ecg_bna_cfg);
            end
            if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_LFP'))
                Rpeak_evoked_lfp_raw.sessions_avg = ...
                    ecg_bna_avg_sessions_Rpeak_evoked_LFP(Rpeak_evoked_lfp_raw, ecg_bna_cfg);
            end
        end
        
        if any(strcmp(ecg_bna_cfg.compute_avg_across, 'sites'))
            if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_TFS'))
                Rpeak_evoked_lfp_tfr.sessions_avg = ...
                    lfp_tfa_avg_tfr_across_sites(Rpeak_evoked_lfp_tfr, ecg_bna_cfg);
            end
            if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_LFP'))
                Rpeak_evoked_lfp_raw.sessions_avg = ...
                    ecg_bna_avg_sites_Rpeak_evoked_LFP(Rpeak_evoked_lfp_raw, ecg_bna_cfg);
            end
        end
        
    end
    
end