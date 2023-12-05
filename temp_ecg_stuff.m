
            
            if cfg.process_ECG
                fprintf('Reading ECG for session %s\n', session_name);
                cfg.session = session_name;
                % read ECG data for each session
                % if ~isfield(sessions_info(i), 'Input_ECG') %% create ECG Rpeaks file!
                %       bsa_ecg_analyze_one_session(session_path,pathExcel,settings_filename,varargin);
                % end
                
                if isfield(sessions_info(i), 'Input_ECG_combined') && ~isempty(sessions_info(i).Input_ECG_combined)
                    session_ecg = ecg_bna_read_combined_ECG(sessions_info(i), cfg.plottrials); %this one isnt fixed at all (?)
                elseif isfield(sessions_info(i), 'Input_ECG_preproc') && ~isempty(sessions_info(i).Input_ECG_preproc)
                    session_ecg = ecg_bna_read_preproc_ECG(sessions_info(i));
                    %% Luba
                    %    session_ecg = ecg_bna_read_preproc_ECG_simple(sessions_info(i), 1);
                end
                session_ecg = ecg_bna_combine_shuffled_Rpeaks(session_ecg, Rpeaks,cfg); %% adding rpeaks (and shuffled rpeaks) CAREFUL: this is with 2k sampling frequency
                
                if isempty(fieldnames(session_ecg))
                    continue;
                end
                
                % folder to which results of analysis of this session should be stored
                cfg.session_ecg_fldr = fullfile(cfg.analyse_ecg_folder, session_name);
                
                % Calculate and plot the session average ECG, evoked response for different conditions
                
                %ecg_bna_compute_session_evoked_ECG( session_ecg,sessions_info(i), ecg_bna_cfg.analyse_states, ecg_bna_cfg );
                %             %% shuffle state onsets?
                %             ecg_bna_compute_session_state_evoked_ECG( session_ecg, sessions_info(i), ecg_bna_cfg.analyse_Rpeak_states, ecg_bna_cfg );
                %
                %             %% add shuffled Rpeaks !!
                %             ecg_bna_compute_session_Rpeak_evoked_state_onsets( session_ecg, sessions_info(i), ecg_bna_cfg.analyse_Rpeak_states, ecg_bna_cfg );
                %             %ecg_bna_compute_session_evoked_ECG_R2Rt( session_ecg, sessions_info(i), ecg_bna_cfg.event_triggers, ecg_bna_cfg );
                %
%             elseif ecg_bna_cfg.process_LFP
%                 session_ecg_filename = fullfile(sessions_info(i).proc_ecg_fldr, ['session_ecg_' sessions_info(i).session '.mat']);
%                 if exist(session_ecg_filename, 'file')
%                     load(session_ecg_filename, 'session_ecg');
%                     
%                 else
%                     seed=rng;
%                     save(seed_filename,'seed');
%                 end
            end