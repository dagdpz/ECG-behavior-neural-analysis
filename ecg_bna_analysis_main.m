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

%Lukas was here

ecg_bna_location     =which('ecg_bna_define_settings');
github_folder        =ecg_bna_location(1:strfind(ecg_bna_location,['ECG-behavior-neural-analysis' filesep 'ecg_bna_define_settings'])-1);
for v = 1:length(versions)
    version = versions{v};
    ecg_bna_cfg = ecg_bna_define_settings(github_folder,project,version);
    
    %% do spike analysis if necessary
    if ecg_bna_cfg.process_spikes && any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_spike_histogram'))
        ecg_bna_spike_analysis_main(project, ecg_bna_cfg)
    end
    
    %% Get info about sessions to be analysed
    % Read the info about sessions to analyse
    sessions_info = ecg_bna_cfg.session_info;
    
    %% per session processing..
    perform_per_session_analysis=1; %% move this one to settings
    only_plotting=0;
    if perform_per_session_analysis
    for i = 1:length(sessions_info)
        java.lang.System.gc() % added by Luba to control the number of graphical device interface handles (hope to handle the problem with freezing plots while creating figures)
            
        session_name = [sessions_info(i).Monkey '_' sessions_info(i).Date];
        load(sessions_info(i).Input_trials);
        
        % First make seed and ecg shuffles, then use those shuffles for all subfunctions
        seed_filename=[ecg_bna_cfg.ECG_root_results_fldr filesep 'seed.mat']; %% not quite sure yet where to put seed
        if exist(seed_filename,'file')
            load(seed_filename);
            rng(seed);
        else
            seed=rng;
            save(seed_filename,'seed');
        end
        if ecg_bna_cfg.process_LFP || ecg_bna_cfg.process_ECG || ecg_bna_cfg.process_spikes
            %Rpeaks=ecg_bna_jitter(sessions_info(i),ecg_bna_cfg);
            Rpeaks=ecg_bna_compute_session_shuffled_Rpeaks(sessions_info(i),ecg_bna_cfg);
        end
            
        if ecg_bna_cfg.process_ECG
            fprintf('Reading ECG for session %s\n', session_name);
            ecg_bna_cfg.session = session_name;
            % read ECG data for each session
            % if ~isfield(sessions_info(i), 'Input_ECG') %% create ECG Rpeaks file!
            %       bsa_ecg_analyze_one_session(session_path,pathExcel,settings_filename,varargin);
            % end
            
            if isfield(sessions_info(i), 'Input_ECG_combined') && ~isempty(sessions_info(i).Input_ECG_combined)
                session_ecg = ecg_bna_read_combined_ECG(sessions_info(i), ecg_bna_cfg.plottrials); %this one isnt fixed at all (?)
            elseif isfield(sessions_info(i), 'Input_ECG_preproc') && ~isempty(sessions_info(i).Input_ECG_preproc)
                session_ecg = ecg_bna_read_preproc_ECG(sessions_info(i));
                %% Luba
                %    session_ecg = ecg_bna_read_preproc_ECG_simple(sessions_info(i), 1);
            end
            session_ecg = ecg_bna_combine_shuffled_Rpeaks(session_ecg, Rpeaks,ecg_bna_cfg); %% adding rpeaks (and shuffled rpeaks) CAREFUL: this is with 2k sampling frequency
            
            if isempty(fieldnames(session_ecg))
                continue;
            end
            
            % folder to which results of analysis of this session should be stored
            ecg_bna_cfg.session_ecg_fldr = fullfile(ecg_bna_cfg.analyse_ecg_folder, session_name);
            
            % Calculate and plot the session average ECG, evoked response for different conditions
            
            %ecg_bna_compute_session_evoked_ECG( session_ecg,sessions_info(i), ecg_bna_cfg.analyse_states, ecg_bna_cfg );
            %             %% shuffle state onsets?
            %             ecg_bna_compute_session_state_evoked_ECG( session_ecg, sessions_info(i), ecg_bna_cfg.analyse_Rpeak_states, ecg_bna_cfg );
            %
            %             %% add shuffled Rpeaks !!
            %             ecg_bna_compute_session_Rpeak_evoked_state_onsets( session_ecg, sessions_info(i), ecg_bna_cfg.analyse_Rpeak_states, ecg_bna_cfg );
            %             %ecg_bna_compute_session_evoked_ECG_R2Rt( session_ecg, sessions_info(i), ecg_bna_cfg.event_triggers, ecg_bna_cfg );
            %
        elseif ecg_bna_cfg.process_LFP
            session_ecg_filename = fullfile(sessions_info(i).proc_ecg_fldr, ['session_ecg_' sessions_info(i).session '.mat']);
            if exist(session_ecg_filename, 'file')
                load(session_ecg_filename, 'session_ecg');
                
            else
                seed=rng;
                save(seed_filename,'seed');
            end
        end
        if ecg_bna_cfg.process_LFP
            
            fprintf('Analysing for session %s\n', session_name);
            
            % Read LFP data
            sessions_info(i).proc_results_fldr=sessions_info(i).proc_lfp_fldr;
            ecg_bna_cfg.session_lfp_fldr = fullfile(ecg_bna_cfg.analyse_lfp_folder, 'Per_Session');
            ecg_bna_cfg.sites_lfp_fldr   = fullfile(ecg_bna_cfg.analyse_lfp_folder, 'Per_Site');
            
            %% this is new
            sitesdir=fileparts(sessions_info(i).Input_LFP{:});
            [sitefiles]=dir(sessions_info(i).Input_LFP{:});
            sr=unique([trials.TDT_LFPx_SR]);
            ts_original=1/sr;
                    
            for s = 1:length(sitefiles)
                load([sitesdir filesep sitefiles(s).name], 'sites');                
                site_LFP= ecg_bna_process_LFP(sites, ecg_bna_cfg, ts_original);
                
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
                
                site_data=ecg_bna_compute_session_Rpeak_triggered_variables( site_LFP,triggers,trials,ecg_bna_cfg );
                triggered_session_data.sites(s) = site_data;
                triggered_session_data.session = site_data.session;
            end
            
            % make a folder to save figures
            session_result_folder = fullfile(ecg_bna_cfg.session_lfp_fldr);
            if ~exist(session_result_folder, 'dir')
                mkdir(session_result_folder);
            end
            save(fullfile(session_result_folder, ['Triggered_session_' triggered_session_data.session '.mat']), 'triggered_session_data');        
            clear session_proc_lfp site_data triggered_session_data;
            
        end
    end
    end
    
    %% per session? (why not per site??)
    if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_LFP'))
        %% eventually, we want 3 things (all of this is target wise):
        session_Rpeak_triggered_raw=load_stuff(sessions_info,'analyse_lfp_fldr','Rpeak_triggered_session_','Per_Session','session_data');
            % a) average across sites for each session (first cosntruct this)
            %% the same function should be able to either 
            % b) average of those session averages (append all session averages, and average)
                     % as a part of this, add scatter plots
            
            % c) GRAND average across all sites <-- (append all sites, and average)
                     % as a part of this, add scatter plots
        
%         %% loop through sessions
%         for s=1:numel(sessions_info)
%         session_Rpeak_triggered_raw(s)=load_stuff(sessions_info(s),'analyse_lfp_fldr','Rpeak_triggered_session_','Per_Session','session_data');
%         session_Rpeak_triggered_raw(s).sites_avg = ecg_bna_avg_site_and_sessions_Rpeak_triggered_results(session_Rpeak_triggered_raw, ecg_bna_cfg, 'sites');
% 
%         end
        
        if any(strcmp(ecg_bna_cfg.compute_avg_across, 'sessions'))
            session_Rpeak_triggered_raw.sessions_avg = ecg_bna_avg_site_and_sessions_Rpeak_triggered_results(session_Rpeak_triggered_raw, ecg_bna_cfg,'sessions');
        end
        if any(strcmp(ecg_bna_cfg.compute_avg_across, 'sites'))
            session_Rpeak_triggered_raw.sites_avg = ecg_bna_avg_site_and_sessions_Rpeak_triggered_results(session_Rpeak_triggered_raw, ecg_bna_cfg, 'sites');
        end
        session_data = session_Rpeak_triggered_raw;
        save(fullfile(ecg_bna_cfg.session_lfp_fldr, ['Rpeak_triggered_session_' session_Rpeak_triggered_raw.session '.mat']), 'session_data');
        clear session_data
        clear session_Rpeak_triggered_raw
    end
    
    %%
    % computing the grand avg over all sites of all sessions in each Target
    
    
    ecg_bna_cfg.session_lfp_fldr = fullfile(ecg_bna_cfg.analyse_lfp_folder, 'Per_Session');
    ecg_bna_cfg.sites_lfp_fldr   = fullfile(ecg_bna_cfg.analyse_lfp_folder, 'Per_Site');
    
%     grand_avg = ecg_bna_compute_grand_avg(ecg_bna_cfg,'w_units');
%     grand_avg = ecg_bna_compute_grand_avg(ecg_bna_cfg,'wo_units');
    grand_avg = ecg_bna_compute_grand_avg(ecg_bna_cfg,'all');
    %%
    %     if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_TFS'))
    %         ecg_bna_cfg.root_results_fldr=ecg_bna_cfg.results_folder; %??
    %         Rpeak_evoked_lfp_tfr=load_stuff(sessions_info,'analyse_lfp_fldr','Rpeak_evoked_TFS','LFP_TFR','Per_Session','session_tfs');
    %         if any(strcmp(ecg_bna_cfg.compute_avg_across, 'sessions'))
    %             Rpeak_evoked_lfp_tfr.sessions_avg = lfp_tfa_avg_tfr_across_sessions(Rpeak_evoked_lfp_tfr.session, ecg_bna_cfg);
    %         end
    %         if any(strcmp(ecg_bna_cfg.compute_avg_across, 'sites'))
    %             Rpeak_evoked_lfp_tfr.sites_avg = lfp_tfa_avg_tfr_across_sites(Rpeak_evoked_lfp_tfr.session, ecg_bna_cfg);
    %         end
    %         clear Rpeak_evoked_lfp_tfr
    %     end
    %     if false
    %
    %         % Average task evoked ECG b2bt - does not exist yet :D
    %
    %         % %     % Average task evoked ECG time frequency spectrogram
    %         %     tfs_ecg.sessions_avg = lfp_tfa_avg_sessions_ECG_tfs(tfs_ecg, lfp_tfa_cfg);
    %         %       Average Rpeak evoked state onset probability
    %
    %         %% these functions here are still not clear!!
    %         % Average task evoked ECG
    %         Rpeak_evoked_ECG                        =load_stuff(sessions_info,'analyse_ecg_fldr','','ECG_evoked_','Rpeak_evoked_ECG','site_evoked_ecg');
    %         Rpeak_evoked_ECG.sessions_avg           = ecg_bna_avg_sessions_Rpeak_evoked_ECG(Rpeak_evoked_ECG, ecg_bna_cfg);
    %         ECG_R2Rt_evoked                         =load_stuff(sessions_info,'analyse_ecg_fldr','','ECG_R2Rt_evoked_','ECG R2Rt','session_evoked_ecg_R2Rt');
    %         ECG_R2Rt_evoked.sessions_avg            = ecg_bna_avg_sessions_ECGb2bt_evoked(ECG_R2Rt_evoked, ecg_bna_cfg);
    %         Rpeak_evoked_event_prob                 =load_stuff(sessions_info,'analyse_ecg_fldr','','Rpeak_Evoked_state_onsets_','Rpeak_Evoked_states','site_evoked_ecg');
    %         Rpeak_evoked_event_prob.sessions_avg    = ecg_bna_avg_sessions_Rpeak_evoked_state_onsets( Rpeak_evoked_event_prob, ecg_bna_cfg);
    %     end
    %
end
end 

function Out = load_stuff(sessions_info,subfolder,namepart,per,varname)
for i = 1:length(sessions_info)
    if ~isempty(subfolder) %% unfortunate inconsistent naming
        monkey=sessions_info(i).Monkey(1:3);
        
        results_folder = fullfile(sessions_info(i).(subfolder),per);
        to_load = load(fullfile(results_folder, [namepart monkey '_' sessions_info(i).Date '.mat']),varname);
        Out=to_load.(varname);
    else
        error('The Per session folder is empty');
    end
end
end
