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
ecg_bna_location     =which('ecg_bna_define_settings');
github_folder        =ecg_bna_location(1:strfind(ecg_bna_location,['ECG-behavior-neural-analysis' filesep 'ecg_bna_define_settings'])-1);
for v = 1:length(versions)
    version = versions{v};
    ecg_bna_cfg = ecg_bna_define_settings(github_folder,project,version);
    
    %% read tuning table?
    keys=struct;
    keys=ph_general_settings(project,keys);
    project_specific_settings=[keys.db_folder 'ph_project_settings.m'];
    run(project_specific_settings)
    version_specific_settings=[keys.db_folder ecg_bna_cfg.spikes_version filesep 'ph_project_version_settings.m'];
    run(version_specific_settings)
    keys.anova_table_file=[keys.basepath_to_save ecg_bna_cfg.spikes_version filesep 'tuning_table_combined_CI.mat'];
    keys.tuning_table=ph_load_tuning_table(keys); %% load tuning table
    keys.tt.tasktypes= {'Fsac_opt';'Vsac_opt'};
%     keys.tt.SNR_rating=keys.cal.SNR_rating;
%     keys.tt.stability_rating=keys.cal.stablity; %:(
%     keys.tt.Single_rating=keys.cal.single_rating; %% :(
    keys.tt.FR=keys.cal.FR; %:(
    keys.tt.n_spikes=keys.cal.n_spikes; %% :(
    keys.monkey=''; %% empty because we ignore which monkey it is basically
    keys=ph_tuning_table_correction(keys);
    ecg_bna_cfg.unit_IDS=keys.tuning_table(2:end,1);
    
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
            Rpeaks=ecg_bna_compute_session_shuffled_Rpeaks(sessions_info(i),ecg_bna_cfg);
        end
        if ecg_bna_cfg.process_spikes && any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_spike_histogram'))            
            ecg_bna_compute_session_spike_histogram(sessions_info(i),Rpeaks,ecg_bna_cfg,trials); %% conditions are still a mess here
            ecg_bna_plot_session_spike_histogram(sessions_info(i),ecg_bna_cfg);
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
            % for  all data before june 2023 use " ecg_bna_process_LFP "
            % and after june "ecg_bna_process_LFP_newTask"
            
            load(sessions_info(i).Input_trials);
            
            ecg_bna_cfg.session_lfp_fldr = fullfile(ecg_bna_cfg.analyse_lfp_folder, 'Per_Session');
            ecg_bna_cfg.sites_lfp_fldr   = fullfile(ecg_bna_cfg.analyse_lfp_folder, 'Per_Site');
            
            %% this is new
            sitesdir=fileparts(sessions_info(i).Input_LFP{:});
            [sitefiles]=dir(sessions_info(i).Input_LFP{:});
            
            for s = 1:length(sitefiles)
                load([sitesdir filesep sitefiles(s).name], 'sites');
                site_LFP= ecg_bna_process_LFP(sites, ecg_bna_cfg, trials);
                if s==1
                    session_ecg = ecg_bna_combine_shuffled_Rpeaks(session_ecg, Rpeaks,ecg_bna_cfg,site_LFP.trials(1).tsample); %% adding rpeaks (and shuffled rpeaks) CAREFUL: this is with 2k sampling frequency
                end
                site_data=ecg_bna_compute_session_Rpeak_triggered_variables( site_LFP,session_ecg,ecg_bna_cfg.analyse_states, ecg_bna_cfg );
                session_data.sites(s) = site_data;
                session_data.session = site_data.session;
            end
            
            % make a folder to save figures
            session_result_folder = fullfile(ecg_bna_cfg.session_lfp_fldr);
            if ~exist(session_result_folder, 'dir')
                mkdir(session_result_folder);
            end
            save(fullfile(session_result_folder, ['Rpeak_triggered_session_' session_data.session '.mat']), 'session_data');        
            clear session_proc_lfp site_data session_data;
            
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
    %
        %% average across sessions
        if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_spike_histogram'))
            SPK_PSTH=load_spikes(sessions_info,'SPK_fldr','','','per_unit','Output');
            ecg_bna_avg_spike_histogram(SPK_PSTH,sessions_info, ecg_bna_cfg);
        end
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

function Out=load_spikes(sessions_info,subfield,subfolder,namepart,per,varname)
condition_labels = {'Rest', 'Task'};
for i = 1:length(sessions_info)
    disp(['Session ' num2str(i) ' out of ' num2str(length(sessions_info))])
    if isempty(subfolder) %% unfortunate inconsistent naming
        monkey=sessions_info(i).Monkey(1:3);
    else
        monkey=['_' sessions_info(i).Monkey(1:3)];
    end
    results_folder = fullfile(sessions_info(i).(subfield),per,subfolder);
    to_load=dir(fullfile(results_folder, [namepart monkey '_' sessions_info(i).Date '*.mat']));
    tmp = arrayfun(@(x) load([x.folder filesep x.name], varname), to_load, 'UniformOutput', false);
    
    % put common parameters in the structure
    Out{i}.unit_ID = cellfun(@(x) x.Output.unit_ID, tmp, 'UniformOutput', false);
    Out{i}.target = cellfun(@(x) x.Output.target, tmp, 'UniformOutput', false);
    Out{i}.quantSNR = cellfun(@(x) x.Output.quantSNR, tmp);
    Out{i}.Single_rating = cellfun(@(x) x.Output.Single_rating, tmp);
    Out{i}.stability_rating = cellfun(@(x) x.Output.stability_rating, tmp);
    
    field_names = fieldnames(tmp{1}.Output.Rest);
    for conditionNum = 1:2
        for fieldNum = 1:length(field_names)
            if strcmp(field_names{fieldNum}, 'Rts') || ...
                    strcmp(field_names{fieldNum}, 'Rds') || ...
                    strcmp(field_names{fieldNum}, 'Rds_perm') || ...
                    strcmp(field_names{fieldNum}, 'raster') || ...
                    ischar(tmp{1}.Output.(condition_labels{conditionNum}).(field_names{fieldNum}))
                Out{i}.(condition_labels{conditionNum}).(field_names{fieldNum}) = ...
                    cellfun(@(x) cat(1,x.Output.(condition_labels{conditionNum}).(field_names{fieldNum})), tmp, 'UniformOutput', false);
            else
                tmp_out = cellfun(@(x) x.Output.(condition_labels{conditionNum}).(field_names{fieldNum}), tmp, 'UniformOutput', false);
                if strcmp(field_names{fieldNum}, 'SDsubstractedSDP') || ...
                        strcmp(field_names{fieldNum}, 'SDsubstractedSDP_normalized')
                    len_tmp_out = cellfun(@length, tmp_out);
                    tmp_out(len_tmp_out == 1) = {single(nan(1, 51))};
                end
                tmp_out = cat(1, tmp_out{:});
                Out{i}.(condition_labels{conditionNum}).(field_names{fieldNum}) = tmp_out;
                clear tmp_out
            end
        end
    end
end
end
