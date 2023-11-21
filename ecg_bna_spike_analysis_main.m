function ecg_bna_spike_analysis_main(project, ecg_bna_cfg)
% This is a function that runs all the spike analysis related functions
% together.

%% setting section - hopefully temporary
perform_per_session_analysis=1; %% move this one to settings
only_plotting=0;

%% apply exclusion criteria, save lists of included and excluded units
keys=struct;
keys=ph_general_settings(project,keys);
project_specific_settings=[keys.db_folder 'ph_project_settings.m'];
run(project_specific_settings)
version_specific_settings=[keys.db_folder ecg_bna_cfg.spikes_version filesep 'ph_project_version_settings.m'];
run(version_specific_settings)
keys.anova_table_file=[keys.basepath_to_save ecg_bna_cfg.spikes_version filesep 'tuning_table_combined_CI.mat'];
keys.tuning_table=ph_load_tuning_table(keys); %% load tuning table

% create unit list before exclusion criteria and save it
unit_ids_before_exclusion = keys.tuning_table(2:end,1);
targets_before_exclusion  = keys.tuning_table(2:end,3);
filename = [ecg_bna_cfg.SPK_root_results_fldr filesep 'unitInfo_before_exclusion.mat'];
save(filename, 'unit_ids_before_exclusion', 'targets_before_exclusion')

% apply exclusion criteria - choose only cells for BOTH task and rest
keys.tt.tasktypes= {'Fsac_opt';'Vsac_opt'};
% keys.tt.SNR_rating=keys.cal.SNR_rating;
% keys.tt.stability_rating=keys.cal.stablity;
%     keys.tt.Single_rating=keys.cal.single_rating; %% :(
keys.tt.FR=keys.cal.FR; %:(
keys.tt.n_spikes=keys.cal.n_spikes; %% :(
keys.monkey=''; %% empty because we ignore which monkey it is basically
keys=ph_tuning_table_correction(keys);
unit_ids_after_exclusion = keys.tuning_table(2:end,1);
targets_after_exclusion  = keys.tuning_table(2:end,3);
filename = [ecg_bna_cfg.SPK_root_results_fldr filesep 'unitInfo_after_exclusion_stableTaskAndRest.mat'];
save(filename, 'unit_ids_after_exclusion', 'targets_after_exclusion')
clear unit_ids_after_exclusion targets_after_exclusion

% apply exclusion criteria - choose cells that had either a task or a rest
% block
% 1. choose units that had a task block
keys.tt.tasktypes = {'Fsac_opt'}; % fixation saccade - for rest
keys.tuning_table=ph_load_tuning_table(keys); %% load tuning table
keys = ph_tuning_table_correction(keys);
unit_ids_after_exclusion_rest = keys.tuning_table(2:end,1);
targets_after_exclusion_rest  = keys.tuning_table(2:end,3);

% 2. choose units that had a rest block
keys.tt.tasktypes= {'Vsac_opt'}; % fixation saccade - for rest
keys.tuning_table=ph_load_tuning_table(keys); %% load tuning table
keys = ph_tuning_table_correction(keys);
unit_ids_after_exclusion_task = keys.tuning_table(2:end,1);
targets_after_exclusion_task  = keys.tuning_table(2:end,3);

% 3. find units that had at least 1 rest block or 1 task block, figure out
% targets for those
unit_ids_after_exclusion = union(unit_ids_after_exclusion_rest, unit_ids_after_exclusion_task);
ids_rest = ismember(unit_ids_after_exclusion, unit_ids_after_exclusion_rest);
ids_task = ismember(unit_ids_after_exclusion, unit_ids_after_exclusion_task);
targets_after_exclusion = cell(length(unit_ids_after_exclusion), 1);
targets_after_exclusion(ids_rest) = targets_after_exclusion_rest;
targets_after_exclusion(ids_task) = targets_after_exclusion_task;

filename = [ecg_bna_cfg.SPK_root_results_fldr filesep 'unitInfo_after_exclusion.mat'];
save(filename, 'unit_ids_after_exclusion', 'targets_after_exclusion')
clear unit_ids_after_exclusion targets_after_exclusion
% ecg_bna_cfg.unit_IDS=keys.tuning_table(2:end,1);

%% Temporary - copy selected units separately
%     ecg_bna_copy_selected_units(ecg_bna_cfg)

%% Get info about sessions to be analysed
% Read the info about sessions to analyse
sessions_info = ecg_bna_cfg.session_info;

%% per session processing..
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
            ecg_bna_compute_session_spike_histogram(sessions_info(i),Rpeaks,ecg_bna_cfg);
            ecg_bna_plot_session_spike_histogram(sessions_info(i),ecg_bna_cfg);
        end
        
        if ecg_bna_cfg.process_spikes && any(strcmp(ecg_bna_cfg.analyses, 'spike_phase_ECG_cycle'))
            %             ecg_bna_compute_session_ECG_related_spikePhase(sessions_info(i),Rpeaks,ecg_bna_cfg)
        end
        
        if ecg_bna_cfg.process_spikes && any(strcmp(ecg_bna_cfg.analyses, 'spike_phase_ECG_cycle'))
            %             ecg_bna_plot_session_ECG_related_spikePhase(sessions_info(i),ecg_bna_cfg)
        end
    end
end

%% average across sessions
% if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_spike_histogram'))
%     SPK_PSTH=load_spikes(sessions_info,'SPK_fldr','','','per_unit','Output');
%     ecg_bna_avg_spike_histogram(SPK_PSTH,sessions_info, ecg_bna_cfg);
% end

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
                tmp_out = cat(1, tmp_out{:});
                Out{i}.(condition_labels{conditionNum}).(field_names{fieldNum}) = tmp_out;
                clear tmp_out
            end
        end
    end
end
end