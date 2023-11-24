function ecg_bna_spike_analysis_main(project, ecg_bna_cfg)
% This is a function that runs all the spike analysis related functions
% together.

%% setting section - hopefully temporary
compute_unit_subsets      = 0;
move_files                = 0;

compute_spike_histograms  = 0;
plot_spike_histograms     = 0;
compute_spike_phase       = 0;
plot_spike_phase          = 0;

population_analysis       = 1;

%% apply exclusion criteria, save lists of included and excluded units
if compute_unit_subsets
    keys=struct;
    keys=ph_general_settings(project,keys);
    project_specific_settings=[keys.db_folder 'ph_project_settings.m'];
    run(project_specific_settings)
    version_specific_settings=[keys.db_folder ecg_bna_cfg.spikes_version filesep 'ph_project_version_settings.m'];
    run(version_specific_settings)
    keys.anova_table_file=[keys.basepath_to_save ecg_bna_cfg.spikes_version filesep 'tuning_table_combined_CI.mat'];
    keys.tuning_table = ph_load_tuning_table(keys); %% load tuning table
    keys.tuning_table(2:end,3) = cellfun(@(x) x(1:end-2), keys.tuning_table(2:end,3), 'UniformOutput', false);
    
    % I. Unit list before exclusion - this one doesn't take number of
    % R-peaks into account
    % create unit list before exclusion criteria and save it
    unit_ids_before_exclusion = keys.tuning_table(2:end,1);
    targets_before_exclusion  = keys.tuning_table(2:end,3);
    filename = [ecg_bna_cfg.SPK_root_results_fldr filesep 'unitInfo_before_exclusion'];
    save(filename, 'unit_ids_before_exclusion', 'targets_before_exclusion')
    % create table of unit numbers by area
    [unique_areas, ~, ic] = unique(targets_before_exclusion);
    unit_counts = accumarray(ic,1);
    T = table(unique_areas, unit_counts);
    writetable(T, [filename '.xls'])
%     clear T unique_areas ic unit_counts filename
    
    % II. Unit list for BOTH task and rest
    % 1. apply exclusion criteria - choose only cells for BOTH task and rest
    keys.tt.tasktypes= {'Fsac_opt';'Vsac_opt'};
    keys.tt.FR=keys.cal.FR;
    keys.tt.n_spikes=keys.cal.n_spikes;
    keys.monkey=''; %% empty because we ignore which monkey it is basically
    keys=ph_tuning_table_correction(keys);
    unit_ids_after_exclusion = keys.tuning_table(2:end,1);
    targets_after_exclusion  = keys.tuning_table(2:end,3);
    filename = [ecg_bna_cfg.SPK_root_results_fldr filesep 'unitInfo_after_exclusion_stableTaskAndRest'];
    save(filename, 'unit_ids_after_exclusion', 'targets_after_exclusion')
    % 2. create table of unit numbers by area
    [unique_areas, ~, ic] = unique(targets_after_exclusion);
    unit_counts = accumarray(ic,1);
    T = table(unique_areas, unit_counts);
    writetable(T, [filename '.xls'])
    % 3. exclude additionally by the number of R-peaks
    [num_Rpeaks_rest, num_Rpeaks_task] = count_Rpeaks_per_unit(unit_ids_after_exclusion);
    ids_enough_Rpeaks = num_Rpeaks_rest > ecg_bna_cfg.unit_exclusion.nCardiacCycles & ...
        num_Rpeaks_task > ecg_bna_cfg.unit_exclusion.nCardiacCycles;
    unit_ids_after_exclusion_Rpeaks = unit_ids_after_exclusion(ids_enough_Rpeaks);
    targets_after_exclusion_Rpeaks  = targets_after_exclusion(ids_enough_Rpeaks);
    save([filename '_600'], 'unit_ids_after_exclusion_Rpeaks', 'targets_after_exclusion_Rpeaks')
    % 4. save the table with unit counts
    [unique_areas, ~, ic] = unique(targets_after_exclusion_Rpeaks);
    unit_counts = accumarray(ic,1);
    T = table(unique_areas, unit_counts);
    writetable(T, [filename '_600.xls'])
%     clear T unique_areas ic unit_counts filename unit_ids_after_exclusion targets_after_exclusion
    
    % III. Unit list of excluded units (that aren't units existing for both task and rest)
    ids_excluded = ~ismember(unit_ids_before_exclusion, unit_ids_after_exclusion);
    unit_ids_excluded = unit_ids_before_exclusion(ids_excluded);
    targets_excluded  = targets_before_exclusion(ids_excluded);
    filename = [ecg_bna_cfg.SPK_root_results_fldr filesep 'unitInfo_excluded_stableTaskAndRest'];
    save(filename, 'unit_ids_excluded', 'targets_excluded')
    [unique_areas, ~, ic] = unique(targets_excluded);
    unit_counts = accumarray(ic,1);
    T = table(unique_areas, unit_counts);
    writetable(T, [filename '.xls'])
%     clear T unique_areas ic unit_counts filename unit_ids_excluded targets_excluded
    
    % IV. Unit list for either task or rest
    % apply exclusion criteria - choose cells that had either a task or a rest
    % block
    % 1. choose units that had a task block
    keys.tt.tasktypes = {'Fsac_opt'}; % fixation saccade - for rest
    keys.tuning_table=ph_load_tuning_table(keys); %% load tuning table
    keys.tuning_table(2:end,3) = cellfun(@(x) x(1:end-2), keys.tuning_table(2:end,3), 'UniformOutput', false);
    keys = ph_tuning_table_correction(keys);
    unit_ids_after_exclusion_rest = keys.tuning_table(2:end,1);
    targets_after_exclusion_rest  = keys.tuning_table(2:end,3);
    
    % 1.1. figure out their R-peak counts
    [num_Rpeaks_rest, ~] = count_Rpeaks_per_unit(unit_ids_after_exclusion_rest);
    ids_enough_Rpeaks = num_Rpeaks_rest > ecg_bna_cfg.unit_exclusion.nCardiacCycles;
    unit_ids_after_exclusion_rest_Rpeaks = unit_ids_after_exclusion_rest(ids_enough_Rpeaks);
    targets_after_exclusion_rest_Rpeaks  = targets_after_exclusion_rest(ids_enough_Rpeaks);
    
    % 2. choose units that had a rest block
    keys.tt.tasktypes= {'Vsac_opt'}; % fixation saccade - for rest
    keys.tuning_table=ph_load_tuning_table(keys); %% load tuning table
    keys.tuning_table(2:end,3) = cellfun(@(x) x(1:end-2), keys.tuning_table(2:end,3), 'UniformOutput', false);
    keys = ph_tuning_table_correction(keys);
    unit_ids_after_exclusion_task = keys.tuning_table(2:end,1);
    targets_after_exclusion_task  = keys.tuning_table(2:end,3);
    
    % 2.1. figure out their R-peak counts
    [~, num_Rpeaks_task] = count_Rpeaks_per_unit(unit_ids_after_exclusion_task);
    ids_enough_Rpeaks = num_Rpeaks_task > ecg_bna_cfg.unit_exclusion.nCardiacCycles;
    unit_ids_after_exclusion_task_Rpeaks = unit_ids_after_exclusion_task(ids_enough_Rpeaks);
    targets_after_exclusion_task_Rpeaks  = targets_after_exclusion_task(ids_enough_Rpeaks);
    
    % 3. find units that had at least 1 rest block or 1 task block, figure out
    % targets for those
    unit_ids_after_exclusion = union(unit_ids_after_exclusion_rest, unit_ids_after_exclusion_task);
    ids_rest = ismember(unit_ids_after_exclusion, unit_ids_after_exclusion_rest);
    ids_task = ismember(unit_ids_after_exclusion, unit_ids_after_exclusion_task);
    targets_after_exclusion = cell(length(unit_ids_after_exclusion), 1);
    targets_after_exclusion(ids_rest) = targets_after_exclusion_rest;
    targets_after_exclusion(ids_task) = targets_after_exclusion_task;
    
    filename = [ecg_bna_cfg.SPK_root_results_fldr filesep 'unitInfo_after_exclusion'];
    save(filename, 'unit_ids_after_exclusion', 'targets_after_exclusion')
    % clear unit_ids_after_exclusion targets_after_exclusion
    [unique_areas, ~, ic] = unique(targets_after_exclusion);
    unit_counts = accumarray(ic,1);
    T = table(unique_areas, unit_counts);
    writetable(T, [filename '.xls'])
%     clear T unique_areas ic unit_counts unit_ids_excluded targets_after_exclusion
    
    % 3.1. find units that had at least 1 rest block or 1 task block AND
    % enough R-peaks for subsequent analysis
    unit_ids_after_exclusion_Rpeaks = ...
        union(unit_ids_after_exclusion_rest_Rpeaks, unit_ids_after_exclusion_task_Rpeaks);
    ids_rest_Rpeaks = ismember(unit_ids_after_exclusion_Rpeaks, unit_ids_after_exclusion_rest_Rpeaks);
    ids_task_Rpeaks = ismember(unit_ids_after_exclusion_Rpeaks, unit_ids_after_exclusion_task_Rpeaks);
    targets_after_exclusion_Rpeaks = cell(length(unit_ids_after_exclusion_Rpeaks), 1);
    targets_after_exclusion(ids_rest_Rpeaks) = targets_after_exclusion_rest_Rpeaks;
    targets_after_exclusion(ids_task_Rpeaks) = targets_after_exclusion_task_Rpeaks;
    save([filename '_600'], 'unit_ids_after_exclusion_Rpeaks', 'targets_after_exclusion_Rpeaks')
    [unique_areas, ~, ic] = unique(targets_after_exclusion);
    unit_counts = accumarray(ic,1);
    if size(unique_areas,1) < size(unique_areas,2)
        unique_areas = unique_areas';
    end
    T = table(unique_areas, unit_counts);
    writetable(T, [filename '_600.xls'])
    
    % V. Unit list of excluded units (that aren't units existing either task or rest)
    ids_excluded = ~ismember(unit_ids_before_exclusion, unit_ids_after_exclusion);
    unit_ids_excluded = unit_ids_before_exclusion(ids_excluded);
    targets_excluded  = targets_before_exclusion(ids_excluded);
    filename = [ecg_bna_cfg.SPK_root_results_fldr filesep 'unitInfo_excluded'];
    save(filename, 'unit_ids_excluded', 'targets_excluded')
    [unique_areas, ~, ic] = unique(targets_excluded);
    unit_counts = accumarray(ic,1);
    T = table(unique_areas, unit_counts);
    writetable(T, [filename '.xls'])
%     clear T unique_areas ic unit_counts filename unit_ids_excluded targets_excluded unit_ids_excluded targets_excluded unit_ids_before_exclusion targets_before_exclusion
end

if move_files
    %% Temporary - copy selected units separately
    ecg_bna_copy_selected_units(ecg_bna_cfg)
end
% ecg_bna_cfg.unit_IDS=keys.tuning_table(2:end,1);

%% Get info about sessions to be analysed
% Read the info about sessions to analyse
sessions_info = ecg_bna_cfg.session_info;

%% per session processing..
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
    
    if compute_spike_histograms
        Rpeaks=ecg_bna_compute_session_shuffled_Rpeaks(sessions_info(i),ecg_bna_cfg);
        ecg_bna_compute_session_spike_histogram(sessions_info(i),Rpeaks,ecg_bna_cfg);
    end
    
    if plot_spike_histograms
        ecg_bna_plot_session_spike_histogram(sessions_info(i),ecg_bna_cfg);
    end
    
    if compute_spike_phase
        Rpeaks=ecg_bna_compute_session_shuffled_Rpeaks(sessions_info(i),ecg_bna_cfg);
        ecg_bna_compute_session_ECG_related_spikePhase(sessions_info(i),Rpeaks,ecg_bna_cfg)
    end
    
    if plot_spike_phase
        ecg_bna_plot_session_ECG_related_spikePhase(sessions_info(i),ecg_bna_cfg)
    end
end

%% average across sessions
if population_analysis
    SPK_PSTH=load_spikes(sessions_info,'SPK_fldr','','','per_unit','Output');
    ecg_bna_avg_spike_histogram(SPK_PSTH,sessions_info, ecg_bna_cfg);
end

end

function [num_Rpeaks_rest, num_Rpeaks_task] = count_Rpeaks_per_unit(unit_list)

num_Rpeaks_rest = nan(length(unit_list),1);
num_Rpeaks_task = nan(length(unit_list),1);

for fileNum = 1:length(unit_list)
    
    currFileInfo = dir(['Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_TaskRest_Magnus_merged\per_unit\' unit_list{fileNum} '_*.mat']);
    
    if ~isempty(currFileInfo)
        
        disp(['Reading file ' num2str(fileNum) ' out of ' num2str(length(unit_list))])
        load([currFileInfo.folder filesep currFileInfo.name], 'Output')
        
        num_Rpeaks_rest(fileNum) = Output.Rest.NrEvents;
        num_Rpeaks_task(fileNum) = Output.Task.NrEvents;
        
        clear Output
    end
        
end

end

function Out=load_spikes(sessions_info,subfield,subfolder,namepart,per,varname)
condition_labels = {'Rest', 'Task'};
Out = cell(length(sessions_info), 1);
for i = 1:length(sessions_info)
    disp(['Session ' num2str(i) ' out of ' num2str(length(sessions_info))])
    if isempty(subfolder) %% unfortunate inconsistent naming
        monkey=sessions_info(i).Monkey(1:3);
    else
        monkey=['_' sessions_info(i).Monkey(1:3)];
    end
    warning('There is a hardcoded filepath')
%     results_folder = fullfile(sessions_info(i).(subfield),per,subfolder);
    results_folder = 'Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_TaskRest_Magnus_merged\per_unit_selected_stableTaskAndRest_600';
    to_load=dir(fullfile(results_folder, [namepart monkey '_' sessions_info(i).Date '*.mat']));
    if ~isempty(to_load)
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
Out = Out(~cellfun('isempty',Out));
end