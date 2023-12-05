function keys=ecg_bna_get_unit_list(cfg,compute_unit_subsets)

%% apply exclusion criteria, save lists of included and excluded units
keys=struct;
keys=ph_general_settings(cfg.project,keys);
project_specific_settings=[keys.db_folder 'ph_project_settings.m'];
run(project_specific_settings)
version_specific_settings=[keys.db_folder cfg.spikes_version filesep 'ph_project_version_settings.m'];
run(version_specific_settings)
keys.anova_table_file=[keys.basepath_to_save cfg.spikes_version filesep 'tuning_table_combined_CI.mat'];
keys.tuning_table = ph_load_tuning_table(keys); %% load tuning table
%keys.tuning_table(2:end,3) = cellfun(@(x) x(1:end-2), keys.tuning_table(2:end,3), 'UniformOutput', false);

% I. Unit list before exclusion - this one doesn't take number of
% R-peaks into account
% create unit list before exclusion criteria and save it
unit_ids_before_exclusion = keys.tuning_table(2:end,1);
targets_before_exclusion  = keys.tuning_table(2:end,3);
filename = [cfg.SPK_root_results_fldr filesep 'unitInfo_before_exclusion'];
save(filename, 'unit_ids_before_exclusion', 'targets_before_exclusion')
% create table of unit numbers by area
[unique_areas, ~, ic] = unique(targets_before_exclusion);
unit_counts = accumarray(ic,1);
T = table(unique_areas, unit_counts);
writetable(T, [filename '.xls'])
%     clear T unique_areas ic unit_counts filename

if compute_unit_subsets
    % II. Unit list for BOTH task and rest
    % 1. apply exclusion criteria - choose only cells for BOTH task and rest
    keys.tt.tasktypes= {'Fsac_opt','Vsac_opt'};
    keys.tt.FR=keys.cal.FR;
    keys.tt.n_spikes=keys.cal.n_spikes;
    keys.monkey=''; %% empty because we ignore which monkey it is basically
    keys=ph_tuning_table_correction(keys);
    unit_ids_after_exclusion = keys.tuning_table(2:end,1);
    targets_after_exclusion  = keys.tuning_table(2:end,3);
    filename = [cfg.SPK_root_results_fldr filesep 'unitInfo_after_exclusion_stableTaskAndRest'];
    save(filename, 'unit_ids_after_exclusion', 'targets_after_exclusion')
    % 2. create table of unit numbers by area
    [unique_areas, ~, ic] = unique(targets_after_exclusion);
    unit_counts = accumarray(ic,1);
    T = table(unique_areas, unit_counts);
    writetable(T, [filename '.xls'])
    % 3. exclude additionally by the number of R-peaks
    [num_Rpeaks_rest, num_Rpeaks_task] = count_Rpeaks_per_unit(cfg,unit_ids_after_exclusion);
    ids_enough_Rpeaks = num_Rpeaks_rest > cfg.spk.unit_exclusion.nCardiacCycles & ...
        num_Rpeaks_task > cfg.spk.unit_exclusion.nCardiacCycles;
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
    filename = [cfg.SPK_root_results_fldr filesep 'unitInfo_excluded_stableTaskAndRest'];
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
    [num_Rpeaks_rest, ~] = count_Rpeaks_per_unit(cfg,unit_ids_after_exclusion_rest);
    ids_enough_Rpeaks = num_Rpeaks_rest > cfg.spk.unit_exclusion.nCardiacCycles;
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
    [~, num_Rpeaks_task] = count_Rpeaks_per_unit(cfg,unit_ids_after_exclusion_task);
    ids_enough_Rpeaks = num_Rpeaks_task > cfg.spk.unit_exclusion.nCardiacCycles;
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
    
    filename = [cfg.SPK_root_results_fldr filesep 'unitInfo_after_exclusion'];
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
    filename = [cfg.SPK_root_results_fldr filesep 'unitInfo_excluded'];
    save(filename, 'unit_ids_excluded', 'targets_excluded')
    [unique_areas, ~, ic] = unique(targets_excluded);
    unit_counts = accumarray(ic,1);
    T = table(unique_areas, unit_counts);
    writetable(T, [filename '.xls'])
    %     clear T unique_areas ic unit_counts filename unit_ids_excluded targets_excluded unit_ids_excluded targets_excluded unit_ids_before_exclusion targets_before_exclusion
end

end


function [num_Rpeaks_rest, num_Rpeaks_task] = count_Rpeaks_per_unit(cfg,unit_list)

num_Rpeaks_rest = nan(length(unit_list),1);
num_Rpeaks_task = nan(length(unit_list),1);

for fileNum = 1:length(unit_list)
    folder=[cfg.SPK_root_results_fldr filesep 'per_unit' filesep];
    currFileInfo = dir([folder unit_list{fileNum} '_*.mat']);
    
    if ~isempty(currFileInfo)
        
        disp(['Reading file ' num2str(fileNum) ' out of ' num2str(length(unit_list))])
        load([folder filesep currFileInfo.name], 'Output')
        
        num_Rpeaks_rest(fileNum) = Output.Rest.NrEvents;
        num_Rpeaks_task(fileNum) = Output.Task.NrEvents;
        
        clear Output
    end
    
end

end
