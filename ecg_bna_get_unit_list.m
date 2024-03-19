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
% find column indices
unit_column_index   = DAG_find_column_index(keys.tuning_table, 'unit_ID');
target_column_index = DAG_find_column_index(keys.tuning_table, 'target');
site_column_index   = DAG_find_column_index(keys.tuning_table, 'site_ID');
keys.tuning_table(2:end, target_column_index) = cellfun(@(x) x(1:end-2), keys.tuning_table(2:end, target_column_index), 'UniformOutput', false);

% I. Unit list before exclusion - this one doesn't take number of
% R-peaks into account
% create unit list before exclusion criteria and save it
unit_ids = keys.tuning_table(2:end, unit_column_index);
targets  = keys.tuning_table(2:end, target_column_index);
sites    = keys.tuning_table(2:end, site_column_index);

filename = [cfg.unit_lists filesep 'unitInfo_before_exclusion'];
save(filename, 'unit_ids', 'targets', 'sites')

% keep unit list before exclusion for later
unit_ids_before_exclusion = unit_ids;
targets_before_exclusion  = targets;
sites_before_exclusion    = sites;

% create table of unit numbers by area
[unique_areas, ~, ic] = unique(targets);
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
    unit_ids = keys.tuning_table(2:end, unit_column_index);
    targets  = keys.tuning_table(2:end, target_column_index);
    sites    = keys.tuning_table(2:end, site_column_index);
    filename = [cfg.unit_lists filesep 'unitInfo_after_exclusion_stableTaskAndRest'];
    save(filename, 'unit_ids', 'targets', 'sites')
    % 2. create table of unit numbers by area
    [unique_areas, ~, ic] = unique(targets);
    unit_counts = accumarray(ic,1);
    T = table(unique_areas, unit_counts);
    writetable(T, [filename '.xls'])
    
    % III. Unit list of excluded units (that aren't units existing for both task and rest)
    ids_excluded = ~ismember(unit_ids_before_exclusion, unit_ids);
    unit_ids_excluded = unit_ids_before_exclusion(ids_excluded);
    targets_excluded  = targets_before_exclusion(ids_excluded);
    sites_excluded    = sites_before_exclusion(ids_excluded);
    filename = [cfg.unit_lists filesep 'unitInfo_excluded_stableTaskAndRest'];
    unit_ids = unit_ids_excluded;
    targets  = targets_excluded;
    sites    = sites_excluded;
    save(filename, 'unit_ids', 'targets', 'sites')
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
    keys.tuning_table(2:end, target_column_index) = cellfun(@(x) x(1:end-2), keys.tuning_table(2:end, target_column_index), 'UniformOutput', false);
    keys = ph_tuning_table_correction(keys);
    unit_ids_after_exclusion_rest = keys.tuning_table(2:end, unit_column_index);
    targets_after_exclusion_rest  = keys.tuning_table(2:end, target_column_index);
    sites_after_exclusion_rest    = keys.tuning_table(2:end, site_column_index);
    
    % 2. choose units that had a rest block
    keys.tt.tasktypes= {'Vsac_opt'}; % fixation saccade - for rest
    keys.tuning_table=ph_load_tuning_table(keys); %% load tuning table
    keys.tuning_table(2:end, target_column_index) = cellfun(@(x) x(1:end-2), keys.tuning_table(2:end, target_column_index), 'UniformOutput', false);
    keys = ph_tuning_table_correction(keys);
    unit_ids_after_exclusion_task = keys.tuning_table(2:end, unit_column_index);
    targets_after_exclusion_task  = keys.tuning_table(2:end, target_column_index);
    sites_after_exclusion_task    = keys.tuning_table(2:end, site_column_index);
    
    % 3. find units that had at least 1 rest block or 1 task block, figure out
    % targets for those
    unit_ids_after_exclusion = union(unit_ids_after_exclusion_rest, unit_ids_after_exclusion_task);
    ids_rest = ismember(unit_ids_after_exclusion, unit_ids_after_exclusion_rest);
    ids_task = ismember(unit_ids_after_exclusion, unit_ids_after_exclusion_task);
    targets_after_exclusion = cell(length(unit_ids_after_exclusion), 1);
    targets_after_exclusion(ids_rest) = targets_after_exclusion_rest;
    targets_after_exclusion(ids_task) = targets_after_exclusion_task;
    sites_after_exclusion   = cell(length(unit_ids_after_exclusion), 1);
    sites_after_exclusion(ids_rest) = sites_after_exclusion_rest;
    sites_after_exclusion(ids_task) = sites_after_exclusion_task;
    
    filename = [cfg.unit_lists filesep 'unitInfo_after_exclusion'];
    unit_ids = unit_ids_after_exclusion;
    targets  = targets_after_exclusion;
    sites    = sites_after_exclusion;
    save(filename, 'unit_ids', 'targets', 'sites')
    % clear unit_ids_after_exclusion targets_after_exclusion
    [unique_areas, ~, ic] = unique(targets_after_exclusion);
    unit_counts = accumarray(ic,1);
    T = table(unique_areas, unit_counts);
    writetable(T, [filename '.xls'])
    %     clear T unique_areas ic unit_counts unit_ids_excluded targets_after_exclusion
    
    % V. Unit list of excluded units (that aren't units existing either task or rest)
    ids_excluded = ~ismember(unit_ids, unit_ids_after_exclusion);
    unit_ids_excluded = unit_ids(ids_excluded);
    targets_excluded  = targets(ids_excluded);
    sites_excluded    = sites(ids_excluded);
    unit_ids = unit_ids_excluded;
    targets  = targets_excluded;
    sites    = sites_excluded;
    filename = [cfg.unit_lists filesep 'unitInfo_excluded'];
    save(filename, 'unit_ids', 'targets', 'sites')
    [unique_areas, ~, ic] = unique(targets_excluded);
    unit_counts = accumarray(ic,1);
    T = table(unique_areas, unit_counts);
    writetable(T, [filename '.xls'])
    %     clear T unique_areas ic unit_counts filename unit_ids_excluded targets_excluded unit_ids_excluded targets_excluded unit_ids_before_exclusion targets_before_exclusion
end

end
