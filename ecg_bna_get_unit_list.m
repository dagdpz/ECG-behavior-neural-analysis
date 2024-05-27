function keys=ecg_bna_get_unit_list(cfg,compute_unit_subsets)
% apply exclusion criteria, save lists of included and excluded units
% requires the processed part of cardioballistic analysis

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
write_unit_table(targets, filename);
clear unit_ids targets sites

if compute_unit_subsets
    % II. Unit list for BOTH task and rest
    % 0. Set up exclusion crtieria
    SNR       = keys.tt.avg_SNR;
    stability = keys.tt.avg_stability;
    FR        = keys.cal.FR;
    
    keys.tt.avg_stability = [-Inf Inf]; % drop stability thresholds, we'll implement it from 'criteria' in the neuronal data
    keys.tt.avg_SNR       = [-Inf Inf]; % drop SNR thresholds, we'll implement it from 'criteria' in the neuronal data
    keys.tt.tasktypes     = {''};       % drop conditions here, we'll implement exclusion by those in the next step
    keys.tt.FR            = [0 Inf];% drop FR, we'll use it on the next step
    keys.monkey           = ''; %% empty because we ignore which monkey it is basically
    % 1. Start with units with at least 1 spike (non-empty clusters)
    keys.tt.n_spikes      = keys.cal.n_spikes; % add 
    keys                  = ph_tuning_table_correction(keys);
    % 1.1. create unit list
    unit_ids = keys.tuning_table(2:end, unit_column_index);
    targets  = keys.tuning_table(2:end, target_column_index);
    sites    = keys.tuning_table(2:end, site_column_index);
    filename = [cfg.unit_lists filesep 'unitInfo_after_spike_exclusion_stable'];
    save(filename, 'unit_ids', 'targets', 'sites')
    % 1.2. [passed] create and save unit count table
    write_unit_table(targets, filename);
    % 1.3. figure out units excluded by number of spikes and save those
    unit_ids_after_spike_exclusion_stable = unit_ids;
    targets_after_spike_exclusion_stable  = targets;
    sites_after_spike_exclusion_stable    = sites;
    % find excluded units
    [unit_ids, ia] = setdiff(unit_ids_before_exclusion, unit_ids_after_spike_exclusion_stable); % 1 - before current exclusion step, 2 - after
    targets  = targets_before_exclusion(ia);
    sites    = sites_before_exclusion(ia);
    filename = [cfg.unit_lists filesep 'unitInfo_excluded_by_spike_stable'];
    save(filename, 'unit_ids', 'targets', 'sites')
    % 1.4. [excluded] create and save unit count table
    write_unit_table(targets, filename);
    clear unit_ids targets sites
    
    % 2. Exclude units not having both conditions
    keys.tt.tasktypes     = {'Fsac_opt','Vsac_opt'}; % WE NEED BOTH CONDITIONS HERE
    keys                  = ph_tuning_table_correction(keys);
    % 2.1. create unit list
    unit_ids = keys.tuning_table(2:end, unit_column_index);
    targets  = keys.tuning_table(2:end, target_column_index);
    sites    = keys.tuning_table(2:end, site_column_index);
    filename = [cfg.unit_lists filesep 'unitInfo_after_condition_exclusion_stable'];
    save(filename, 'unit_ids', 'targets', 'sites')
    % 2.2. [passed] create and save unit count table
    write_unit_table(targets, filename);
    % 2.3. figure out units excluded by number of spikes and save those
    unit_ids_after_condition_exclusion_stable = unit_ids;
%     targets_after_condition_exclusion_stable  = targets;
%     sites_after_condition_exclusion_stable    = sites;
    % find excluded units
    [unit_ids, ia] = setdiff(unit_ids_after_spike_exclusion_stable, unit_ids_after_condition_exclusion_stable); % 1 - before current exclusion step, 2 - after
    targets  = targets_after_spike_exclusion_stable(ia);
    sites    = sites_after_spike_exclusion_stable(ia);
    filename = [cfg.unit_lists filesep 'unitInfo_excluded_by_condition_stable'];
    save(filename, 'unit_ids', 'targets', 'sites')
    % 2.4. [excluded] create and save unit count table
    write_unit_table(targets, filename);
    clear unit_ids targets sites
    
    % 3. Exclude by FR
    keys.tt.FR            = FR;
    keys                  = ph_tuning_table_correction(keys);
    % 3.1. create unit list
    unit_ids = keys.tuning_table(2:end, unit_column_index);
    targets  = keys.tuning_table(2:end, target_column_index);
    sites    = keys.tuning_table(2:end, site_column_index);
    filename = [cfg.unit_lists filesep 'unitInfo_after_FR_exclusion_stable'];
    save(filename, 'unit_ids', 'targets', 'sites')
    % 3.2. [passed] create and save unit count table
    write_unit_table(targets, filename);
    % 3.3. figure out units excluded by FR and save those
    unit_ids_after_FR_exclusion_stable = unit_ids;
    targets_after_FR_exclusion_stable  = targets;
    sites_after_FR_exclusion_stable    = sites;
    % find excluded units
    [unit_ids, ia] = setdiff(unit_ids_after_spike_exclusion_stable, unit_ids_after_FR_exclusion_stable);
    targets  = targets_before_exclusion(ia);
    sites    = sites_before_exclusion(ia);
    filename = [cfg.unit_lists filesep 'unitInfo_excluded_by_FR_stable'];
    save(filename, 'unit_ids', 'targets', 'sites')
    % 3.4. [excluded] create and save unit count table
    write_unit_table(targets, filename);
    clear unit_ids targets sites
    
    % 4. Exclude by SNR
    % 4.0. load 'criteria' structure for units that passed FR cirterion
    curr_unit_ids = unit_ids_after_FR_exclusion_stable;
    dt = ecg_bna_load_variables(cfg, curr_unit_ids, 'cardioballistic', 'data', {});
    % 4.1. [passed] exclude by SNR
    task_SNR_ids  = dt.criteria.SNR_F >= SNR(1) & dt.criteria.SNR_F <= SNR(2);
    rest_SNR_ids  = dt.criteria.SNR_V >= SNR(1) & dt.criteria.SNR_V <= SNR(2);
    SNR_ids       = task_SNR_ids & rest_SNR_ids;
    unit_ids      = unit_ids_after_FR_exclusion_stable(SNR_ids);
    targets       = targets_after_FR_exclusion_stable(SNR_ids);
    sites         = sites_after_FR_exclusion_stable(SNR_ids);
    filename = [cfg.unit_lists filesep 'unitInfo_after_SNR_exclusion_stable'];
    save(filename, 'unit_ids', 'targets', 'sites')
    % 4.2. [passed] save table
    write_unit_table(targets, filename);
    % 4.3. figure out units excluded by SNR and save those
    unit_ids_after_SNR_exclusion_stable = unit_ids;
    targets_after_SNR_exclusion_stable  = targets;
    sites_after_SNR_exclusion_stable    = sites;
    
    % find excluded units
    [unit_ids, ia] = setdiff(unit_ids_after_FR_exclusion_stable, unit_ids_after_SNR_exclusion_stable);
    targets  = targets_after_FR_exclusion_stable(ia);
    sites    = sites_after_FR_exclusion_stable(ia);
    filename = [cfg.unit_lists filesep 'unitInfo_excluded_by_SNR_stable'];
    save(filename, 'unit_ids', 'targets', 'sites')
    % 4.4. [excluded] create and save unit count table
    write_unit_table(targets, filename);
    clear unit_ids targets sites
    
    % 5. Exclude by stability
    % 5.0. load 'criteria' structure for units that passed FR crIterion
%     curr_unit_ids = unit_ids_after_SNR_exclusion_stable;
%     dt = ecg_bna_load_variables(cfg, curr_unit_ids, 'cardioballistic', 'data', {});
    % 5.1. [passed] exclude by stability
    task_stab_ids = dt.criteria.stability_F >= stability(1) & dt.criteria.stability_F <= stability(2);
    rest_stab_ids = dt.criteria.stability_V >= stability(1) & dt.criteria.stability_V <= stability(2);
    stability_ids = task_stab_ids & rest_stab_ids;
    unit_ids      = unit_ids_after_FR_exclusion_stable(SNR_ids & stability_ids); % take list after FR exclusion and excluded both SNR and stability to load criteria only once
    targets       = targets_after_FR_exclusion_stable(SNR_ids & stability_ids);
    sites         = sites_after_FR_exclusion_stable(SNR_ids & stability_ids);
    filename = [cfg.unit_lists filesep 'unitInfo_after_stability_exclusion_stable'];
    save(filename, 'unit_ids', 'targets', 'sites')
    % 5.2. write_unit_table(targets, filename);
    write_unit_table(targets, filename);
    % 5.3. figure out units excluded by stability and save those
    unit_ids_after_stability_exclusion_stable = unit_ids;
    % find excluded units
    [unit_ids, ia] = setdiff(unit_ids_after_SNR_exclusion_stable, unit_ids_after_stability_exclusion_stable);
    targets  = targets_after_SNR_exclusion_stable(ia);
    sites    = sites_after_SNR_exclusion_stable(ia);
    filename = [cfg.unit_lists filesep 'unitInfo_excluded_by_stability_stable'];
    save(filename, 'unit_ids', 'targets', 'sites')
    % 5.4. [excluded] create and save unit count table
    write_unit_table(targets, filename);
    clear unit_ids targets sites
    
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

function write_unit_table(targets, filename)
[unique_areas, ~, ic] = unique(targets);
unit_counts = accumarray(ic,1);
T = table(unique_areas, unit_counts);
writetable(T, [filename '.xls'])
end
