function ecg_bna_get_unit_list_ecg_params(cfg)

list_of_lists = dir([cfg.unit_lists filesep '*.mat']);

list_of_lists_600 = dir([cfg.unit_lists filesep '*_600.mat']); % find files after exclusion

lists2take = ~ismember({list_of_lists.name}, {list_of_lists_600.name});

list_of_lists = list_of_lists(lists2take);

for listNum = 1:length(list_of_lists)
    
    filename = [list_of_lists(listNum).folder filesep list_of_lists(listNum).name];
    
    if contains(filename, 'unitInfo_after_exclusion.')
        % stable for at least one task block OR one rest block
        load(filename, 'unit_ids', 'targets')
    
        [num_Rpeaks_rest, num_Rpeaks_task] = count_Rpeaks_per_unit(cfg,unit_ids);
        ids_enough_Rpeaks = num_Rpeaks_rest > cfg.spk.unit_exclusion.nCardiacCycles | ...
            num_Rpeaks_task > cfg.spk.unit_exclusion.nCardiacCycles;
        
    elseif contains(filename, 'unitInfo_after_exclusion_stableTaskAndRest')
        % stable for both task and rest
        load(filename, 'unit_ids', 'targets')
    
        [num_Rpeaks_rest, num_Rpeaks_task] = count_Rpeaks_per_unit(cfg,unit_ids);
        ids_enough_Rpeaks = num_Rpeaks_rest > cfg.spk.unit_exclusion.nCardiacCycles & ...
            num_Rpeaks_task > cfg.spk.unit_exclusion.nCardiacCycles;
        
    else
        
        continue
        
    end
    
    unit_ids = unit_ids(ids_enough_Rpeaks);
    targets  = targets(ids_enough_Rpeaks);
    save([filename(1:end-4) '_600.mat'], 'unit_ids', 'targets')
    % 4. save the table with unit counts
    [unique_areas, ~, ic] = unique(targets);
    unit_counts = accumarray(ic,1);
    T = table(unique_areas, unit_counts);
    writetable(T, [filename(1:end-4) '_600.xls'])
    %     clear T unique_areas ic unit_counts filename unit_ids_after_exclusion targets_after_exclusion
    
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