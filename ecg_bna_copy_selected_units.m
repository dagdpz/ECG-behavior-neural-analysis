function ecg_bna_copy_selected_units(cfg)

exclusion_lists = dir([cfg.SPK_root_results_fldr filesep '*unitInfo*.mat']);

from_folder = [cfg.SPK_root_results_fldr filesep 'per_unit\'];

%%
for listNum = 1:length(exclusion_lists)
    
    unit_list = load([exclusion_lists(listNum).folder filesep exclusion_lists(listNum).name]);
    
    fns = fieldnames(unit_list);
    
    if contains(exclusion_lists(listNum).name, 'stableTaskAndRest')
        if contains(exclusion_lists(listNum).name, '_600')
            to_folder = [from_folder(1:end-1) '_selected_stableTaskAndRest_600'];
        elseif contains(exclusion_lists(listNum).name, 'after_exclusion')
            to_folder = [from_folder(1:end-1) '_selected_stableTaskAndRest'];
        else
            to_folder = [from_folder(1:end-1) '_excluded_stableTaskAndRest'];
        end
    elseif contains(exclusion_lists(listNum).name, 'after_exclusion')
        to_folder = [from_folder(1:end-1) '_selected'];
    elseif contains(exclusion_lists(listNum).name, 'excluded')
        to_folder = [from_folder(1:end-1) '_excluded'];
    elseif contains(exclusion_lists(listNum).name, 'before_exclusion')
        continue
    end
    
    if ~exist(to_folder, 'dir')
        mkdir(to_folder)
    end
    
    for unitNum = 1:length(unit_list.(fns{1}))
        
        fileList = dir([from_folder '*' unit_list.(fns{1}){unitNum} '*']);
        
        if isempty(fileList)
            disp(['No files found for ' unit_list.(fns{1}){unitNum}])
        end
        
        for flNum = 1:length(fileList)
        
            [success,message] = ...
                copyfile([fileList(flNum).folder filesep fileList(flNum).name], to_folder);
            if success
                disp(['Success: ' fileList(flNum).name '  --->  ' to_folder])
                fprintf('\n')
            else
                disp('Copying failed ...')
            end
        end
    end
end
