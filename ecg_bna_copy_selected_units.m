function ecg_bna_copy_selected_units(unit_list, destination_folder, cfg)

load([cfg.unit_lists filesep unit_list], 'unit_ids');

if ~exist(destination_folder, 'dir')
    mkdir(destination_folder)
end

for unitNum = 1:length(unit_ids)
    
    fileList = dir([cfg.per_session_folder filesep '*' unit_ids{unitNum} '*']);
    
    if isempty(fileList)
        disp(['No files found for ' unit_ids{unitNum}])
    end
    
    for flNum = 1:length(fileList)
        
        [success,message] = ...
            copyfile([fileList(flNum).folder filesep fileList(flNum).name], destination_folder);
        if success
            disp(['Success: ' fileList(flNum).name '  --->  ' destination_folder])
            fprintf('\n')
        else
            disp('Copying failed ...')
        end
    end
end
