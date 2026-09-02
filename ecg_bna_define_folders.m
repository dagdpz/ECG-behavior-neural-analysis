function cfg = ecg_bna_define_folders(cfg)
%% simply defines fodlers and subfolders and creates them

cfg.fldr.ECG_root      = fullfile(cfg.fldr.root, 'ECG', num2str(cfg.version));
cfg.fldr.LFP_root      = fullfile(cfg.fldr.root, 'LFP', num2str(cfg.version));
cfg.fldr.MUA_root      = fullfile(cfg.fldr.root, 'MUA', num2str(cfg.version));
cfg.fldr.SPK_root      = fullfile(cfg.fldr.root, 'SPK', num2str(cfg.version));

fldr_FNs=fieldnames(cfg.fldr);
for f=1:numel(fldr_FNs)
    FN=fldr_FNs{f};
    if ~exist(cfg.fldr.(FN), 'dir')
        mkdir(cfg.fldr.(FN));
    end
end



cfg.fldr.seeds         = fullfile(cfg.fldr.ECG_root,'seeds');
cfg.fldr.LFP_sessions  = fullfile(cfg.fldr.LFP_root, 'Per_Session');
cfg.fldr.LFP_population= fullfile(cfg.fldr.LFP_root, 'Population');
cfg.fldr.LFP_sites     = fullfile(cfg.fldr.LFP_root, 'Per_Site');
cfg.fldr.MUA_sessions  = fullfile(cfg.fldr.MUA_root, 'Per_Session');
cfg.fldr.MUA_sites     = fullfile(cfg.fldr.MUA_root, 'Per_Site');
cfg.fldr.MUA_population= fullfile(cfg.fldr.MUA_root, 'Population');


cfg.fldr.unit_lists                 = fullfile(cfg.fldr.SPK_root, 'unit_lists'); % store unit lists before and after exclusion criteria
cfg.fldr.per_session_folder         = fullfile(cfg.fldr.SPK_root, 'per_unit'); % r-peak triggered psths per unit (data and plots)
cfg.fldr.per_session_selected       = fullfile(cfg.fldr.SPK_root, 'per_unit_selected_600'); % all units with 1 block of either task or rest, included by recording quality and # of R-peaks
cfg.fldr.per_session_stable         = fullfile(cfg.fldr.SPK_root, 'per_unit_stable_600'); % units that have both task and rest, included by recording quality and # of R-peaks
cfg.fldr.cardioballistic_folder     = fullfile(cfg.fldr.SPK_root, 'cardioballistic'); % cardioballistic analysis (data and plots)
cfg.fldr.cardioballistic_selected   = fullfile(cfg.fldr.SPK_root, 'cardioballistic_selected_600'); % cardioballistic: all units with 1 block of either task or rest, included by recording quality and # of R-peaks
cfg.fldr.cardioballistic_stable     = fullfile(cfg.fldr.SPK_root, 'cardioballistic_stable_600'); % cardioballistic: units that have both task and rest, included by recording quality and # of R-peaks
cfg.fldr.population_all             = fullfile(cfg.fldr.SPK_root, 'Population_AllUnits'); % population folder for all units with 1 block of either task or rest
cfg.fldr.population_bothTaskRest    = fullfile(cfg.fldr.SPK_root, 'Population_bothTaskRest'); % population folder for units that have both task and rest
cfg.fldr.population_cardioballistic = fullfile(cfg.fldr.SPK_root, 'Population_cardioballistic');


cfg.fldr.LFP_noise= cfg.fldr.LFP_root;                            % folder to save noise rejection results
cfg.fldr.ECG_processed     = [cfg.fldr.ECG_root filesep 'Processed ECG'];  % folder to save ecg processing results
cfg.fldr.LFP_processed     = [cfg.fldr.LFP_root filesep 'Processed LFP'];  % folder to save lfp processing results

fldr_FNs=fieldnames(cfg.fldr);
for f=1:numel(fldr_FNs)
    FN=fldr_FNs{f};
    if ~exist(cfg.fldr.(FN), 'dir')
        mkdir(cfg.fldr.(FN));
    end
end


% save settings struct
save(fullfile(cfg.fldr.ECG_root, ['settings_' num2str(cfg.version) '.mat']), 'cfg');
end

