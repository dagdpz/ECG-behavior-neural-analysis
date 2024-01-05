function cfg = ecg_bna_define_folders(cfg)
%ecg_bna_define_settings - Function to define ECG related bahvior and
%neural analysis settings
%
% USAGE:
%	ecg_bna_cfg = ecg_bna_define_settings(settings_filepath)
%
% INPUTS:
%       settings_filepath         - absolute path to the matlab script file
%       where LFP TFA settings are defined, for example,
%       [path_to_pipeline]\settings\Magnus\ecg_bna_settings_Magnus_dPul_ECG_LFP.m
%
% OUTPUTS:
%		ecg_bna_cfg               - structure containing all settings
%
% REQUIRES:	lfp_tfa_global_define_states, lfp_tfa_read_info_file,
% lfp_tfa_compare_conditions
%
% See also lfp_tfa_define_settings
%
% Author(s):	S.Nair, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%

% define state IDs
%lfp_tfa_global_define_states;

% load the specified settings file

cfg.ECG_root_results_fldr      = fullfile(cfg.results_folder, 'ECG', num2str(cfg.version));
cfg.LFP_root_results_fldr      = fullfile(cfg.results_folder, 'LFP', num2str(cfg.version));
cfg.SPK_root_results_fldr      = fullfile(cfg.results_folder, 'ECG_triggered_spikes', num2str(cfg.version));
cfg.unit_lists                 = fullfile(cfg.SPK_root_results_fldr, 'unit_lists'); % store unit lists before and after exclusion criteria
cfg.per_session_folder         = fullfile(cfg.SPK_root_results_fldr, 'per_unit'); % r-peak triggered psths per unit (data and plots)
cfg.per_session_selected       = fullfile(cfg.SPK_root_results_fldr, 'per_unit_selected_600'); % all units with 1 block of either task or rest, included by recording quality and # of R-peaks
cfg.per_session_stable         = fullfile(cfg.SPK_root_results_fldr, 'per_unit_stable_600'); % units that have both task and rest, included by recording quality and # of R-peaks
cfg.cardioballistic_folder     = fullfile(cfg.SPK_root_results_fldr, 'cardioballistic'); % cardioballistic analysis (data and plots)
cfg.cardioballistic_selected   = fullfile(cfg.SPK_root_results_fldr, 'cardioballistic_selected_600'); % cardioballistic: all units with 1 block of either task or rest, included by recording quality and # of R-peaks
cfg.cardioballistic_stable     = fullfile(cfg.SPK_root_results_fldr, 'cardioballistic_stable_600'); % cardioballistic: units that have both task and rest, included by recording quality and # of R-peaks
cfg.population_all             = fullfile(cfg.SPK_root_results_fldr, 'Population_AllUnits'); % population folder for all units with 1 block of either task or rest
cfg.population_bothTaskRest    = fullfile(cfg.SPK_root_results_fldr, 'Population_bothTaskRest'); % population folder for units that have both task and rest
cfg.population_cardioballistic = fullfile(cfg.SPK_root_results_fldr, 'Population_cardioballistic');
if ~exist(cfg.ECG_root_results_fldr, 'dir')
    mkdir(cfg.ECG_root_results_fldr);
end
if ~exist(cfg.LFP_root_results_fldr, 'dir')
    mkdir(cfg.LFP_root_results_fldr);
end
if ~exist(cfg.SPK_root_results_fldr, 'dir')
    mkdir(cfg.SPK_root_results_fldr);
end
if ~exist(cfg.unit_lists, 'dir')
    mkdir(cfg.unit_lists);
end
if ~exist(cfg.per_session_folder, 'dir')
    mkdir(cfg.per_session_folder);
end
if ~exist(cfg.cardioballistic_folder, 'dir')
    mkdir(cfg.cardioballistic_folder);
end
if ~exist(cfg.population_all, 'dir')
    mkdir(cfg.population_all);
end
if ~exist(cfg.population_bothTaskRest, 'dir')
    mkdir(cfg.population_bothTaskRest);
end

% get conditions to be included in the analysis
%ecg_bna_cfg.conditions          = lfp_tfa_compare_conditions(ecg_bna_cfg);
cfg.noise.results_folder= cfg.LFP_root_results_fldr;                            % folder to save noise rejection results
cfg.results_folder      = cfg.LFP_root_results_fldr;                            % folder to save baseline results
cfg.proc_ecg_folder     = [cfg.ECG_root_results_fldr filesep 'Processed ECG'];  % folder to save ecg processing results
cfg.analyse_ecg_folder  = [cfg.ECG_root_results_fldr];                          % folder to save ecg analysis results
cfg.proc_lfp_folder     = [cfg.LFP_root_results_fldr filesep 'Processed LFP'];  % folder to save lfp processing results
cfg.analyse_lfp_folder  = [cfg.LFP_root_results_fldr];                          % folder to save lfp analysis results

% % Folder to save results of LFP processing
% if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_LFP')) || any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_TFS'))
%     if ecg_bna_cfg.process_LFP || ~exist(ecg_bna_cfg.proc_lfp_folder, 'dir')
%         ecg_bna_cfg.process_LFP = true;
%     end
% % end

% save settings struct
save(fullfile(cfg.ECG_root_results_fldr, ['settings_' num2str(cfg.version) '.mat']), 'cfg');
end

