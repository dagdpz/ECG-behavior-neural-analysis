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

cfg.ECG_root_results_fldr = fullfile(cfg.results_folder, 'ECG', num2str(cfg.version));
cfg.LFP_root_results_fldr = fullfile(cfg.results_folder, 'LFP', num2str(cfg.version));
cfg.SPK_root_results_fldr = fullfile(cfg.results_folder, 'ECG_triggered_spikes', num2str(cfg.version));
if ~exist(cfg.ECG_root_results_fldr, 'dir')
    mkdir(cfg.ECG_root_results_fldr);
end
if ~exist(cfg.LFP_root_results_fldr, 'dir')
    mkdir(cfg.LFP_root_results_fldr);
end
if ~exist(cfg.SPK_root_results_fldr, 'dir')
    mkdir(cfg.SPK_root_results_fldr);
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
% end

% folder to store session-wise analysis results - let's try to please NOT do this
% for i = 1:length(cfg.session_info)
%     cfg.session_info(i).session             = [cfg.session_info(i).Monkey, '_', cfg.session_info(i).Date];
%     cfg.session_info(i).proc_ecg_fldr       = cfg.proc_ecg_folder;
%     cfg.session_info(i).analyse_ecg_fldr    = cfg.analyse_ecg_folder;
%     cfg.session_info(i).analyse_lfp_fldr    = cfg.analyse_lfp_folder;
%     cfg.session_info(i).proc_lfp_fldr       = cfg.proc_lfp_folder;
%     cfg.session_info(i).SPK_fldr            = cfg.SPK_root_results_fldr;
% end

% save settings struct
save(fullfile(cfg.ECG_root_results_fldr, ['settings_' num2str(cfg.version) '.mat']), 'cfg');
end

