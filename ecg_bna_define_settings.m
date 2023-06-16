function ecg_bna_cfg = ecg_bna_define_settings(github_folder,project,version)
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
lfp_tfa_global_define_states;

% load the specified settings file
run([github_folder filesep 'Settings' filesep project filesep 'ECG_bna' filesep version '.m']);

ecg_bna_cfg.ECG_root_results_fldr = fullfile(ecg_bna_cfg.results_folder, 'ECG', num2str(ecg_bna_cfg.version));
ecg_bna_cfg.LFP_root_results_fldr = fullfile(ecg_bna_cfg.results_folder, 'LFP', num2str(ecg_bna_cfg.LFP_version));
ecg_bna_cfg.SPK_root_results_fldr = fullfile(ecg_bna_cfg.results_folder, 'ECG_triggered_spikes', num2str(ecg_bna_cfg.version));
if ~exist(ecg_bna_cfg.ECG_root_results_fldr, 'dir')
    mkdir(ecg_bna_cfg.ECG_root_results_fldr);
end
if ~exist(ecg_bna_cfg.LFP_root_results_fldr, 'dir')
    mkdir(ecg_bna_cfg.LFP_root_results_fldr);
end
if ~exist(ecg_bna_cfg.SPK_root_results_fldr, 'dir')
    mkdir(ecg_bna_cfg.SPK_root_results_fldr);
end

% get conditions to be included in the analysis
ecg_bna_cfg.conditions          = lfp_tfa_compare_conditions(ecg_bna_cfg);
ecg_bna_cfg.noise.results_folder= ecg_bna_cfg.LFP_root_results_fldr;                            % folder to save noise rejection results
ecg_bna_cfg.results_folder      = ecg_bna_cfg.LFP_root_results_fldr;                            % folder to save baseline results
ecg_bna_cfg.proc_ecg_folder     = [ecg_bna_cfg.ECG_root_results_fldr filesep 'Processed ECG'];  % folder to save ecg processing results
ecg_bna_cfg.analyse_ecg_folder  = [ecg_bna_cfg.ECG_root_results_fldr];                          % folder to save ecg analysis results
ecg_bna_cfg.proc_lfp_folder     = [ecg_bna_cfg.LFP_root_results_fldr filesep 'Processed LFP'];  % folder to save lfp processing results
ecg_bna_cfg.analyse_lfp_folder  = [ecg_bna_cfg.LFP_root_results_fldr];                          % folder to save lfp analysis results

% Folder to save results of LFP processing
if any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_LFP')) || any(strcmp(ecg_bna_cfg.analyses, 'Rpeak_evoked_TFS'))
    if ecg_bna_cfg.process_LFP || ~exist(ecg_bna_cfg.proc_lfp_folder, 'dir')
        ecg_bna_cfg.process_LFP = true;
    end
end

% folder to store session-wise analysis results - let's try to please NOT do this
for i = 1:length(ecg_bna_cfg.session_info)
    ecg_bna_cfg.session_info(i).session             = [ecg_bna_cfg.session_info(i).Monkey, '_', ecg_bna_cfg.session_info(i).Date];
    ecg_bna_cfg.session_info(i).proc_ecg_fldr       = ecg_bna_cfg.proc_ecg_folder;
    ecg_bna_cfg.session_info(i).analyse_ecg_fldr    = ecg_bna_cfg.analyse_ecg_folder;
    ecg_bna_cfg.session_info(i).analyse_lfp_fldr    = ecg_bna_cfg.analyse_lfp_folder;
    ecg_bna_cfg.session_info(i).proc_lfp_fldr       = ecg_bna_cfg.proc_lfp_folder;
    ecg_bna_cfg.session_info(i).SPK_fldr            = ecg_bna_cfg.SPK_root_results_fldr;
end
% save settings struct
save(fullfile(ecg_bna_cfg.ECG_root_results_fldr, ['ecg_bna_settings_' num2str(ecg_bna_cfg.version) '.mat']), 'ecg_bna_cfg');
end

