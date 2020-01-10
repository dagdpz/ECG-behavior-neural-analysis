function ecg_bna_cfg = ecg_bna_define_settings(settings_filepath)
%lfp_tfa_define_settings - Function to define LFP time frequency analysis settings 
%
% USAGE:
%	lfp_tfa_cfg = lfp_tfa_define_settings(settings_filepath, maxsites)
%
% INPUTS:
%       settings_filepath         - absolute path to the matlab script file
%       where LFP TFA settings are defined, see settings/lfp_tfa_settings
%       maxsites                  - maximum number of sites to be analysed per session,
%       set to infinity to analyse all sites
%
% OUTPUTS:
%		lfp_tfa_cfg               - structure containing all settings
%
% REQUIRES:	lfp_tfa_read_info_file, lfp_tfa_compare_conditions,
% lfp_tfa_define_states, lfp_tfa_define_epochs
%
% See also settings/lfp_tfa_settings, lfp_tfa_read_info_file, lfp_tfa_compare_conditions,
% lfp_tfa_define_states, lfp_tfa_define_epochs
%
% Author(s):	S.Nair, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2019-02-15:	Created function (Sarath Nair)
% 2019-03-05:	First Revision
% ...
% $Revision: 1.0 $  $Date: 2019-03-05 17:18:00 $

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%

    % add functions from lfp analysis pipeline
    addpath('C:\Source\MATLAB\LFP_timefrequency_analysis');
    
    % add external functions to path
    addpath(genpath('./external'));
    
    % define state IDs
    lfp_tfa_global_define_states;    

    % load the specified settings file
    run(settings_filepath);
    
    % read info excel file (Sorted neurons file)
    ecg_bna_cfg.sites_info = lfp_tfa_read_info_file(ecg_bna_cfg);
     
    % create a root folder to save results of the analysis
    % root_results_folder = [ecg_tfa_cfg.results_folder,
    % ecg_tfa_cfg.version];
    ecg_bna_cfg.root_results_fldr = fullfile(ecg_bna_cfg.results_folder, ...
        num2str(ecg_bna_cfg.version));
    if ~exist(ecg_bna_cfg.root_results_fldr, 'dir')
        mkdir(ecg_bna_cfg.root_results_fldr);
    end
    
    % get conditions to be included in the analysis
    ecg_bna_cfg.conditions = lfp_tfa_compare_conditions(ecg_bna_cfg);
    
    % folder to save noise rejection results
    ecg_bna_cfg.noise.results_folder = ecg_bna_cfg.root_results_fldr;
    % folder to save baseline results
    ecg_bna_cfg.results_folder = ecg_bna_cfg.root_results_fldr;
    
    % folder to save ecg processing results
    ecg_bna_cfg.proc_ecg_folder = [ecg_bna_cfg.root_results_fldr filesep 'Processed ECG'];
    % folder to save ecg analysis results
    ecg_bna_cfg.analyse_ecg_folder = [ecg_bna_cfg.root_results_fldr filesep 'ECG Analysis'];
    % folder to store session-wise analysis results
    for i = 1:length(ecg_bna_cfg.session_info)
        ecg_bna_cfg.session_info(i).session = ...
            [ecg_bna_cfg.session_info(i).Monkey, '_', ecg_bna_cfg.session_info(i).Date];
        ecg_bna_cfg.session_info(i).proc_results_fldr = ...
                fullfile(ecg_bna_cfg.proc_ecg_folder, ecg_bna_cfg.session_info(i).session);
    end

    % save settings struct
    save(fullfile(ecg_bna_cfg.root_results_fldr, ['ecg_bna_settings_' num2str(ecg_bna_cfg.version) '.mat']), ...
        'ecg_bna_cfg');

end

