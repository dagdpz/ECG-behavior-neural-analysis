function sessions_avg = ecg_bna_avg_sessions_Rpeak_evoked_ECG(Rpeak_evoked_ECG, ecg_bna_cfg)
%ecg_bna_avg_sessions_Rpeak_evoked_ECG  - Rpeak evoked LFP response
% average across many sessions
%
% USAGE:
%	sessions_avg = ecg_bna_avg_sessions_Rpeak_evoked_ECG(Rpeak_evoked_ECG, ecg_bna_cfg)
%
% INPUTS:
%		Rpeak_evoked_ECG    - struct containing the Rpeak evoked ECG
%       response average for indiviual sessions, output of 
%       ecg_bna_compute_session_evoked_ECG.m
%		ecg_bna_cfg         - struct containing the required settings
%           Required Fields:
%               conditions          - trial conditions to compare, see
%               lfp_tfa_settings.m and lfp_tfa_compare_conditions.m
%               root_results_fldr   - root folder where results are saved
%               compare.targets     - targets to compare, see lfp_tfa_settings.m
% OUTPUTS:
%		sessions_avg    - structure containing Rpeak evoked ECG
%		response averaged across multiple sessions
%
% REQUIRES:	lfp_tfa_plot_evoked_lfp
%
% See also lfp_tfa_settings, lfp_tfa_define_settings, lfp_tfa_compare_conditions, 
% lfp_tfa_plot_site_evoked_LFP
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


    % results folder
    results_fldr = fullfile(ecg_bna_cfg.analyse_ecg_folder, 'Evoked ECG');
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % Average Evoked LFP response across sessions
    sessions_avg = struct();
    t = 1;
    
    for cn = 1:length(ecg_bna_cfg.conditions)
        fprintf('Condition %s\n', ecg_bna_cfg.conditions(cn).label);
        sessions_avg(t).condition(cn).hs_tuned_evoked = struct();
        sessions_avg(t).condition(cn).cfg_condition = ecg_bna_cfg.conditions(cn);
        sessions_avg(t).condition(cn).label = ecg_bna_cfg.conditions(cn).label;
        
        % initialize number of site pairs for each handspace
        % label
        for st = 1:size(Rpeak_evoked_ECG.session(1).condition(cn).hs_tuned_evoked, 1)
            for hs = 1:size(Rpeak_evoked_ECG.session(1).condition(cn).hs_tuned_evoked, 2)
                sessions_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsessions = 0;
                sessions_avg(t).condition(cn).hs_tuned_evoked(st, hs).trial = {};                
            end
        end  

        for i = 1:length(Rpeak_evoked_ECG.session)
            if isempty(Rpeak_evoked_ECG.session(i).condition)
                continue;
            end
            if isfield(Rpeak_evoked_ECG.session(i).condition(cn).hs_tuned_evoked, 'mean')
                for st = 1:size(Rpeak_evoked_ECG.session(i).condition(cn).hs_tuned_evoked, 1)
                    for hs = 1:size(Rpeak_evoked_ECG.session(i).condition(cn).hs_tuned_evoked, 2)
                        if isfield(Rpeak_evoked_ECG.session(i).condition(cn).hs_tuned_evoked(st, hs), 'mean') ...
                                && ~isempty(Rpeak_evoked_ECG.session(i).condition(cn).hs_tuned_evoked(st, hs).mean)
                            sessions_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsessions = ...
                                sessions_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsessions + 1;
                            if sessions_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsessions == 1
                                sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).time ...
                                    = Rpeak_evoked_ECG.session(i).condition(cn).hs_tuned_evoked(st, hs).time;
                                sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).hs_label ...
                                    = Rpeak_evoked_ECG.session(i).condition(cn).hs_tuned_evoked(st, hs).hs_label;
                                if isfield(Rpeak_evoked_ECG.session(i).condition(cn).hs_tuned_evoked(st, hs), 'state') && ...
                                        isfield(Rpeak_evoked_ECG.session(i).condition(cn).hs_tuned_evoked(st, hs), 'state_name')
                                    sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).state ...
                                        = Rpeak_evoked_ECG.session(i).condition(cn).hs_tuned_evoked(st, hs).state;
                                    sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).state_name ...
                                        = Rpeak_evoked_ECG.session(i).condition(cn).hs_tuned_evoked(st, hs).state_name;
                                end
                            end
                            sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial ...
                                = [sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial, ...
                                Rpeak_evoked_ECG.session(i).condition(cn).hs_tuned_evoked(st, hs).mean];  
                        else
                            continue; %sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs) = struct();
                        end
                    end
                end
            end                               
        end

        % compute average
        if isfield(sessions_avg(t).condition(cn).hs_tuned_evoked, 'trial')
            for st = 1:size(sessions_avg(t).condition(cn).hs_tuned_evoked, 1)
                for hs = 1:size(sessions_avg(t).condition(cn).hs_tuned_evoked, 2)
                    if ~isempty(sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial)
                        sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial = ...
                            cat(1, sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial{:});
                        sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).dimord = 'nsessions_time';
%                         sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt(...
%                             :, isnan(sum(sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt, 1))) = nan;
                        sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).std = ...
                            std(sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial, 0, 1);
                        sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).mean = ...
                            mean(sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial, 1);                           
                    end
                end
            end
        end


        if ~isempty(sessions_avg(t).condition(cn).hs_tuned_evoked)
            if isfield(sessions_avg(t).condition(cn).hs_tuned_evoked,... 
                    'mean')
                plottitle = [ecg_bna_cfg.compare.targets{t},...
                     ecg_bna_cfg.conditions(cn).label];
                result_file = fullfile(results_fldr, ...
                                ['ECG_b2bt_Evoked_' ecg_bna_cfg.conditions(cn).label]);
                ecg_bna_plot_evoked_lfp(sessions_avg(t).condition(cn).hs_tuned_evoked, ...
                            ecg_bna_cfg, plottitle, result_file, 'ylabel', 'ECG amplitude');
            end
        end
        
    end

    % difference between conditions
    sessions_avg(t).difference = [];
    for diff = 1:size(ecg_bna_cfg.diff_condition, 2)
        diff_condition = ecg_bna_cfg.diff_condition{diff};
        diff_color = []; diff_legend = [];
        if isfield(ecg_bna_cfg, 'diff_color')
            diff_color = ecg_bna_cfg.diff_color{diff};
        end
        if isfield(ecg_bna_cfg, 'diff_legend')
            diff_legend = ecg_bna_cfg.diff_legend{diff};
        end
        sessions_avg(t).difference = [sessions_avg(t).difference, ...
            ecg_bna_compute_diff_condition_average('Event_evoked_ECG_R2Rt', ...
            sessions_avg(t).condition, diff_condition, diff_color,diff_legend)];
    end
    % plot Difference TFR
    for dcn = 1:length(sessions_avg(t).difference)
        if ~isempty(sessions_avg(t).difference(dcn).hs_tuned_evoked)
            if isfield(sessions_avg(t).difference(dcn).hs_tuned_evoked,... 
                    'mean')
                plottitle = ['Target ', ecg_bna_cfg.compare.targets{t}, ...
                    sessions_avg(t).difference(dcn).label];
                result_file = fullfile(results_fldr, ...
                    ['ECG_b2bt_DiffEvoked_' 'diff_condition' num2str(dcn)]);
                    %sessions_avg(t).difference(dcn).label '.png']);
                ecg_bna_plot_evoked_lfp(sessions_avg(t).difference(dcn).hs_tuned_evoked, ...
                    ecg_bna_cfg, plottitle, result_file, 'ylabel', 'ECG amplitude');
            end
        end
    end
    
    % save session average tfs
    save(fullfile(results_fldr, 'sessions_avg_Rpeak_evoked_LFP.mat'), 'sessions_avg');
        
    close all;
end