function sessions_avg = ecg_bna_avg_sessions_ECGb2bt_evoked(ecg_b2bt_evoked, ecg_bna_cfg)
% ecg_bna_avg_sessions_ECGb2bt_evoked  - Event evoked ECG b2bt average
% across sessions
%
% USAGE:
%	sessions_avg = ecg_bna_avg_sessions_ECGb2bt_evoked(ecg_b2bt_evoked, ecg_bna_cfg)
%
% INPUTS:
%		ecg_b2bt_evoked		- struct containing the ECG b2bt for
%		individual sessions, output of ecg_bna_compute_session_evoked_ECG_R2Rt.m
%		ecg_bna_cfg         - struct containing the required settings
%           Required Fields:
%               conditions          - trial conditions to compare,
%               condition is a combination of hand-space tuning,
%               control/inactivation, choice/trial, type-effector, and
%               success
%               root_results_fldr   - root folder where results are saved
%               compare.targets     - targets to compare, see lfp_tfa_settings.m
% OUTPUTS:
%		sessions_avg    - structure containing condition-based evoked LFP
%		response averaged across multiple sessions
%           Fields:
%              condition.hs_tuned_evoked - structure containing relevant
%              results for event evoked ECG b2bt average across sessions
%
% REQUIRES:	ecg_bna_plot_evoked_R2Rt
%
% See also ecg_bna_define_settings, ecg_bna_plot_evoked_R2Rt
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
results_fldr = fullfile(ecg_bna_cfg.analyse_ecg_folder, 'ECG R2Rt');
if ~exist(results_fldr, 'dir')
    mkdir(results_fldr);
end

% Average Evoked LFP response across sessions
sessions_avg = struct();
% this variable was used for target in lfp, can be used for something
% else here. for eg. baseline vs inactivation
t = 1;

for cn = 1:length(ecg_bna_cfg.conditions)
    fprintf('Condition %s\n', ecg_bna_cfg.conditions(cn).label);
    sessions_avg(t).condition(cn).hs_tuned_evoked = struct();
    sessions_avg(t).condition(cn).cfg_condition = ecg_bna_cfg.conditions(cn);
    sessions_avg(t).condition(cn).label = ecg_bna_cfg.conditions(cn).label;
    
    for st = 1:size(ecg_bna_cfg.analyse_Rpeak_states, 1)
        for hs = 1:size(ecg_bna_cfg.conditions(cn).hs_labels, 2)
            sessions_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsessions = 0;
            sessions_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg_b2bt = {};
        end
    end
    
    for i = 1:length(ecg_b2bt_evoked.session)
        
        if isempty(ecg_b2bt_evoked.session(i).condition) || ...
                ~isfield(ecg_b2bt_evoked.session(i).condition(cn).hs_tuned_evoked, 'mean') || ...
                all(isnan([ecg_b2bt_evoked.session(i).condition(cn).hs_tuned_evoked.mean]));
            continue;
        end
        for st = 1:size(ecg_b2bt_evoked.session(i).condition(cn).hs_tuned_evoked, 1)
            for hs = 1:size(ecg_b2bt_evoked.session(i).condition(cn).hs_tuned_evoked, 2)
                if all(isnan(ecg_b2bt_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).mean))
                    continue;
                end
                sessions_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsessions = sessions_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsessions + 1;
                sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).time = ecg_b2bt_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).time;
                sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).hs_label = ecg_b2bt_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).hs_label;
                if isfield(ecg_b2bt_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs), 'state') && isfield(ecg_b2bt_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs), 'state_name')
                    sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).state = ecg_b2bt_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).state;
                    sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).state_name = ecg_b2bt_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).state_name;
                end
                sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt = [sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt, ecg_b2bt_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).mean];
               
            end
        end
    end
    
    % compute average
    if isfield(sessions_avg(t).condition(cn).hs_tuned_evoked, 'ecg_b2bt')
        for st = 1:size(sessions_avg(t).condition(cn).hs_tuned_evoked, 1)
            for hs = 1:size(sessions_avg(t).condition(cn).hs_tuned_evoked, 2)
                if ~isempty(sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt)
                    sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt = ...
                        cat(1, sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt{:});
                    sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).dimord = 'nsessions_time';
                    %                         sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt(...
                    %                             :, isnan(sum(sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt, 1))) = nan;
                    sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).std = ...
                        std(sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt, 0, 1);
                    sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).mean = ...
                        mean(sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt, 1);
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
            ecg_bna_plot_evoked_R2Rt(sessions_avg(t).condition(cn).hs_tuned_evoked, ...
                ecg_bna_cfg, plottitle, result_file);
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
            ecg_bna_plot_evoked_R2Rt(sessions_avg(t).difference(dcn).hs_tuned_evoked, ...
                ecg_bna_cfg, plottitle, result_file);
        end
    end
end

% save session average tfs
save(fullfile(results_fldr, 'sessions_evoked_ECG_R2Rt.mat'), 'sessions_avg');

close all;
end