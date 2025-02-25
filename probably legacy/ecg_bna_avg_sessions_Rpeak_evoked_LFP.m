function sessions_avg = ecg_bna_avg_sessions_Rpeak_evoked_LFP(Rpeak_evoked_LFP, ecg_bna_cfg)
%ecg_bna_avg_sessions_Rpeak_evoked_LFP  - Rpeak evoked LFP response
% average across many session averages
%
% USAGE:
%	sessions_avg = ecg_bna_avg_sessions_Rpeak_evoked_LFP(Rpeak_evoked_LFP,
%	ecg_bna_cfg)
%
% INPUTS:
%		Rpeak_evoked_LFP	- struct containing the average Rpeak evoked
%		LFP for individual sessions. See
%       ecg_bna_compute_session_Rpeak_evoked_LFP.m
%           Required Fields:
%               1. session.session_avg - 1xN struct containing
%               average evoked LFP response for N sessions (session_avg =
%               Average of site averages for one session)
%		ecg_bna_cfg     - struct containing the required settings
%           Required Fields:
%               conditions          - trial conditions to compare, see
%               lfp_tfa_settings.m and lfp_tfa_compare_conditions.m
%               root_results_fldr   - root folder where results are saved
%               compare.targets     - targets to compare, see lfp_tfa_settings.m
%               diff_condition      - conditions to compare, the plot
%               for compared conditions would be shown one on top of the
%               other
%           Optional Fields:
%               diff_color          - color to be used for plotting the
%               compared conditions
%               diff_legend         - legend to be used while plotting the
%               compared conditions
%
% OUTPUTS:
%		sessions_avg    - structure containing condition-based evoked LFP
%		response averaged across multiple sessions
%
% REQUIRES:	ecg_bna_plot_evoked_lfp
%
% See also ecg_bna_compute_session_Rpeak_evoked_LFP,
% ecg_bna_plot_evoked_lfp
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
results_fldr = fullfile(ecg_bna_cfg.analyse_lfp_folder, 'ECG_triggered_avg_across_sessions');
if ~exist(results_fldr, 'dir')
    mkdir(results_fldr);
end

% Average Evoked LFP response across sessions
sessions_avg = struct();

for t = 1:length(ecg_bna_cfg.compare.targets)
    sessions_avg(t).target = ecg_bna_cfg.compare.targets{t};
    for cn = 1:length(ecg_bna_cfg.conditions)
        fprintf('Condition %s\n', ecg_bna_cfg.conditions(cn).label);
        sessions_avg(t).condition(cn).hs_tuned_evoked = struct();
        sessions_avg(t).condition(cn).cfg_condition = ecg_bna_cfg.conditions(cn);
        sessions_avg(t).condition(cn).label = ecg_bna_cfg.conditions(cn).label;
        
        % initialize number of site pairs for each handspace label
        for st = 1:size(ecg_bna_cfg.analyse_states, 1)
            for hs = 1:size(ecg_bna_cfg.conditions(cn).hs_labels, 2)
                sessions_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsessions = 0;
                sessions_avg(t).condition(cn).hs_tuned_evoked(st, hs).trial = {};
                sessions_avg(t).condition(cn).hs_tuned_evoked(st, hs).trial_shuffled = {};
            end
        end
        
        for i = 1:length(Rpeak_evoked_LFP.session)
            for k = 1:length(Rpeak_evoked_LFP.session(i).session_avg)
                if isempty(Rpeak_evoked_LFP.session(i).session_avg(k).condition) || ...
                        ~strcmp(Rpeak_evoked_LFP.session(i).session_avg(k).target, ecg_bna_cfg.compare.targets{t}) || ...
                        ~isfield(Rpeak_evoked_LFP.session(i).session_avg(k).condition(cn).hs_tuned_evoked, 'mean') || ...
                        all(isnan([Rpeak_evoked_LFP.session(i).session_avg(k).condition(cn).hs_tuned_evoked.mean]));
                    continue;
                end
                for st = 1:size(Rpeak_evoked_LFP.session(i).session_avg(k).condition(cn).hs_tuned_evoked, 1)
                    for hs = 1:size(Rpeak_evoked_LFP.session(i).session_avg(k).condition(cn).hs_tuned_evoked, 2)
                        if all(isnan(Rpeak_evoked_LFP.session(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).mean))
                            continue;
                        end
                        sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).nsessions  = sessions_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsessions + 1;
                        sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).time       = Rpeak_evoked_LFP.session(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).time;
                        sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).hs_label   = Rpeak_evoked_LFP.session(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).hs_label;
                        sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial      = [sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial, Rpeak_evoked_LFP.session(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).mean];
                        
                        %% subtract or plot extra ?
                        sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial_shuffled = [sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial_shuffled, Rpeak_evoked_LFP.session(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).shuffled_mean];
                    end
                end
            end
        end
        
        % compute average
        for st = 1:size(sessions_avg(t).condition(cn).hs_tuned_evoked, 1)
            for hs = 1:size(sessions_avg(t).condition(cn).hs_tuned_evoked, 2)
                sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial = cat(1, sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial{:});
                sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).dimord = 'nsessions_time';
                %                         sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt(...
                %                             :, isnan(sum(sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt, 1))) = nan;
                sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).std  = std(sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial, 0, 1);
                sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).mean = mean(sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial, 1);
                
                %% remove if subtraction
                sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).shuffled_std  = std(vertcat(sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial_shuffled{:}), 0, 1);              
                sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).shuffled_mean = mean(vertcat(sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial_shuffled{:}), 1);
                
            end
        end
        
        plottitle = [ecg_bna_cfg.compare.targets{t},ecg_bna_cfg.conditions(cn).label];
        result_file = fullfile(results_fldr, ['Rpeak_evoked_LFP_' ecg_bna_cfg.compare.targets{t} '_' ecg_bna_cfg.conditions(cn).label]);
        ecg_bna_plot_evoked_lfp(sessions_avg(t).condition(cn).hs_tuned_evoked, ecg_bna_cfg, plottitle, result_file, 'ylabel', 'LFP amplitude');
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
        sessions_avg(t).difference = [sessions_avg(t).difference, ecg_bna_compute_diff_condition_average('Rpeak_evoked_LFP', sessions_avg(t).condition, diff_condition, diff_color,diff_legend)];
    end
    % plot Difference LFP
    for dcn = 1:length(sessions_avg(t).difference)
        if ~isempty(sessions_avg(t).difference(dcn).hs_tuned_evoked)
            if isfield(sessions_avg(t).difference(dcn).hs_tuned_evoked,'mean')
                plottitle = ['Target ', ecg_bna_cfg.compare.targets{t}, sessions_avg(t).difference(dcn).label];
                result_file = fullfile(results_fldr, ['Rpeak_diffevoked_LFP_' 'diff_condition' num2str(dcn)]);
                %sessions_avg(t).difference(dcn).label '.png']);
                ecg_bna_plot_evoked_lfp(sessions_avg(t).difference(dcn).hs_tuned_evoked, ecg_bna_cfg, plottitle, result_file, 'ylabel', 'LFP amplitude');
            end
        end
    end
end

% save session average tfs
save(fullfile(results_fldr, 'sessions_avg_evoked_LFP.mat'), 'sessions_avg');

close all;
end