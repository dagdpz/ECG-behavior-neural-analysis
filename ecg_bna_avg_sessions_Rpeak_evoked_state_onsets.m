function sessions_avg = ecg_bna_avg_sessions_Rpeak_evoked_state_onsets(Rpeak_state_onset, ecg_bna_cfg)
%ecg_bna_avg_sessions_Rpeak_evoked_state_onsets  - Rpeak evoked event onset
%probability average across sessions
%
% USAGE:
%	sessions_avg =
%	ecg_bna_avg_sessions_Rpeak_evoked_state_onsets(Rpeak_state_onset, lfp_tfa_cfg)
%
% INPUTS:
%		Rpeak_state_onset	- struct containing the Rpeak evoked event onset probability for
%		indiviual sessions, output of
%		ecg_bna_compute_session_Rpeak_evoked_state_onsets.m
%           Required Fields:
%               1. session.session_avg - 1xN struct containing condition-based
%               Rpeak evoked event onset probability for N sessions
%		ecg_bna_cfg         - struct containing the required settings
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
% OUTPUTS:
%		sessions_avg    - structure containing condition-based Rpeak evoked
%		state onset probability averaged across multiple sessions
%
% REQUIRES:	ecg_bna_plot_Rpeak_ref_state_onsets
%
% See also ecg_bna_compute_session_Rpeak_evoked_state_onsets,
% ecg_bna_plot_Rpeak_ref_state_onsets
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
results_fldr = fullfile(ecg_bna_cfg.analyse_ecg_folder, ...
    'Rpeak_evoked_states');
if ~exist(results_fldr, 'dir')
    mkdir(results_fldr);
end

% Average Evoked LFP response across sessions
sessions_avg = struct();
t = 1;

for cn = 1:length(ecg_bna_cfg.conditions)
    fprintf('Condition %s\n', ecg_bna_cfg.conditions(cn).label);
    sessions_avg(t).condition(cn).Rpeak_evoked = struct();
    sessions_avg(t).condition(cn).cfg_condition = ecg_bna_cfg.conditions(cn);
    sessions_avg(t).condition(cn).label = ecg_bna_cfg.conditions(cn).label;
    
    % initialize number of site pairs for each handspacelabel
    for st = 1:size(ecg_bna_cfg.analyse_Rpeak_states, 1)
        for hs = 1:size(ecg_bna_cfg.conditions(cn).hs_labels, 2)
            sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).nsessions = 0;
            sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).abs_timefromRpeak = {};
            sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).rel_timefromRpeak = {};
            sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).abs_timeprob.prob = [];
            sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).rel_timeprob.prob = [];
            sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).abs_timeprob.timebins = [];
            sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).rel_timeprob.timebins = [];
        end
    end
    
    for i = 1:length(Rpeak_state_onset.session)
        if isempty(Rpeak_state_onset.session(i).condition) || ...
                ~isfield(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked, 'abs_timefromRpeak') || ...
                ~isfield(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked, 'rel_timefromRpeak')
            continue;
        end
        for st = 1:size(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked, 1)
            for hs = 1:size(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked, 2)
                if ~isfield(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs), 'abs_timefromRpeak') ...
                        || isempty(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).abs_timefromRpeak) ...
                        || ~isfield(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs), 'rel_timefromRpeak') ...
                        || isempty(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).rel_timefromRpeak)
                    continue
                end
                sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).nsessions = sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).nsessions + 1;
                %if sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).nsessions == 1
                sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).hs_label = Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).hs_label;
                if isfield(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs), 'state') && isfield(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs), 'state_name')
                    sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).state = Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).state;
                    sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).state_name = Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).state_name;
                end
                %end
                sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timefromRpeak = [sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timefromRpeak, Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).abs_timefromRpeak];
                sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).rel_timefromRpeak = [sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).rel_timefromRpeak, Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).rel_timefromRpeak];
                % accumulate histogram counts
                sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.timebins = Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).abs_timeprob.timebins;
                sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.prob = [sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.prob; Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).abs_timeprob.prob];
                sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).rel_timeprob.timebins = Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).rel_timeprob.timebins;
                sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).rel_timeprob.prob = [sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).rel_timeprob.prob; Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).rel_timeprob.prob];
                
            end
        end
        
    end
    
    % concatenate probability histograms
    %         for st = 1:size(sessions_avg(t).condition(cn).Rpeak_evoked, 1)
    %             for hs = 1:size(sessions_avg(t).condition(cn).Rpeak_evoked, 2)
    %                 if isfield(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs), 'abs_timeprob')
    %                     nsessions = length(sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.timebins);
    %                     [max_ntimebins, idx] = max(cellfun('size', ...
    %                         sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.timebins, 2));
    %                     sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.timebins = ...
    %                         sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.timebins{idx};
    %                     state_probability = nan*ones(nsessions, max_ntimebins-1);
    %                     for i = 1:nsessions
    %                         nprob = length(sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.prob{i});
    %                         state_probability(i, 1:nprob) = ...
    %                             sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.prob{i};
    %                     end
    %                     sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.prob = ...
    %                         state_probability;
    %                 end
    %             end
    %         end
    
    if ~isempty(sessions_avg(t).condition(cn).Rpeak_evoked)
        if isfield(sessions_avg(t).condition(cn).Rpeak_evoked,'abs_timefromRpeak') ...
                && isfield(sessions_avg(t).condition(cn).Rpeak_evoked, 'rel_timefromRpeak') ...
                && ~isempty([sessions_avg(t).condition(cn).Rpeak_evoked.abs_timefromRpeak]) ...
                && ~isempty([sessions_avg(t).condition(cn).Rpeak_evoked.rel_timefromRpeak])
            plottitle = [ecg_bna_cfg.compare.targets{t},ecg_bna_cfg.conditions(cn).label];
            result_file = fullfile(results_fldr, ['Rpeak_trig_P_event_' ecg_bna_cfg.conditions(cn).label]);
            ecg_bna_plot_Rpeak_ref_state_onsets (sessions_avg(t).condition(cn).Rpeak_evoked, ecg_bna_cfg, plottitle, result_file);
        end
    end
    % save session average tfs
    %save(fullfile(results_fldr, 'sessions_Rpeak_evoked_onsets.mat'), 'sessions_avg');
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
        ecg_bna_compute_diff_condition_average('Rpeak_evoked_states_onset', ...
        sessions_avg(t).condition, diff_condition, diff_color, diff_legend)];
end
% plot Difference TFR
for dcn = 1:length(sessions_avg(t).difference)
    if ~isempty(sessions_avg(t).difference(dcn).Rpeak_evoked)
        if isfield(sessions_avg(t).difference(dcn).Rpeak_evoked,'abs_timefromRpeak') ...
                && isfield(sessions_avg(t).difference(dcn).Rpeak_evoked, 'rel_timefromRpeak') ...
                && ~isempty([sessions_avg(t).difference(dcn).Rpeak_evoked.abs_timefromRpeak]) ...
                && ~isempty([sessions_avg(t).difference(dcn).Rpeak_evoked.rel_timefromRpeak])
            plottitle = ['Target ', ecg_bna_cfg.compare.targets{t}, sessions_avg(t).difference(dcn).label];
            result_file = fullfile(results_fldr, ['Diff_Rpeak_trig_P_event_' 'diff_condition' num2str(dcn)]);
            %sessions_avg(t).difference(dcn).label '.png']);
            ecg_bna_plot_Rpeak_ref_state_onsets(sessions_avg(t).difference(dcn).Rpeak_evoked, ecg_bna_cfg, plottitle, result_file);
        end
    end
end
% save session average tfs
save(fullfile(results_fldr, 'sessions_avg_Rpeak_evoked_onsets.mat'), 'sessions_avg');
close all;
end