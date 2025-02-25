function sites_avg = ecg_bna_avg_sites_Rpeak_evoked_LFP(Rpeak_evoked_LFP, ecg_bna_cfg)
%ecg_bna_avg_sites_Rpeak_evoked_LFP  - Condition-based evoked LFP response
% average across many sites from a single session or multiple sessions
%
% USAGE:
%	sites_avg = ecg_bna_avg_sites_Rpeak_evoked_LFP(Rpeak_evoked_LFP,
%	ecg_bna_cfg)
%
% INPUTS:
%		Rpeak_evoked_LFP	- struct containing the condition-based evoked
%		LFP response for indiviual sites, output of
%		ecg_bna_compute_session_Rpeak_evoked_LFP.m
%           Required Fields:
%               session.sites - 1xM struct containing condition-based
%               average evoked LFP response for M sites
%		ecg_bna_cfg         - struct containing the required settings
%           Required Fields:
%               conditions          - trial conditions to compare, see
%               lfp_tfa_settings.m and lfp_tfa_compare_conditions.m
%               root_results_fldr   - root folder where results are saved
%               compare.targets     - targets to compare, see lfp_tfa_settings.m
%               ref_hemisphere      - reference hemisphere for ipsi and
%               contra labeling
%               diff_condition      - conditions to compare, the plot
%               for compared conditions would be shown one on top of the
%               other
%           Optional Fields:
%               diff_color          - color to be used for plotting the
%               compared conditions
%               diff_legend         - legend to be used while plotting the
%               compared conditions
% OUTPUTS:
%		sites_avg           - structure containing condition-based evoked
%		LFP response averaged across multiple sites
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
results_fldr = fullfile(ecg_bna_cfg.analyse_lfp_folder, 'ECG_triggered_avg_across_sites');
if ~exist(results_fldr, 'dir')
    mkdir(results_fldr);
end

% Average Evoked LFP response across sessions
sites_avg = struct();

for t = 1:length(ecg_bna_cfg.compare.targets)
    sites_avg(t).target = ecg_bna_cfg.compare.targets{t};
    for cn = 1:length(ecg_bna_cfg.conditions)
        fprintf('Condition %s\n', ecg_bna_cfg.conditions(cn).label);
        sites_avg(t).condition(cn).hs_tuned_evoked = struct();
        sites_avg(t).condition(cn).cfg_condition = ecg_bna_cfg.conditions(cn);
        sites_avg(t).condition(cn).label = ecg_bna_cfg.conditions(cn).label;
        
        % initialize number of site pairs for each handspace label
        for st = 1:size(ecg_bna_cfg.analyse_states, 1)
            for hs = 1:size(ecg_bna_cfg.conditions(cn).hs_labels, 2)
                sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = 0;
                sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).trial = {};
                sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).shuffled_trial = {};
            end
        end
        
        for i = 1:length(Rpeak_evoked_LFP.session)
            for k = 1:length(Rpeak_evoked_LFP.session(i).sites)
                if isempty(Rpeak_evoked_LFP.session(i).sites(k).condition) || ...
                        ~strcmp(Rpeak_evoked_LFP.session(i).sites(k).target, ecg_bna_cfg.compare.targets{t}) || ...
                        ~isfield(Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked, 'mean') || ...
                        all(isnan([Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked.mean]));
                    continue;
                end
                for st = 1:size(Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked, 1)
                    for hs = 1:size(Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked, 2)
                        if all(isnan(Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked(st, hs).mean))
                            continue;
                        end
                        sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites + 1;
                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).time = Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked(st, hs).time;
                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).hs_label = Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked(st, hs).hs_label;
                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial = [sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial, Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked(st, hs).mean];
                        % subtract or plot extra?
                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).shuffled_trial = [sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).shuffled_trial, Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked(st, hs).shuffled_mean];
                    end
                end
            end
        end
        
        % compute average
        for st = 1:size(sites_avg(t).condition(cn).hs_tuned_evoked, 1)
            for hs = 1:size(sites_avg(t).condition(cn).hs_tuned_evoked, 2)
                sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial = cat(1, sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial{:});
                sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).shuffled_trial = cat(1, sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).shuffled_trial{:});
                sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).dimord = 'nsites_time';
                %                         sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt(...
                %                             :, isnan(sum(sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt, 1))) = nan;
                sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).std = std(sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial, 0, 1);
                sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).mean = mean(sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial, 1);
                % subtract or plot extra ?
                sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).shuffled_std = std(sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).shuffled_trial, 0, 1);
                sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).shuffled_mean = mean(sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).shuffled_trial, 1);
            end
        end
        
        plottitle = [ecg_bna_cfg.compare.targets{t},ecg_bna_cfg.conditions(cn).label];
        result_file = fullfile(results_fldr, ['Rpeak_evoked_LFP_' ecg_bna_cfg.compare.targets{t} '_' ecg_bna_cfg.conditions(cn).label]);
        ecg_bna_plot_evoked_lfp(sites_avg(t).condition(cn).hs_tuned_evoked, ecg_bna_cfg, plottitle, result_file, 'ylabel', 'ECG amplitude');
    end
    
    % difference between conditions
    sites_avg(t).difference = [];
    for diff = 1:size(ecg_bna_cfg.diff_condition, 2)
        diff_condition = ecg_bna_cfg.diff_condition{diff};
        diff_color = []; diff_legend = [];
        if isfield(ecg_bna_cfg, 'diff_color')
            diff_color = ecg_bna_cfg.diff_color{diff};
        end
        if isfield(ecg_bna_cfg, 'diff_legend')
            diff_legend = ecg_bna_cfg.diff_legend{diff};
        end
        sites_avg(t).difference = [sites_avg(t).difference, ...
            ecg_bna_compute_diff_condition_average('Event_evoked_ECG_R2Rt', sites_avg(t).condition, diff_condition, diff_color,diff_legend)];
    end
    % plot Difference TFR
    for dcn = 1:length(sites_avg(t).difference)
        if ~isempty(sites_avg(t).difference(dcn).hs_tuned_evoked)
            if isfield(sites_avg(t).difference(dcn).hs_tuned_evoked,'mean')
                plottitle = ['Target ', ecg_bna_cfg.compare.targets{t}, sites_avg(t).difference(dcn).label];
                result_file = fullfile(results_fldr, ['Rpeak_diffevoked_LFP_' 'diff_condition' num2str(dcn)]);
                %sessions_avg(t).difference(dcn).label '.png']);
                ecg_bna_plot_evoked_lfp(sites_avg(t).difference(dcn).hs_tuned_evoked, ecg_bna_cfg, plottitle, result_file, 'ylabel', 'ECG amplitude');
            end
        end
    end
end

% save session average tfs
save(fullfile(results_fldr, 'sites_avg_Rpeak_evoked_LFP.mat'), 'sites_avg');

close all;
end