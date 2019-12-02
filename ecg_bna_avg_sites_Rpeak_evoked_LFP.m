function sites_avg = ecg_bna_avg_sites_Rpeak_evoked_LFP(Rpeak_evoked_LFP, ecg_bna_cfg)
%lfp_tfa_avg_evoked_LFP_across_sessions  - Condition-based evoked LFP response
% average across many session averages
%
% USAGE:
%	sessions_avg = lfp_tfa_avg_sessions_ECG_evoked(evoked_ecg, lfp_tfa_cfg)
%
% INPUTS:
%		lfp_evoked		- struct containing the condition-based evoked LFP response for
%		indiviual sites, output of lfp_tfa_plot_site_evoked_LFP.m
%           Required Fields:
%               1. session.session_avg - 1xN struct containing condition-based
%               average evoked LFP response for N sessions (session_avg =
%               Average of site averages for one session)
%		lfp_tfa_cfg     - struct containing the required settings
%           Required Fields:
%               1. conditions          - trial conditions to compare, see
%               lfp_tfa_settings.m and lfp_tfa_compare_conditions.m
%               2. root_results_fldr   - root folder where results are saved
%               3. compare.targets     - targets to compare, see lfp_tfa_settings.m
%               4. ref_hemisphere      - reference hemisphere for ipsi and
%               contra labeling
% OUTPUTS:
%		sessions_avg    - structure containing condition-based evoked LFP
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
    results_fldr = fullfile(ecg_bna_cfg.root_results_fldr, 'ECG analysis');
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

            % initialize number of site pairs for each handspace
            % label
            for st = 1:size(Rpeak_evoked_LFP.session(1).sites(1).condition(cn).hs_tuned_evoked, 1)
                for hs = 1:size(Rpeak_evoked_LFP.session(1).sites(1).condition(cn).hs_tuned_evoked, 2)
                    sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = 0;
                    sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).trial = {};                
                end
            end  

            for i = 1:length(Rpeak_evoked_LFP.session)
                for k = 1:length(Rpeak_evoked_LFP.session(i).sites)
                    if isempty(Rpeak_evoked_LFP.session(i).sites(k).condition)
                        continue;
                    end
                    if ~strcmp(Rpeak_evoked_LFP.session(i).sites(k).target, ecg_bna_cfg.compare.targets{t})
                        continue;
                    end
                    if isfield(Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked, 'mean')
                        for st = 1:size(Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked, 1)
                            for hs = 1:size(Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked, 2)
                                if isfield(Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked(st, hs), 'mean') ...
                                        && ~isempty(Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked(st, hs).mean)
                                    sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = ...
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites + 1;
                                    if sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites == 1
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).time ...
                                            = Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked(st, hs).time;
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).hs_label ...
                                            = Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked(st, hs).hs_label;
                                        if isfield(Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked(st, hs), 'state') && ...
                                                isfield(Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked(st, hs), 'state_name')
                                            sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).state ...
                                                = Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked(st, hs).state;
                                            sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).state_name ...
                                                = Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked(st, hs).state_name;
                                        end
                                    end
                                    sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial ...
                                        = [sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial, ...
                                        Rpeak_evoked_LFP.session(i).sites(k).condition(cn).hs_tuned_evoked(st, hs).mean];  
                                else
                                    continue; %sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs) = struct();
                                end
                            end
                        end
                    end     
                end
            end

            % compute average
            if isfield(sites_avg(t).condition(cn).hs_tuned_evoked, 'trial')
                for st = 1:size(sites_avg(t).condition(cn).hs_tuned_evoked, 1)
                    for hs = 1:size(sites_avg(t).condition(cn).hs_tuned_evoked, 2)
                        if ~isempty(sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial)
                            sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial = ...
                                cat(1, sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial{:});
                            sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).dimord = 'nsites_time';
    %                         sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt(...
    %                             :, isnan(sum(sessions_avg(t).condition(cn).hs_tuned_evoked(st,hs).ecg_b2bt, 1))) = nan;
                            sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).std = ...
                                std(sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial, 0, 1);
                            sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).mean = ...
                                mean(sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).trial, 1);                           
                        end
                    end
                end
            end


            if ~isempty(sites_avg(t).condition(cn).hs_tuned_evoked)
                if isfield(sites_avg(t).condition(cn).hs_tuned_evoked,... 
                        'mean')
                    plottitle = [ecg_bna_cfg.compare.targets{t},...
                         ecg_bna_cfg.conditions(cn).label];
                    result_file = fullfile(results_fldr, ...
                                    ['ECG_b2bt_Evoked_' ecg_bna_cfg.conditions(cn).label]);
                    ecg_bna_plot_evoked_lfp(sites_avg(t).condition(cn).hs_tuned_evoked, ...
                                ecg_bna_cfg, plottitle, result_file, 'ylabel', 'ECG amplitude');
                end
            end

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
                ecg_bna_compute_diff_condition_average('Event_evoked_ECG_R2Rt', ...
                sites_avg(t).condition, diff_condition, diff_color,diff_legend)];
        end
        % plot Difference TFR
        for dcn = 1:length(sites_avg(t).difference)
            if ~isempty(sites_avg(t).difference(dcn).hs_tuned_evoked)
                if isfield(sites_avg(t).difference(dcn).hs_tuned_evoked,... 
                        'mean')
                    plottitle = ['Target ', ecg_bna_cfg.compare.targets{t}, ...
                        sites_avg(t).difference(dcn).label];
                    result_file = fullfile(results_fldr, ...
                        ['ECG_b2bt_DiffEvoked_' 'diff_condition' num2str(dcn)]);
                        %sessions_avg(t).difference(dcn).label '.png']);
                    ecg_bna_plot_evoked_lfp(sites_avg(t).difference(dcn).hs_tuned_evoked, ...
                        ecg_bna_cfg, plottitle, result_file, 'ylabel', 'ECG amplitude');
                end
            end
        end
    end
    
    % save session average tfs
    save(fullfile(results_fldr, 'sites_Rpeak_evoked_LFP.mat'), 'sites_avg');
        
    close all;
end