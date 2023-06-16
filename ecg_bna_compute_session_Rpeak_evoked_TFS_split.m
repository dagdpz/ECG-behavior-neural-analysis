function [session_tfs] = ecg_bna_compute_session_Rpeak_evoked_TFS_split( session_proc_lfp,session_ecg, analyse_states, ecg_bna_cfg )
% ecg_bna_compute_session_Rpeak_evoked_LFP  - compute average LFP time
% frequency spectrograms for specified time windows around Rpeak from
% trials of  given different conditions for each site of a session.
% A condition is a combination of possibly hand-space tuning (for
% consitency with LFP analysis), control/inactivation, choice/instructed,
% type-effector values.
%
% USAGE:
%	[ session_tfs ] = ecg_bna_compute_session_Rpeak_evoked_TFS(
%	session_proc_lfp, analyse_states, ecg_bna_cfg )
%
% INPUTS:
%		session_proc_lfp  	- 1xN struct containing raw LFP data for a
%		session,  output from ecg_bna_process_combined_LFP_ECG
%       analyse_states  - cell array containing states to be
%       analysed and corresponding time windows
%       ecg_bna_cfg     - struct containing configuration settings
%           Required fields:
%               session_results_fldr        - folder to which the
%               results of the session should be saved
%               mintrials_percondition          - minimum number of trials
%               required per condition for considering the site for
%               averaging
%               diff_condition      - conditions to compare, the plot
%               for compared conditions would be shown one on top of the
%               other
%               ref_hemisphere      - reference hemisphere for ipsi- and
%               contra- hand and space labeling
%
% OUTPUTS:
%		session_tfs	- output structure which saves the average
%       LFP TFS in a time window around Rpeak for trials of given
%       conditions
%
% REQUIRES:	lfp_tfa_compare_conditions, lfp_tfa_get_condition_trials,
% ecg_bna_get_ECG_triggered_tfs, lfp_tfa_compute_difference_condition_tfr,
% lfp_tfa_plot_hs_tuned_tfr_multiple_img
%
% See also ecg_bna_compute_session_evoked_ECG,
% ecg_bna_compute_session_Rpeak_evoked_LFP,
% ecg_bna_compute_session_evoked_ECG_R2Rt,
% ecg_bna_compute_session_Rpeak_evoked_state_onsets


% suppress warning for xticklabel
warning ('off', 'MATLAB:hg:willberemoved');

% make a folder to save figures
results_folder_tfr = fullfile(ecg_bna_cfg.session_lfp_fldr, 'Rpeak_evoked_TFS');
if ~exist(results_folder_tfr, 'dir')
    mkdir(results_folder_tfr);
end

% condition based TFS
% struct to store TFR for each site
sites_tfr = struct();
% struct to accumulate TFR for each site and store session average
session_tfs = struct();
% session name
session_tfs.session = session_proc_lfp(1).session;
% get the trial conditions for this session
site_conditions = lfp_tfa_compare_conditions(ecg_bna_cfg, {0, 1});

% loop through each site
for i = 1:length(session_proc_lfp)
    
    %rng(ecg_bna_cfg.random_seed); % set random seed for reproducibility
    
    site_lfp = session_proc_lfp(i);
    
    % folder to save sitewise results
    site_results_folder = fullfile(ecg_bna_cfg.sites_lfp_fldr, 'Rpeak_evoked_TFS');
    if ~exist(site_results_folder, 'dir')
        mkdir(site_results_folder);
    end
    
    % structure to store condition-wise tfs
    sites_tfr(i).condition = struct();
    % info about session and site
    sites_tfr(i).site_ID = site_lfp.site_ID;
    sites_tfr(i).session = site_lfp.session;
    sites_tfr(i).target = site_lfp.target;
    % flag to indicate if this site should be used for
    % averaging based on minimum no:of trials per condition
    sites_tfr(i).use_for_avg = 1;
    
    % loop through each condition
    for cn = 1:length(site_conditions)
        
        % hand-space tuning of LFP
        hs_labels = site_conditions(cn).hs_labels;
        
        % store details of condition analysed
        sites_tfr(i).condition(cn).label = site_conditions(cn).label;
        sites_tfr(i).condition(cn).cfg_condition = site_conditions(cn);
        sites_tfr(i).condition(cn).ntrials = zeros(1,length(hs_labels));
        sites_tfr(i).condition(cn).hs_tuned_tfs = struct();
        
        % loop through hand space labels
        for hs = 1:length(hs_labels)
            % get the trial indices which satisfy the given condition
            cond_trials = lfp_tfa_get_condition_trials(site_lfp, site_conditions(cn));
            
            % get the trial indices which satisfy the given hand-space
            % label for the given condition
            if ~strcmp(site_conditions(cn).reach_hands{hs}, 'any')
                cond_trials = cond_trials & ...
                    strcmp({session_proc_lfp(i).trials.reach_hand}, ...
                    site_conditions(cn).reach_hands{hs});
            end
            if ~strcmp(site_conditions(cn).reach_spaces{hs}, 'any')
                cond_trials = cond_trials & ...
                    strcmp({session_proc_lfp(i).trials.reach_space}, ...
                    site_conditions(cn).reach_spaces{hs});
            end
            
            sites_tfr(i).condition(cn).ntrials(hs) = sum(cond_trials);
            
            fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
            fprintf('Total number of trials %g\n', sum(cond_trials));
            
            sites_tfr(i).condition(cn).noisytrials(hs) = ...
                sum(cond_trials & [site_lfp.trials.noisy]);
            
            % consider only non noisy trials
            fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                & [session_proc_lfp(i).trials.noisy]));
            cond_trials = cond_trials & ~[site_lfp.trials.noisy];
            
            % check if the site contains a specified minimum number
            % of trials for all conditions
            if sum(cond_trials) < ecg_bna_cfg.mintrials_percondition
                sites_tfr(i).use_for_avg = 0;
            end
            % loop through states to analyse
            
            cond_trials=find(cond_trials);
            trialsblocks_site=[session_proc_lfp(i).trials(cond_trials).block; session_proc_lfp(i).trials(cond_trials).n];
            trialsblocks_session=[session_ecg.trials.block; session_ecg.trials.n];
            if sum(cond_trials) == 0
                trial_idx=[];
                cix=[];
            else
                trial_idx=ismember(trialsblocks_session',trialsblocks_site','rows');
                cix=ismember(trialsblocks_site',trialsblocks_session','rows');
            end
            
            
            cond_ecg = session_ecg.trials(trial_idx);
            
            for st = 1:size(analyse_states, 1)
                
                if strcmp(analyse_states{st, 1}, 'ecg')
                    cond_LFP=site_lfp;
                    [cond_LFP.trials]=site_lfp.trials(cond_trials(cix));
                    state_tfs = ecg_bna_get_ECG_triggered_tfs_split(cond_LFP, cond_ecg, analyse_states(st, :), ecg_bna_cfg);
                end
                
                if ~isempty(state_tfs.powspctrm)
                    
                    % save average tfs for this condition, hand-space
                    % label, and state
                    sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = state_tfs.powspctrm_normmean;
                    sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm_raw = state_tfs.powspctrm;
                    sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.time = state_tfs.time;
                    sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.freq = state_tfs.freq;
                    sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.cfg = state_tfs.cfg;
                    sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).hs_label = hs_labels(hs);
                    if isfield(state_tfs, 'state') && isfield(state_tfs, 'state_name')
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).state = state_tfs.state;
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).state_name = state_tfs.state_name;
                    end
                    %                         sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).trials = find(cond_trials);
                    %                         sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).ntrials = length(find(cond_trials));
                    sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).trials = cond_trials(cix);
                    sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).ntrials = sum(cix);
                end
                
            end
            
        end
        
        
        % plot TFR
        if ~isempty(fieldnames(sites_tfr(i).condition(cn).hs_tuned_tfs))
            if site_conditions(cn).perturbation == 0
                injection = 'Pre';
            elseif site_conditions(cn).perturbation == 1
                injection = 'Post';
            else
                injection = 'Any';
            end
            plottitle = ['LFP TFR (' injection '): Site ' sites_tfr(i).site_ID ', Target ' sites_tfr(i).target '(ref_' ecg_bna_cfg.ref_hemisphere '), '  site_conditions(cn).label];
            if site_conditions(cn).choice == 0
                plottitle = [plottitle 'Instructed trials'];
            elseif site_conditions(cn).choice == 1
                plottitle = [plottitle 'Choice trials'];
            end
            
            result_file = fullfile(site_results_folder, ['LFP_TFR_' sites_tfr(i).site_ID '_' site_conditions(cn).label '.png']);
            lfp_tfa_plot_hs_tuned_tfr_multiple_img(sites_tfr(i).condition(cn).hs_tuned_tfs, ecg_bna_cfg, plottitle, result_file);
        end
        
    end
    
    sites_tfr(i).difference = [];
    % difference between conditions
    for diff = 1:size(ecg_bna_cfg.diff_condition, 2)
        diff_condition = ecg_bna_cfg.diff_condition{diff};
        sites_tfr(i).difference = [sites_tfr(i).difference, ...
            lfp_tfa_compute_difference_condition_tfr(sites_tfr(i).condition, diff_condition)];
    end
    % Plot TFR difference
    for dcn = 1:length(sites_tfr(i).difference)
        if ~isempty(fieldnames(sites_tfr(i).difference(dcn).hs_tuned_tfs))
            plottitle = [' Target ' sites_tfr(i).target, ' (ref_', ecg_bna_cfg.ref_hemisphere, '), Site ', sites_tfr(i).site_ID sites_tfr(i).difference(dcn).label ];
            result_file = fullfile(site_results_folder, ['LFP_DiffTFR_' sites_tfr(i).site_ID '_' 'diff_condition' num2str(dcn) '.png']);
            lfp_tfa_plot_hs_tuned_tfr_multiple_img(sites_tfr(i).difference(dcn).hs_tuned_tfs, ecg_bna_cfg, plottitle, result_file, 'bluewhitered');
        end
    end
    %end
    close all;
    % save mat file for each site
    site_tfr = sites_tfr(i);
    save(fullfile(site_results_folder, ...
        ['LFP_TFR_' site_tfr.site_ID '.mat']), 'site_tfr');
    % save into a mother struct
    session_tfs.sites(i) = site_tfr;
    
end


% Calculate average TFR across all sites
session_avg = struct();
% targets for this session
targets = unique({session_proc_lfp.target});
% average each target separately
for t = 1:length(targets)
    session_avg(t).target = targets{t};
    session_avg(t).session = session_proc_lfp(1).session;
    % loop through conditions
    for cn = 1:length(site_conditions)
        % condition-wise session average tfs
        session_avg(t).condition(cn).hs_tuned_tfs = struct();
        session_avg(t).condition(cn).cfg_condition = site_conditions(cn);
        session_avg(t).condition(cn).label = site_conditions(cn).label;
        isite = zeros(size(ecg_bna_cfg.analyse_states, 1),size(ecg_bna_cfg.conditions(cn).hs_labels, 2));
        for i = 1:length(session_proc_lfp)
            site_lfp = session_proc_lfp(i);
            if ~strcmp(site_lfp.target, targets{t}) || ...
                    ~sites_tfr(i).use_for_avg || ...
                    isempty(sites_tfr(i).condition(cn).hs_tuned_tfs) && ...
                    ~isfield(sites_tfr(i).condition(cn).hs_tuned_tfs, 'freq')
                continue;
            end
            % check if this site should be used for averaging
            % calculate the average across sites for this condition
            
            for hs = 1:size(sites_tfr(i).condition(cn).hs_tuned_tfs, 2)
                for st = 1:size(sites_tfr(i).condition(cn).hs_tuned_tfs, 1)
                    if  ~isfield(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs),'freq') || ...
                            ~isfield(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq,'powspctrm') || ...
                            isempty(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm)
                        %%? make struct?
                        continue;
                    end
                    isite(st,hs) = isite(st,hs) + 1;
                    if isite(st,hs) == 1
                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
                            sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm ;
                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.time = ...
                            sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.time;
                    else
                        ntimebins = size(session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm, 3);
                        % average same number of time bins
                        if ntimebins > length(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.time)
                            ntimebins = length(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.time);
                            session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
                                session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm(1,:,1:ntimebins) + ...
                                sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm(1,:,1:ntimebins) ;
                        else
                            if ~isempty(session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm)
                                session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
                                    session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm + ...
                                    sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm(1,:,1:ntimebins) ;
                            else
                                session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
                                    sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm(1,:,1:ntimebins) ;
                            end
                        end
                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.time = ...
                            sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.time(1:ntimebins);
                    end
                    % store session tfs
                    session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.freq = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.freq;
                    session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.cfg  = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.cfg;
                    session_avg(t).condition(cn).hs_tuned_tfs(st, hs).hs_label  = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).hs_label;
                    if isfield(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs), 'state') && isfield(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs), 'state_name')
                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).state = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).state;
                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).state_name = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).state_name;
                    end
                    session_avg(t).condition(cn).cfg_condition = site_conditions(cn);
                    session_avg(t).condition(cn).label = site_conditions(cn).label;
                    session_avg(t).condition(cn).session = site_lfp.session;
                    session_avg(t).condition(cn).target = site_lfp.target;
                end
            end
            
        end
        % average TFR across sites for a session
        if isfield(session_avg(t).condition(cn).hs_tuned_tfs, 'freq')
            for hs = 1:size(session_avg(t).condition(cn).hs_tuned_tfs, 2)
                for st = 1:size(session_avg(t).condition(cn).hs_tuned_tfs, 1)
                    if isfield(session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq, 'powspctrm') %% field not present ?
                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
                            session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm / isite(st,hs);
                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsites = isite(st,hs);
                    end
                end
            end
        end
        
        % plot average TFR for this condition and target
        if ~isempty(session_avg(t).condition(cn).hs_tuned_tfs)
            if isfield(session_avg(t).condition(cn).hs_tuned_tfs, 'freq')
                plottitle = ['LFP TFR: Target = ' session_avg(t).target ', (ref_', ecg_bna_cfg.ref_hemisphere, ') ',  'Session ', site_lfp.session, ' ', site_conditions(cn).label];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                elseif site_conditions(cn).choice == 1
                    plottitle = [plottitle 'Choice trials'];
                end
                result_file = fullfile(results_folder_tfr, ...
                    ['LFP_TFR_' session_avg(t).target '_' site_lfp.session '_' site_conditions(cn).label '.png']);
                lfp_tfa_plot_hs_tuned_tfr_multiple_img(session_avg(t).condition(cn).hs_tuned_tfs, ecg_bna_cfg, plottitle, result_file);
            end
        end
        
    end
    
    % Difference TFR for session
    % check if both pre- and post- injection blocks exist
    session_avg(t).difference = [];
    % difference between conditions
    for diff = 1:size(ecg_bna_cfg.diff_condition, 2)
        diff_condition = ecg_bna_cfg.diff_condition{diff};
        session_avg(t).difference = [session_avg(t).difference, lfp_tfa_compute_difference_condition_tfr(session_avg(t).condition, diff_condition)];
    end
    
    % plot average TFR difference across sites for this session
    for dcn = 1:length(session_avg(t).difference)
        if ~isempty(fieldnames(session_avg(t).difference(dcn).hs_tuned_tfs))
            plottitle = ['LFP Diff TFR: Target ' session_avg(t).target '(ref_' ecg_bna_cfg.ref_hemisphere '), Session ', session_avg(t).difference(dcn).session, ' ', session_avg(t).difference(dcn).label];
            result_file = fullfile(results_folder_tfr, ['LFP_DiffTFR_' session_avg(t).target '_' session_avg(t).difference(dcn).session '_' 'diff_condition' num2str(dcn) '.png']);
            lfp_tfa_plot_hs_tuned_tfr_multiple_img(session_avg(t).difference(dcn).hs_tuned_tfs, ecg_bna_cfg, plottitle, result_file, 'bluewhitered');
        end
    end
    %end
    
end

session_tfs.session_avg = session_avg;

% close figures
close all;

% save session average tfs
save(fullfile(results_folder_tfr, ['LFP_TFR_' session_tfs.session '.mat']), 'session_tfs');

end
