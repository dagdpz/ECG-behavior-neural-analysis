function [ session_evoked ] = ecg_bna_compute_session_Rpeak_evoked_LFP( session_proc_lfp, analyse_states, ecg_bna_cfg ) 
% ecg_bna_compute_session_Rpeak_evoked_LFP  - compute average Rpeak evoked 
% ECG for different conditions for each site of a session. 
% A condition is a combination of
% possibly hand-space tuning (for consitency with LFP analysis),
% control/inactivation, choice/instructed, type-effector values. 
%
% USAGE:
%	[ session_evoked ] = ecg_bna_compute_session_Rpeak_evoked_LFP(
%	session_proc_lfp, analyse_states, ecg_bna_cfg )  
%
% INPUTS:
%		session_proc_lfp  	- 1xN struct containing raw LFP data for a
%		session,  output from ecg_bna_process_combined_LFP_ECG 
%       analyse_states  - cell array containing states to be
%       analysed and corresponding time windows
%       ecg_bna_cfg     - struct containing configuration settings
%           Required fields:
%               random_seed                 - random seed for
%               reproducibility of random shuffling of Rpeaks
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
%           Optional Fields:
%               diff_color          - color to be used for plotting the
%               compared conditions
%               diff_legend         - legend to be used while plotting the
%               compared conditions
%               random_permute_triggers     - flag indicating whether to
%               randomly permute the Rpeak triggers
%               n_shuffles                  - integer which specifies how
%               many times the Rpeaks need to be randomly shuffled to
%               compute statistics
%               
%
% OUTPUTS:
%		session_evoked	- output structure which saves the average 
%       evoked LFP in a time window around Rpeak for trials of given 
%       conditions
% 
% REQUIRES:	lfp_tfa_compare_conditions, lfp_tfa_get_condition_trials,
% ecg_bna_get_Rpeak_evoked_LFP, ecg_bna_get_shuffled_Rpeak_evoked_LFP,
% ecg_bna_plot_evoked_lfp, ecg_bna_compute_diff_condition_average
%
% See also ecg_bna_compute_session_evoked_ECG,
% ecg_bna_compute_session_Rpeak_evoked_TFS,
% ecg_bna_compute_session_evoked_ECG_R2Rt,
% ecg_bna_compute_session_Rpeak_evoked_state_onsets 

    % suppress warning for xticklabel
    warning ('off', 'MATLAB:hg:willberemoved');

    % make a folder to save figures
    results_folder_evoked = fullfile(ecg_bna_cfg.session_lfp_fldr, 'Rpeak_evoked_LFP');
    if ~exist(results_folder_evoked, 'dir')
        mkdir(results_folder_evoked);
    end
       
    % condition based Evoked
    sites_evoked = struct();
    session_evoked = struct();
    session_evoked.session = session_proc_lfp(1).session;

    % get trial conditions for this session
    site_conditions = lfp_tfa_compare_conditions(ecg_bna_cfg, {0, 1});
    
    % num sites
    nsites = length(session_proc_lfp);
    
    % loop through each site
    for i = 1:nsites
        rng(ecg_bna_cfg.random_seed); % set random seed for reproducibility
        
        % folder to save sitewise results
        site_results_folder = fullfile(ecg_bna_cfg.sites_lfp_fldr, 'Rpeak_evoked_LFP');
        if ~exist(site_results_folder, 'dir')
            mkdir(site_results_folder);
        end
        % struct to store condition-wise evoked
        sites_evoked(i).condition = struct();
        sites_evoked(i).site_ID = session_proc_lfp(i).site_ID;
        sites_evoked(i).session = session_proc_lfp(i).session;
        sites_evoked(i).target = session_proc_lfp(i).target;
        % flag to indicate if this site should be used for
        % averaging based on minimum no:of trials per condition
        sites_evoked(i).use_for_avg = 1;
        
        % loop through conditions
        for cn = 1:length(site_conditions)
            condition_valid=false;
            % hand-space tuning of LFP
            hs_labels = site_conditions(cn).hs_labels;
             
            % store details of analysed condition
            sites_evoked(i).condition(cn).label = site_conditions(cn).label;
            sites_evoked(i).condition(cn).cfg_condition = site_conditions(cn);
            sites_evoked(i).condition(cn).ntrials = zeros(1,length(hs_labels));
            sites_evoked(i).condition(cn).hs_tuned_evoked = struct();         

            % loop through hand space labels
            for hs = 1:length(hs_labels)
                % get trial indices for the given condition
                cond_trials = lfp_tfa_get_condition_trials(session_proc_lfp(i), site_conditions(cn));
                % filter trials by hand-space labels
                if ~strcmp(site_conditions(cn).reach_hands{hs}, 'any')
                    cond_trials = cond_trials & strcmp({session_proc_lfp(i).trials.reach_hand}, site_conditions(cn).reach_hands{hs});
                end
                if ~strcmp(site_conditions(cn).reach_spaces{hs}, 'any')
                    cond_trials = cond_trials & strcmp({session_proc_lfp(i).trials.reach_space}, site_conditions(cn).reach_spaces{hs});
                end
                
                sites_evoked(i).condition(cn).ntrials(hs)     = sum(cond_trials);
                sites_evoked(i).condition(cn).noisytrials(hs) = sum(cond_trials & [session_proc_lfp(i).trials.noisy]); 
                fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
                fprintf('Total number of trials %g\n', sum(cond_trials));
                fprintf('Number of noisy trials %g\n', sum(cond_trials & [session_proc_lfp(i).trials.noisy]));
                cond_trials = cond_trials & ~[session_proc_lfp(i).trials.noisy];

                % check if the site contains a specified minimum number
                % of trials for all conditions
                if sum(cond_trials) < ecg_bna_cfg.mintrials_percondition
                    sites_evoked(i).use_for_avg = 0;
                end
                
                if sum(cond_trials) == 0
                    continue;
                else
                    condition_valid=true;
                end


                % loop through time windows around the states to analyse
                for st = 1:size(analyse_states, 1)
                    cond_trials_lfp = session_proc_lfp(i).trials(cond_trials);
                    if strcmp(analyse_states{st, 1}, 'ecg')
                        state_evoked = ecg_bna_get_Rpeak_evoked_LFP(cond_trials_lfp, analyse_states(st, :));
                        if isfield(ecg_bna_cfg, 'random_permute_triggers') && ecg_bna_cfg.random_permute_triggers
                            nshuffles = 100;
                            if isfield(ecg_bna_cfg, 'n_shuffles')
                                nshuffles = ecg_bna_cfg.n_shuffles;
                            end
                            shuffled_Rpeak_evoked = ecg_bna_get_shuffled_Rpeak_evoked_LFP(cond_trials_lfp, analyse_states(st, :), nshuffles);
                        end
                    
                    end                        


                    if ~isempty(state_evoked.lfp)

                        % save evoked LFP
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).lfp = state_evoked.lfp;
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean = state_evoked.mean;
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).std = state_evoked.std;                          
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time = state_evoked.lfp_time;
                        
                        if isfield(ecg_bna_cfg, 'random_permute_triggers') && ecg_bna_cfg.random_permute_triggers
                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).shuffled_lfp = shuffled_Rpeak_evoked.lfp;
                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).shuffled_mean = shuffled_Rpeak_evoked.mean;
                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).shuffled_std = shuffled_Rpeak_evoked.std;
                        end
                        
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).trials = find(cond_trials);
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).hs_label = hs_labels(hs);
                        if isfield(state_evoked, 'state_id') && isfield(state_evoked, 'state_name')
                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state = state_evoked.state_id;
                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state_name = state_evoked.state_name;
                        end

                    end

                end

            end
            
            % plots
            % Evoked LFP
            if condition_valid && ~isempty([sites_evoked(i).condition(cn).hs_tuned_evoked.lfp])
                plottitle = ['Site ID: ', sites_evoked(i).site_ID ', Target = ' sites_evoked(i).target '(ref_' ecg_bna_cfg.ref_hemisphere '), '  site_conditions(cn).label '), '];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                elseif site_conditions(cn).choice == 1
                    plottitle = [plottitle 'Choice trials'];
                end
                result_file = fullfile(site_results_folder, ['LFP_Evoked_' sites_evoked(i).site_ID '_' site_conditions(cn).label '.png']);
                ecg_bna_plot_evoked_lfp (sites_evoked(i).condition(cn).hs_tuned_evoked, ecg_bna_cfg, plottitle, result_file);
            end

        end
        
        % difference between conditions
        sites_evoked(i).difference = [];
        for diff = 1:size(ecg_bna_cfg.diff_condition, 2)
            diff_condition = ecg_bna_cfg.diff_condition{diff};
            diff_color = []; diff_legend = [];
            if isfield(ecg_bna_cfg, 'diff_color')
                diff_color = ecg_bna_cfg.diff_color{diff};
            end
            if isfield(ecg_bna_cfg, 'diff_legend')
                diff_legend = ecg_bna_cfg.diff_legend{diff};
            end
            sites_evoked(i).difference = [sites_evoked(i).difference, ecg_bna_compute_diff_condition_average('Rpeak_evoked_LFP', ...
                sites_evoked(i).condition, diff_condition, diff_color, diff_legend)];
        end
        % plot Difference TFR
        for dcn = 1:length(sites_evoked(i).difference)
            if ~isempty(sites_evoked(i).difference(dcn).hs_tuned_evoked)
                if isfield(sites_evoked(i).difference(dcn).hs_tuned_evoked,'mean')
                    plottitle = ['Site ID: ' sites_evoked(i).site_ID ', Target = ' sites_evoked(i).target '(ref_' ecg_bna_cfg.ref_hemisphere '), ' sites_evoked(i).difference(dcn).label];
                    result_file = fullfile(site_results_folder, ['LFP_DiffEvoked_' sites_evoked(i).site_ID '_' 'diff_condition' num2str(dcn) '.png']);
                    ecg_bna_plot_evoked_lfp(sites_evoked(i).difference(dcn).hs_tuned_evoked, ecg_bna_cfg, plottitle, result_file);
                end
            end
        end
        
        site_evoked_lfp = sites_evoked(i);
        % save mat file for site
        save(fullfile(site_results_folder, ['Rpeak_evoked_LFP_' sites_evoked(i).site_ID '.mat']), 'site_evoked_lfp');
        % save to a mother struct
        session_evoked.sites(i) = site_evoked_lfp;
        
        close all;
    end
        
    % Average across sites for a session
    session_avg = struct();
    % targets for this session
    targets = unique({session_proc_lfp.target});
    % average each target separately
    for t = 1:length(targets)
        session_avg(t).target = targets{t};
        
        % conditions
        for cn = 1:length(site_conditions)
            session_avg(t).condition(cn).hs_tuned_evoked = struct();
            isite = 0;
            session_avg(t).condition(cn).cfg_condition = site_conditions(cn);
            session_avg(t).condition(cn).label = site_conditions(cn).label;
            session_avg(t).condition(cn).condition = site_conditions(cn);
            session_avg(t).condition(cn).label = site_conditions(cn).label;
            session_avg(t).condition(cn).session = session_proc_lfp(1).session;
            session_avg(t).condition(cn).target = session_proc_lfp(1).target;
            
            for st = 1:size(analyse_states, 1)
                for hs = 1:length(hs_labels)
                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = 0;
                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp = [];
                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).mean = [];
                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std = [];
                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time = [];
                end
            end
            for i = 1:length(sites_evoked)
                % if the site's target is same as target being considered
                if ~strcmp(session_proc_lfp(i).target, targets{t})
                    continue;
                end
                if sites_evoked(i).use_for_avg
                    % calculate the average evoked LFP across sites for this condition 
                    if ~isempty(sites_evoked(i).condition(cn).hs_tuned_evoked) && ...
                        isfield(sites_evoked(i).condition(cn).hs_tuned_evoked, 'mean')
                        isite = isite + 1;
                        for hs = 1:size(sites_evoked(i).condition(cn).hs_tuned_evoked, 2)
                            for st = 1:size(sites_evoked(i).condition(cn).hs_tuned_evoked, 1)
                                if ~isempty(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean)
                                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = ...
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites + 1;
                                    if session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites == 1
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp = ...
                                            [session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp; ...
                                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean] ;
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time = ...
                                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time;

                                    else % nsamples differs (hopefully only slightly (?) across sites
                                        nsamples = length(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time);
                                        if nsamples > length(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time)
                                            nsamples = length(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time);
                                        end
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp = ...
                                            [session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp(:,1:nsamples); ...
                                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean(1:nsamples)] ;
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time = ...
                                            session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time(1:nsamples) ;

                                    end
                                    % struct to store average evoked LFP across sites
                                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).hs_label = ...
                                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).hs_label;
                                    if isfield(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs), 'state') && ...
                                            isfield(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs), 'state_name')
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).state = ...
                                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state;
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).state_name = ...
                                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state_name;
                                    end
                                end

                            end
                        end
                    end
                end
            end
            % average TFR across sites for a session                
            for hs = 1:size(session_avg(t).condition(cn).hs_tuned_evoked, 2)
                for st = 1:size(session_avg(t).condition(cn).hs_tuned_evoked, 1)
                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).mean = ...
                        nanmean(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp, 1);
                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std = ...
                        nanstd(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp, 0, 1);
                end
            end
            % plot average evoked LFP across sites for this session
            plottitle = ['Session: ', session_avg(t).condition(cn).session ', Target = ' ...
                session_avg(t).condition(cn).target ' (ref_' ecg_bna_cfg.ref_hemisphere '), '  ...
                site_conditions(cn).label ', '];
            if site_conditions(cn).choice == 0
                plottitle = [plottitle 'Instructed trials'];
            elseif site_conditions(cn).choice == 1
                plottitle = [plottitle 'Choice trials'];
            end
            result_file = fullfile(results_folder_evoked, ['LFP_Evoked_' ...
                session_avg(t).condition(cn).session '_' ...
                session_avg(t).condition(cn).target '_' site_conditions(cn).label '.png']);
            ecg_bna_plot_evoked_lfp (session_avg(t).condition(cn).hs_tuned_evoked, ecg_bna_cfg, ...
                plottitle, result_file);
        end 
        
        
        % difference between conditions
        session_avg(t).difference = [];
        for diff = 1:size(ecg_bna_cfg.diff_condition, 2)
            diff_condition = ecg_bna_cfg.diff_condition{diff};
            diff_color = []; diff_legend = [];
            if isfield(ecg_bna_cfg, 'diff_color')
                diff_color = ecg_bna_cfg.diff_color{diff};
            end
            if isfield(ecg_bna_cfg, 'diff_legend')
                diff_legend = ecg_bna_cfg.diff_legend{diff};
            end
            session_avg(t).difference = [session_avg(t).difference, ...
                ecg_bna_compute_diff_condition_average('Rpeak_evoked_LFP', session_avg(t).condition, diff_condition, diff_color, diff_legend)];
        end
        % plot Difference TFR
        for dcn = 1:length(session_avg(t).difference)
            if ~isempty(session_avg(t).difference(dcn).hs_tuned_evoked)
                if isfield(session_avg(t).difference(dcn).hs_tuned_evoked,'mean')
                    plottitle = ['Target ', ecg_bna_cfg.compare.targets{t}, ' (ref_', ecg_bna_cfg.ref_hemisphere, ') ', session_avg(t).difference(dcn).label];
                    result_file = fullfile(results_folder_evoked, ['LFP_DiffEvoked_' ecg_bna_cfg.compare.targets{t} '_' 'diff_condition' num2str(dcn) '.png']);
                        %session_avg(t).difference(dcn).label '.png']);
                    ecg_bna_plot_evoked_lfp(session_avg(t).difference(dcn).hs_tuned_evoked, ecg_bna_cfg, plottitle, result_file);
                end
            end
        end

    end
        
    close all;
    
    % store session average data
    session_evoked.session_avg = session_avg;
    
    % save mat files
    save(fullfile(results_folder_evoked, ['Rpeak_evoked_LFP_' session_evoked.session '.mat']), 'session_evoked');

end
        
