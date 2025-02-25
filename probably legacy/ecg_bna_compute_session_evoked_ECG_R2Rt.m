function [ session_R2Rt ] = ecg_bna_compute_session_evoked_ECG_R2Rt( session_ecg, session_info, analyse_states, ecg_bna_cfg ) 

% ecg_bna_compute_session_evoked_ECG_R2Rt - compute average evoked ECG for
% different conditions for each session. A condition is a combination of
% possibly hand-space tuning (for consitency with LFP analysis),
% control/inactivation, choice/instructed, type-effector values. 
%
% USAGE:
%	[ session_evoked ] = ecg_bna_compute_session_evoked_ECG( session_ecg,
%	session_info, analyse_states, ecg_bna_cfg )  
%
% INPUTS:
%		session_ecg  	- struct containing raw ECG data for a session, 
%       output from ecg_bna_read_preproc_ECG or
%       ecg_bna_read_combined_ECG
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
%		session_evoked_ecg	- output structure which saves the average 
%       evoked ECG in a time window around for trials of given conditions
% 
% REQUIRES:	lfp_tfa_compare_conditions, lfp_tfa_get_condition_trials,
% ecg_bna_get_Rpeak_based_STA, ecg_bna_get_shuffled_Rpeak_evoked_ECG,
% ecg_bna_plot_evoked_lfp, ecg_bna_compute_diff_condition_average
%
% See also ecg_bna_compute_session_Rpeak_evoked_LFP,
% ecg_bna_compute_session_Rpeak_evoked_TFS,
% ecg_bna_compute_session_evoked_ECG,
% ecg_bna_compute_session_Rpeak_evoked_state_onsets 
    
    % suppress warning for xticklabel
    warning ('off', 'MATLAB:hg:willberemoved');

    % make a folder to save figures
    results_folder_ecg = fullfile(session_info.analyse_ecg_fldr, 'ECG R2Rt');
    if ~exist(results_folder_ecg, 'dir')
        mkdir(results_folder_ecg);
    end
       
    % condition based Evoked
    session_R2Rt = struct();
    
    % get trial conditions for this session
    site_conditions = lfp_tfa_compare_conditions(ecg_bna_cfg, {0, 1});
    
    % loop through each site
    for i = 1:length(session_ecg) 

        %rng(ecg_bna_cfg.random_seed); % set random seed for reproducibility
        
        % folder to save sitewise results
        session_results_folder = fullfile(results_folder_ecg);
        if ~exist(session_results_folder, 'dir')
            mkdir(session_results_folder);
        end
        % struct to store condition-wise evoked
        session_R2Rt(i).condition = struct();
        session_R2Rt(i).session = session_ecg(i).session;
        % flag to indicate if this site should be used for
        % averaging based on minimum no:of trials per condition
        session_R2Rt(i).use_for_avg = 1;
        
        % loop through conditions
        for cn = 1:length(site_conditions)

            % hand-space tuning of LFP
            hs_labels = site_conditions(cn).hs_labels;

            % num sites
            nsites = length(session_ecg);                 
            
            % store details of analysed condition
            session_R2Rt(i).condition(cn).label = site_conditions(cn).label;
            session_R2Rt(i).condition(cn).cfg_condition = site_conditions(cn);
            session_R2Rt(i).condition(cn).hs_tuned_evoked = struct(); 
            session_R2Rt(i).condition(cn).ntrials = zeros(1,length(hs_labels));        

            % loop through hand space labels
            for hs = 1:length(hs_labels)
                % get trial indices for the given condition
                cond_trials = lfp_tfa_get_condition_trials(session_ecg(i), site_conditions(cn));
                % filter trials by hand-space labels
                if ~strcmp(site_conditions(cn).reach_hands{hs}, 'any')
                    cond_trials = cond_trials & ...
                        strcmp({session_ecg(i).trials.reach_hand}, ...
                        site_conditions(cn).reach_hands{hs});
                end
                if ~strcmp(site_conditions(cn).reach_spaces{hs}, 'any')
                    cond_trials = cond_trials & ...
                        strcmp({session_ecg(i).trials.reach_space}, ...
                        site_conditions(cn).reach_spaces{hs});
                end
                
                session_R2Rt(i).condition(cn).ntrials(hs) = sum(cond_trials);

                fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
                fprintf('Total number of trials %g\n', sum(cond_trials));

                session_R2Rt(i).condition(cn).noisytrials(hs) = ...
                    sum(cond_trials & [session_ecg(i).trials.noisy]); 

                % consider only non noisy trials
                fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                    & [session_ecg(i).trials.noisy]));
                cond_trials = cond_trials & ~[session_ecg(i).trials.noisy];

                % check if the site contains a specified minimum number
                % of trials for all conditions
                if sum(cond_trials) < ecg_bna_cfg.mintrials_percondition
                    session_R2Rt(i).use_for_avg = 0;
                end


                % loop through time windows around the states to analyse
                for st = 1:size(analyse_states, 1)
                    
                    cond_trials_ecg = session_ecg(i).trials(cond_trials);                    
                     
                    state_evoked = ecg_bna_get_state_evoked_ECG_R2Rt(cond_trials_ecg, ...
                        analyse_states(st, :), ecg_bna_cfg);                    


                    if ~isempty(state_evoked.ecg_b2bt)

                        % save evoked ECG
                        session_R2Rt(i).condition(cn).hs_tuned_evoked(st, hs).ecg_b2bt = state_evoked.ecg_b2bt;
                        session_R2Rt(i).condition(cn).hs_tuned_evoked(st, hs).trials_mean = state_evoked.mean_ecg_b2bt;
                        session_R2Rt(i).condition(cn).hs_tuned_evoked(st, hs).mean = state_evoked.mean;
                        session_R2Rt(i).condition(cn).hs_tuned_evoked(st, hs).std = state_evoked.std; 
                        session_R2Rt(i).condition(cn).hs_tuned_evoked(st, hs).time = state_evoked.ecg_time;
                        session_R2Rt(i).condition(cn).hs_tuned_evoked(st, hs).trials = find(cond_trials);
                        session_R2Rt(i).condition(cn).hs_tuned_evoked(st, hs).valid_trials = ...
                            session_R2Rt(i).condition(cn).hs_tuned_evoked(st, hs).trials(state_evoked.valid_trials);
                        session_R2Rt(i).condition(cn).hs_tuned_evoked(st, hs).ntrials = state_evoked.ntrials;
                        session_R2Rt(i).condition(cn).hs_tuned_evoked(st, hs).hs_label = hs_labels(hs);
                        if isfield(state_evoked, 'state_id') && isfield(state_evoked, 'state_name')
                            session_R2Rt(i).condition(cn).hs_tuned_evoked(st, hs).state = state_evoked.state_id;
                            session_R2Rt(i).condition(cn).hs_tuned_evoked(st, hs).state_name = state_evoked.state_name;
                        end

                    end

                end

            end
            
            % plots
            % Evoked LFP
            if ~isempty(fieldnames(session_R2Rt(i).condition(cn).hs_tuned_evoked))
                plottitle = ['Session: ', session_R2Rt(i).session ...
                    site_conditions(cn).label ')'];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                elseif site_conditions(cn).choice == 1
                    plottitle = [plottitle 'Choice trials'];
                end
                result_file = fullfile(session_results_folder, ...
                    ['ECG_b2bt_Evoked_' session_R2Rt(i).session '_' site_conditions(cn).label]);
                   
                ecg_bna_plot_evoked_R2Rt (session_R2Rt(i).condition(cn).hs_tuned_evoked, ecg_bna_cfg, ...
                    plottitle, result_file);
            end

        end
        
        % difference between conditions
        session_R2Rt(i).difference = [];
        for diff = 1:size(ecg_bna_cfg.diff_condition, 2)
            diff_condition = ecg_bna_cfg.diff_condition{diff};
            diff_color = []; diff_legend = [];
            if isfield(ecg_bna_cfg, 'diff_color')
                diff_color = ecg_bna_cfg.diff_color{diff};
            end
            if isfield(ecg_bna_cfg, 'diff_legend')
                diff_legend = ecg_bna_cfg.diff_legend{diff};
            end
            session_R2Rt(i).difference = [session_R2Rt(i).difference, ...
                ecg_bna_compute_diff_condition_average('Event_evoked_ECG_R2Rt', ...
                session_R2Rt(i).condition, diff_condition, diff_color, diff_legend)];
        end
        % plot Difference TFR
        for dcn = 1:length(session_R2Rt(i).difference)
            if ~isempty(session_R2Rt(i).difference(dcn).hs_tuned_evoked)
                if isfield(session_R2Rt(i).difference(dcn).hs_tuned_evoked,... 
                        'mean')
                    plottitle = ['Session: ', session_R2Rt(i).session ...
                        session_R2Rt(i).difference(dcn).label];
                    result_file = fullfile(session_results_folder, ...
                        ['ECG_b2bt_DiffEvoked_' session_R2Rt(i).session ...
                        '_' 'diff_condition' num2str(dcn)]);
                        %sites_avg(t).difference(dcn).label '.png']);
                    ecg_bna_plot_evoked_R2Rt(session_R2Rt(i).difference(dcn).hs_tuned_evoked, ...
                        ecg_bna_cfg, plottitle, result_file);
                end
            end
        end
        
        session_evoked_ecg_R2Rt = session_R2Rt(i);
        % save mat file for site
        save(fullfile(session_results_folder, ['ECG_R2Rt_evoked_' ...
            session_R2Rt(i).session '.mat']), 'session_evoked_ecg_R2Rt');

    end
        
%     % Average across sites for a session
%     session_avg = struct();
%     % targets for this session
%     targets = unique({session_ecg.target});
%     % average each target separately
%     for t = 1:length(targets)
%         session_avg(t).target = targets{t};
%         
%         % conditions
%         for cn = 1:length(site_conditions)
%             session_avg(t).condition(cn).hs_tuned_evoked = struct();
%             isite = 0;
%             session_avg(t).condition(cn).cfg_condition = lfp_tfa_cfg.conditions(cn);
%             session_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
%             for st = 1:size(sites_evoked(1).condition(1).hs_tuned_evoked, 1)
%                 for hs = 1:size(sites_evoked(1).condition(1).hs_tuned_evoked, 2)
%                     session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = 0;
%                     session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg = [];
%                     session_avg(t).condition(cn).hs_tuned_evoked(st, hs).mean = [];
%                     session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std = [];
%                     session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time = [];
%                 end
%             end
%             for i = 1:length(sites_evoked)
%                 % if the site's target is same as target being considered
%                 if ~strcmp(session_ecg(i).target, targets{t})
%                     continue;
%                 end
%                 if sites_evoked(i).use_for_avg
%                     % calculate the average evoked ECG across sites for this condition 
%                     if ~isempty(sites_evoked(i).condition(cn).hs_tuned_evoked) && ...
%                         isfield(sites_evoked(i).condition(cn).hs_tuned_evoked, 'mean')
%                         isite = isite + 1;
%                         for hs = 1:size(sites_evoked(i).condition(cn).hs_tuned_evoked, 2)
%                             for st = 1:size(sites_evoked(i).condition(cn).hs_tuned_evoked, 1)
%                                 if ~isempty(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean)
%                                     session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = ...
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites + 1;
%                                     if session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites == 1
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg = ...
%                                             [session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg; ...
%                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean] ;
% %                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std = ...
% %                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).std ;
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time = ...
%                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time;
% 
%                                     else
%                                         nsamples = length(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time);
%                                         if nsamples > length(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time)
%                                             nsamples = length(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time);
%                                         end
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg = ...
%                                             [session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg(:,1:nsamples); ...
%                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean(1:nsamples)] ;
% %                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std = ...
% %                                             session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std(1:nsamples) + ...
% %                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).std(1:nsamples) ;
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time = ...
%                                             session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time(1:nsamples) ;
% 
%                                     end
%                                     % struct to store average evoked LFP across sites
%                                     session_avg(t).condition(cn).hs_tuned_evoked(st, hs).hs_label = ...
%                                         sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).hs_label;
%                                     if isfield(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs), 'state') && ...
%                                             isfield(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs), 'state_name')
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).state = ...
%                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state;
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).state_name = ...
%                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state_name;
%                                     end
%                                     session_avg(t).condition(cn).condition = site_conditions(cn);
%                                     session_avg(t).condition(cn).label = site_conditions(cn).label;
%                                     session_avg(t).condition(cn).session = session_ecg(i).session;
%                                     session_avg(t).condition(cn).target = session_ecg(i).target;
%                                     %session_avg(t).condition(cn).nsites = nsites;
%                                 end
% 
%                             end
%                         end
%                     end
%                 end
%             end
%             % average TFR across sites for a session
%             if isfield(session_avg(t).condition(cn).hs_tuned_evoked, 'ecg') && ...
%                     ~isempty([session_avg(t).condition(cn).hs_tuned_evoked.ecg])
%                 
%                 for hs = 1:size(session_avg(t).condition(cn).hs_tuned_evoked, 2)
%                     for st = 1:size(session_avg(t).condition(cn).hs_tuned_evoked, 1)
%                         if ~isempty(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg)
%                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).mean = ...
%                             nanmean(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg, 1);
%                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std = ...
%                             nanstd(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg, 0, 1);
%                         %session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = isite;
%                         end
%                     end
%                 end
%             end           
%             % plot average evoked ECG across sites for this session
%             if ~isempty([session_avg(t).condition(cn).hs_tuned_evoked.ecg])
%                 plottitle = ['Session: ', session_avg(t).condition(cn).session ', Target = ' ...
%                     session_avg(t).condition(cn).target ' (ref_' lfp_tfa_cfg.ref_hemisphere '), '  ...
%                     site_conditions(cn).label ', '];
%                 if site_conditions(cn).choice == 0
%                     plottitle = [plottitle 'Instructed trials'];
%                 elseif site_conditions(cn).choice == 1
%                     plottitle = [plottitle 'Choice trials'];
%                 end
%                 result_file = fullfile(results_folder_evoked, ['ECG_Evoked_' ...
%                     session_avg(t).condition(cn).session '_' ...
%                     session_avg(t).condition(cn).target '_' site_conditions(cn).label '.png']);
%                 lfp_tfa_plot_evoked_lfp (session_avg(t).condition(cn).hs_tuned_evoked, lfp_tfa_cfg, ...
%                     plottitle, result_file);
%             end
%         end 
%         
%         
%         % difference between conditions
%         session_avg(t).difference = [];
%         for diff = 1:size(lfp_tfa_cfg.diff_condition, 2)
%             diff_condition = lfp_tfa_cfg.diff_condition{diff};
%             session_avg(t).difference = [session_avg(t).difference, ...
%                 lfp_tfa_compute_diff_condition_evoked(session_avg(t).condition, diff_condition)];
%         end
%         % plot Difference TFR
%         for dcn = 1:length(session_avg(t).difference)
%             if ~isempty(session_avg(t).difference(dcn).hs_tuned_evoked)
%                 if isfield(session_avg(t).difference(dcn).hs_tuned_evoked,... 
%                         'mean')
%                     plottitle = ['Target ', lfp_tfa_cfg.compare.targets{t}, ...
%                     ' (ref_', lfp_tfa_cfg.ref_hemisphere, ') ', ...
%                     session_avg(t).difference(dcn).label];
%                     result_file = fullfile(results_fldr, ...
%                         ['ECG_DiffEvoked_' lfp_tfa_cfg.compare.targets{t} ...
%                         '_' 'diff_condition' num2str(dcn) '.png']);
%                         %session_avg(t).difference(dcn).label '.png']);
%                     lfp_tfa_plot_evoked_lfp(session_avg(t).difference(dcn).hs_tuned_evoked, ...
%                         lfp_tfa_cfg, plottitle, result_file);
%                 end
%             end
%         end
% 
%     end
%         
%     close all;
%     
%     % store session average data
%     session_evoked_ecg.session_avg = session_avg;
%     
%     % save mat files
%     save(fullfile(results_folder_evoked, ['ECG_evoked_' session_evoked_ecg.session '.mat']), 'session_evoked');
%     % save settings file
%     save(fullfile(results_folder_evoked, 'lfp_tfa_settings.mat'), 'lfp_tfa_cfg');
end
        
