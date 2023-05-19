function [ session_evoked_ecg ] = ...
    ecg_bna_compute_session_state_evoked_ECG( session_ecg, session_info, analyse_states, ecg_bna_cfg )

% ecg_bna_compute_session_evoked_ECG  - compute average evoked ECG for
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
%       evoked ECG in a time window around Rpeak for trials of given
%       trial conditions
%
% REQUIRES:	lfp_tfa_compare_conditions, lfp_tfa_get_condition_trials,
% ecg_bna_get_Rpeak_based_STA, ecg_bna_get_shuffled_Rpeak_evoked_ECG,
% ecg_bna_plot_evoked_lfp, ecg_bna_compute_diff_condition_average
%
% See also ecg_bna_compute_session_Rpeak_evoked_LFP,
% ecg_bna_compute_session_Rpeak_evoked_TFS,
% ecg_bna_compute_session_evoked_ECG_R2Rt,
% ecg_bna_compute_session_Rpeak_evoked_state_onsets

% suppress warning for xticklabel
warning ('off', 'MATLAB:hg:willberemoved');

% make a folder to save figures
results_folder_evoked = fullfile(session_info.analyse_ecg_fldr, 'State_evoked_ECG');
if ~exist(results_folder_evoked, 'dir')
    mkdir(results_folder_evoked);
end

% condition based Evoked
sites_evoked = struct();
session_evoked_ecg = struct();
session_evoked_ecg.session = session_ecg(1).session;

% get trial conditions for this session
site_conditions = lfp_tfa_compare_conditions(ecg_bna_cfg, {0, 1});

% loop through each session (?? this always is only one session!)
% for i = 1:length(session_ecg)
i=1;
rng(ecg_bna_cfg.random_seed); % set random seed for reproducibility

% folder to save sitewise results
site_results_folder = fullfile(results_folder_evoked);
if ~exist(site_results_folder, 'dir')
    mkdir(site_results_folder);
end
% struct to store condition-wise evoked
sites_evoked(i).condition = struct();
sites_evoked(i).session = session_ecg(i).session;
% flag to indicate if this site should be used for
% averaging based on minimum no:of trials per condition
sites_evoked(i).use_for_avg = 1;


%% shuffle STATE ONSETS across all trials!
% if isfield(ecg_bna_cfg, 'random_permute_triggers') && ecg_bna_cfg.random_permute_triggers
%     nshuffles = 100;
%     if isfield(ecg_bna_cfg, 'n_shuffles')
%         nshuffles = ecg_bna_cfg.n_shuffles;
%     end
%     shuffled_Rpeaks = ecg_bna_get_shuffled_Rpeak(session_ecg(i).trials, nshuffles);
% end

% loop through conditions
for cn = 1:length(site_conditions)
    
    % hand-space tuning of LFP
    hs_labels = site_conditions(cn).hs_labels;
    
    % num sites
    nsites = length(session_ecg);
    
    % store details of analysed condition
    sites_evoked(i).condition(cn).label = site_conditions(cn).label;
    sites_evoked(i).condition(cn).cfg_condition = site_conditions(cn);
    sites_evoked(i).condition(cn).hs_tuned_evoked = struct();
    sites_evoked(i).condition(cn).ntrials = zeros(1,length(hs_labels));
    
    % loop through hand space labels
    for hs = 1:length(hs_labels)
        % get trial indices for the given condition
        cond_trials = lfp_tfa_get_condition_trials(session_ecg(i), site_conditions(cn));
        % filter trials by hand-space labels
        if ~strcmp(site_conditions(cn).reach_hands{hs}, 'any')
            cond_trials = cond_trials & strcmp({session_ecg(i).trials.reach_hand}, site_conditions(cn).reach_hands{hs});
        end
        if ~strcmp(site_conditions(cn).reach_spaces{hs}, 'any')
            cond_trials = cond_trials & strcmp({session_ecg(i).trials.reach_space}, site_conditions(cn).reach_spaces{hs});
        end
        
        sites_evoked(i).condition(cn).ntrials(hs) = sum(cond_trials);
        
        fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
        fprintf('Total number of trials %g\n', sum(cond_trials));
        
        cond_trials = cond_trials & ~[session_ecg(i).trials.noisy];
        sites_evoked(i).condition(cn).noisytrials(hs) = sum(cond_trials);
        
        % consider only non noisy trials
        fprintf('Number of noisy trials %g\n',  sum(cond_trials));
        
        % check if the site contains a specified minimum number
        % of trials for all conditions
        if sum(cond_trials) < ecg_bna_cfg.mintrials_percondition
            sites_evoked(i).use_for_avg = 0;
        end
        
       if any(cond_trials)
            cond_trials_ecg = session_ecg(i).trials(cond_trials);
        else
            continue;
       end
        
        
        % loop through time windows around the states to analyse
        for st = 1:size(analyse_states, 1)
                %cond_trials_ecg = ecg_bna_get_state_onset_indexes(cond_trials_ecg, analyse_states{st, 1});
                state_evoked = ecg_bna_get_state_evoked_ECG(cond_trials_ecg, session_ecg.session, analyse_states(st, :), 'ecg',analyse_states{st, 1});
                
%                 %shuffled_Rpeak_evoked.trial =zeros(size(state_evoked.mean));
%                 shuffled_Rpeak_evoked.mean=zeros(size(state_evoked.mean));
%                 shuffled_Rpeak_evoked.std =zeros(size(state_evoked.std));
%                 if isfield(ecg_bna_cfg, 'random_permute_triggers') && ecg_bna_cfg.random_permute_triggers
%                     for shuff=1:nshuffles
%                         cond_trials_shuffled=shuffled_Rpeaks(shuff,cond_trials);
%                         [cond_trials_ecg.ECG_spikes]=cond_trials_shuffled.ECG_spikes;
%                         shuffled_evoked(shuff) = ecg_bna_get_Rpeak_based_STA(cond_trials_ecg, session_ecg.session, analyse_states(st, :), 'ecg');
%                     end
%                     
%                     %shuffled_Rpeak_evoked.lfp = {shuffled_evoked.mean}; %% this aint right...
%                     shuffled_Rpeak_evoked.mean = nanmean(cat(1, shuffled_evoked.mean), 1);
%                     shuffled_Rpeak_evoked.std = nanstd(cat(1, shuffled_evoked.mean), 0, 1);
%                 end
            
            if ~isempty(state_evoked.trial)
                
                % save evoked ECG
                sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).trial = state_evoked.trial;
                sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean = state_evoked.mean;
                sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).std = state_evoked.std;
                sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time = state_evoked.time;
                                    
                %%%sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).shuffled_trial = shuffled_Rpeak_evoked.trial;
%                 sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).shuffled_mean  = shuffled_Rpeak_evoked.mean;
%                 sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).shuffled_std   = shuffled_Rpeak_evoked.std;
                
                sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).trial_idx = find(cond_trials);
                if strfind(state_evoked.dimord, 'npeaks')
                    sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).npeaks = size(state_evoked.trial, 1);
                elseif strfind(state_evoked.dimord, 'nshuffles')
                    sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).nshuffles = size(state_evoked.trial, 1);
                end
                sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).dimord = state_evoked.dimord;
                sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).hs_label = hs_labels(hs);
                if isfield(state_evoked, 'state_id') && isfield(state_evoked, 'state_name')
                    sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state = state_evoked.state_id;
                    sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state_name = state_evoked.state_name;
                end
            end
        end
    end
    
    % plots
    % Evoked ECG
    if ~isempty(fieldnames(sites_evoked(i).condition(cn).hs_tuned_evoked))
        plottitle = ['Session: ', sites_evoked(i).session ...
            site_conditions(cn).label '), '];
        if site_conditions(cn).choice == 0
            plottitle = [plottitle 'Instructed trials'];
        elseif site_conditions(cn).choice == 1
            plottitle = [plottitle 'Choice trials'];
        end
        result_file = fullfile(site_results_folder, ['State_evoked_ECG_' sites_evoked(i).session '_' site_conditions(cn).label]);        
        ecg_bna_plot_evoked_ECG (sites_evoked(i).condition(cn).hs_tuned_evoked, ecg_bna_cfg, plottitle, result_file, 'ylabel', 'ECG amplitude');
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
    sites_evoked(i).difference = [sites_evoked(i).difference, ecg_bna_compute_diff_condition_average('Rpeak_evoked_ECG', sites_evoked(i).condition, diff_condition, diff_color, diff_legend)];
end
% plot Difference evoked
for dcn = 1:length(sites_evoked(i).difference)
    if ~isempty(sites_evoked(i).difference(dcn).hs_tuned_evoked)
        if isfield(sites_evoked(i).difference(dcn).hs_tuned_evoked,'mean')
            plottitle = ['Session: ', sites_evoked(i).session sites_evoked(i).difference(dcn).label];
            result_file = fullfile(site_results_folder, ['ECG_DiffEvoked_' sites_evoked(i).session '_' 'diff_condition' num2str(dcn)]);
            %sites_avg(t).difference(dcn).label '.png']);
            ecg_bna_plot_evoked_lfp(sites_evoked(i).difference(dcn).hs_tuned_evoked, ecg_bna_cfg, plottitle, result_file, 'ylabel', 'ECG amplitude');
        end
    end
end

site_evoked_ecg = sites_evoked(i);
% save mat file for site
save(fullfile(site_results_folder, ['ECG_evoked_' sites_evoked(i).session '.mat']), 'site_evoked_ecg');
% save to a mother struct
session_evoked_ecg = site_evoked_ecg;
end

