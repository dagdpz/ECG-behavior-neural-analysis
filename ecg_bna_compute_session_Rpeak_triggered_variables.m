function [ site_Rpeak_triggered ] = ecg_bna_compute_session_Rpeak_triggered_variables( site_proc , session_ecg,analyse_states, ecg_bna_cfg )
% ecg_bna_compute_session_Rpeak_triggered_variables  - compute evoked_LFP,
% Powspctrm, ITPC, and phaseBP variables for different conditions for each site of a session.
% A condition is a combination of
% possibly hand-space tuning (for consitency with LFP analysis),
% control/inactivation, choice/instructed, type-effector values.
%
% USAGE:
%	[ session_data ] = ecg_bna_compute_session_Rpeak_triggered_variables(
%	site_proc, analyse_states, ecg_bna_cfg )
%
% INPUTS:
%		site_proc  	- 1xN struct containing raw LFP data for a
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
%		session_data	- output structure which saves the average
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



% folder to save sitewise results
% site_results_folder = fullfile(ecg_bna_cfg.sites_lfp_fldr, 'Rpeak_evoked_LFP');
site_results_folder = fullfile(ecg_bna_cfg.sites_lfp_fldr);
if ~exist(site_results_folder, 'dir')
    mkdir(site_results_folder);
end

% % condition based Evoked
% sites_data = struct();
% session_data = struct();
% session_data.session = site_proc(1).session;

% get trial conditions for this session
site_conditions = lfp_tfa_compare_conditions(ecg_bna_cfg, {0, 1});

% % num sites
% nsites = length(site_proc);
% 
% % loop through each site
% for i = 1:nsites
% struct to store condition-wise evoked
trig.condition = struct();
trig.site_ID = site_proc.site_ID;
trig.session = site_proc.session;
trig.target  = site_proc.target;

% loop through conditions
for cn = 1:length(site_conditions)
    % hand-space tuning of LFP
    hs_labels = site_conditions(cn).hs_labels;
    
    % store details of analysed condition
    trig.condition(cn).label    = site_conditions(cn).label;
    trig.condition(cn).cfg_condition = site_conditions(cn);
    trig.condition(cn).ntrials  = zeros(1,length(hs_labels));
    trig.condition(cn).state_hs = struct();
    
    % loop through hand space labels
    for hs = 1:length(hs_labels)
        % get trial indices for the given condition
        cond_trials = lfp_tfa_get_condition_trials(site_proc, site_conditions(cn));
        % filter trials by hand-space labels
        if ~strcmp(site_conditions(cn).reach_hands{hs}, 'any')
            cond_trials = cond_trials & strcmp({site_proc.trials.reach_hand}, site_conditions(cn).reach_hands{hs});
        end
        if ~strcmp(site_conditions(cn).reach_spaces{hs}, 'any')
            cond_trials = cond_trials & strcmp({site_proc.trials.reach_space}, site_conditions(cn).reach_spaces{hs});
        end
        
        trig.condition(cn).ntrials(hs)     = sum(cond_trials);
        trig.condition(cn).noisytrials(hs) = sum(cond_trials & [site_proc.trials.noisy]);
        fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
        fprintf('Total number of trials %g\n', sum(cond_trials));
        fprintf('Number of noisy trials %g\n', sum(cond_trials & [site_proc.trials.noisy]));
        cond_trials = cond_trials & ~[site_proc.trials.noisy];
        
        if sum(cond_trials) == 0
            continue;
        end
        cond_trials=find(cond_trials);
        trialsblocks_site=[site_proc.trials(cond_trials).block; site_proc.trials(cond_trials).n];
        trialsblocks_session=[session_ecg.trials.block; session_ecg.trials.n];
        trial_idx=ismember(trialsblocks_session',trialsblocks_site','rows');
        cix=ismember(trialsblocks_site',trialsblocks_session(:,trial_idx)','rows');
        if sum(cix) == 0
            continue;
        end
        cond_ecg = session_ecg.trials(trial_idx);
        fprintf('Condition Trials = %i , Num Trial_idx = %i, CIX = %i\n',numel(cond_trials),sum(trial_idx), sum(cix))
        
        % loop through time windows around the states to analyse
        for st = 1:size(analyse_states, 1)
            cond_LFP=site_proc;
            [cond_LFP.trials]=cond_LFP.trials(cond_trials(cix));
            [cond_LFP.trials.ECG_spikes] = cond_ecg.ECG_spikes;
            
            %if strcmp(analyse_states{st, 1}, 'ecg') % honestly, this needs to go inside the triggering functions
            %end
            real = ecg_bna_get_triggered_parameters(cond_LFP, analyse_states(st, :), ecg_bna_cfg);
            if isfield(ecg_bna_cfg, 'random_permute_triggers') && ecg_bna_cfg.random_permute_triggers
                %% compute shuffled power spectra, ITPC spectra, lfp, and bandpassed ITPC:
                [cond_LFP.trials.ECG_spikes]=cond_ecg.ECG_spikes_shuffled;
                shuffled = ecg_bna_get_triggered_parameters(cond_LFP, analyse_states(st, :), ecg_bna_cfg);
            else % some sort of dummies
                shuffled_evoked.lfp.mean=zeros(size(real.lfp.mean));
                shuffled_evoked.lfp.std =zeros(size(real.lfp.std));
                shuffled_evoked.itpcbp.mean=zeros(size(real.itpcbp.mean));
                shuffled_evoked.itpcbp.std =zeros(size(real.itpcbp.std));
                shuffled_evoked.powbp.mean=zeros(size(real.powbp.mean));
                shuffled_evoked.powbp.std =zeros(size(real.powbp.std));
            end
            if ~isempty(real.pow.mean) || ~isempty(real.lfp.mean)
                if isfield(real, 'state') && isfield(real, 'state_name')
                    trig.condition(cn).state_hs(st, hs).state = real.state;
                    trig.condition(cn).state_hs(st, hs).state_name = real.state_name;
                end
                trig.condition(cn).state_hs(st, hs).trials = cond_trials(cix); %find(cond_trials);
                trig.condition(cn).state_hs(st, hs).ntrials = sum(cix);
                trig.condition(cn).state_hs(st, hs).hs_label = hs_labels(hs);
            end
            if ~isempty(shuffled)
                normalized = ecg_bna_compute_shufflePredictor_normalization_general(real,shuffled,ecg_bna_cfg);
                significance = ecg_bna_compute_significance(real,shuffled,ecg_bna_cfg);
            end
            if ~isempty(real.pow.mean) && ~isempty(real.lfp.mean)
                trig.condition(cn).state_hs(st, hs).real=real;
                trig.condition(cn).state_hs(st, hs).shuffled=shuffled;
                trig.condition(cn).state_hs(st, hs).normalized=normalized;
                trig.condition(cn).state_hs(st, hs).significance=significance;
            end
            methods= {'real','shuffled','normalized'};
            for mt = 1: numel(methods)
                trig.condition(cn).state_hs(st, hs).(methods{mt}).time=real.time;
                trig.condition(cn).state_hs(st, hs).(methods{mt}).tfr_time=real.tfr_time;
                trig.condition(cn).state_hs(st, hs).(methods{mt}).state=real.state;
                trig.condition(cn).state_hs(st, hs).(methods{mt}).state_name=real.state_name;
                trig.condition(cn).state_hs(st, hs).(methods{mt}).freq=real.freq;
                %trig.condition(cn).state_hs(st, hs).(methods{mt}).hs_label=real_tfs.hs_label;
            end
        end
    end
end

% plots
methods= {'real','shuffled','normalized'};
for mt = 1: numel(methods)
    close all
    ecg_bna_plots_per_session( trig, site_conditions, ecg_bna_cfg, methods{mt}) % per site!
    % Notice: ===> last input could be 'real', 'shuffled', or 'normalized'
end
%
site_Rpeak_triggered = trig;
save(fullfile(site_results_folder, [trig.site_ID '.mat']), 'site_Rpeak_triggered');
close all;
end