function site_lfp= ecg_bna_process_LFP( sites, lfp_tfa_cfg,trials )

% lfp_tfa_process_LFP - function to read in the trial-wise LFP data for all
% sites recorded in a session, compute the LFP time frequency spectrogram,
% detect the noisy trials, and compute site-wise baseline power
%
% USAGE:
%	session_info = lfp_tfa_process_LFP( session_info, lfp_tfa_cfg )
%
% INPUTS:
%       session_info        - structure containing information about the
%       session to be processed, see settings/lfp_tfa_settings_example
%       Required fields:
%           Input               - path to the file which contains the
%                               trial-wise LFP data for all sites recorded
%                               in a session
%           proc_results_folder - folder where the results has to be
%           stored, see lfp_tfa_define_settings
%       lfp_tfa_cfg         - structure containing configurations for
%       reading LFP data and calculating spectrograms, see
%       settings/lfp_tfa_settings_example
%       Required fields:
%           ref_hemisphere          - reference hemisphere ('R'/'L') for ipsi-
%                                   and contra- hand and space labeling
%           trialinfo.start_state   - the ID of the state(event) which
%                                   marks the beginning of a trial, see
%                                   lfp_tfa_global_states
%           trialinfo.ref_tstart    - the offset from the onset of
%                                   trialinfo.start_state to be considered
%                                   as start time of a trial
%           trialinfo.end_state     - the ID of the state(event) which
%                                   marks the end of a trial, see
%                                   lfp_tfa_global_states
%           trialinfo.ref_tend      - the offset from the onset of
%                                   trialinfo.start_state to be considered
%                                   as end time of a trial
%           noise                   - settings for noisy trial detection,
%                                   see lfp_tfa_reject_noisy_lfp_trials for
%                                   more details
%           analyses                - which kind of analyses should be
%                                   performed on LFP, (Should be a
%                                   combination of 'tfs', 'evoked', 'pow',
%                                   'sync' and 'syncsp')
%
% OUTPUTS:
%		session_info            - same as input structure session_info
%
% REQUIRES: lfp_tfa_compute_site_tfr, lfp_tfa_reject_noisy_lfp_trials,
% lfp_tfa_compute_site_baseline, lfp_tfa_compute_sitepair_csd
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_define_settings,
% lfp_tfa_global_states, lfp_tfa_reject_noisy_lfp_trials,
% lfp_tfa_compute_site_baseline

% struct to save data for a site
site_lfp = struct();
fprintf('=============================================================\n');
fprintf('Processing site, %s\n', sites.site_ID);
site_lfp.session = sites.site_ID(1:12);
site_lfp.site_ID = sites.site_ID;
site_lfp.target = sites.target;
site_lfp.recorded_hemisphere = upper(sites.target(end));
site_lfp.xpos = sites.grid_x;
site_lfp.ypos = sites.grid_y;
site_lfp.zpos = sites.electrode_depth;

sitetrials=ph_get_unit_trials(sites,trials);
positions=[sitetrials.tar_pos]-[sitetrials.fix_pos];
hemifields=num2cell(sign(real(positions)));
fixations=num2cell([sitetrials.fix_pos]);
positions=num2cell(positions);
[sitetrials.position]=deal(positions{:});
[sitetrials.hemifield]=deal(hemifields{:});
[sitetrials.fixation]=deal(fixations{:});
sitetrials=ph_LR_to_CI(lfp_tfa_cfg,sites,sitetrials);  %% convert...

[sites.trial]=sitetrials;

lfp_samples=[0 cumsum(sites.LFP_samples)];

%% now loop through each trial for this site to do what exactly?
for t = 1:length(sites.trial)
    
    %% retrieve LFP data
    
    start_time = (sites.trial(t).TDT_LFPx_tStart); % trial start time
    fs = sites.trial(t).TDT_LFPx_SR; % sample rate
    %LFP = sites.trial(t).LFP; % LFP data
    LFP = sites.LFP((lfp_samples(t)+1):lfp_samples(t+1)); % LFP data
    ts = (1/fs); % sample time
    nsamples = numel(LFP);
    end_time = start_time + (ts*(nsamples-1));
    timestamps = linspace(start_time, end_time, nsamples);
    
    site_lfp.trials(t).lfp_data    = LFP;
    site_lfp.trials(t).time        = timestamps;
    site_lfp.trials(t).fsample     = fs;
    site_lfp.trials(t).tsample     = ts;
    
    %% need information about which trial it was(originally?)
    site_lfp.trials(t).n  = sites.trial(t).n;
    
    % save retrieved data into struct
    
    perturbation = sites.trial(t).perturbation; % 0 = control
    if isnan(perturbation)
        perturbation = 0;
    end
    site_lfp.trials(t).perturbation  = perturbation;
    site_lfp.trials(t).completed    = sites.trial(t).completed;
    site_lfp.trials(t).success      = sites.trial(t).success;
    site_lfp.trials(t).type         = sites.trial(t).type;
    site_lfp.trials(t).effector     = sites.trial(t).effector;
    site_lfp.trials(t).run          = sites.trial(t).run;
    site_lfp.trials(t).block        = sites.trial(t).block;
    site_lfp.trials(t).choice_trial = sites.trial(t).choice;
    
    % flag to mark noisy trials, default False, filled in by lfp_tfa_reject_noisy_lfp.m
    site_lfp.trials(t).noisy = 0;
    
    % get state onset times and onset samples - test and delete
    site_lfp.trials(t).states = struct();
    
    for s = 1:length(sites.trial(t).states)
        % get state ID
        state_id = sites.trial(t).states(s);
        % get state onset time
        state_onset = sites.trial(t).states_onset(sites.trial(t).states == state_id);
        % get sample number of state onset time
        state_onset_sample = find(abs(timestamps - state_onset) == min(abs(timestamps - state_onset)));
        % save into struct
        site_lfp.trials(t).states(s).id = state_id;
        site_lfp.trials(t).states(s).onset_t  = state_onset;
        site_lfp.trials(t).states(s).onset_s  = state_onset_sample;
    end
    
    trial_start_t = site_lfp.trials(t).states([site_lfp.trials(t).states.id] == lfp_tfa_cfg.trialinfo.start_state).onset_t + lfp_tfa_cfg.trialinfo.ref_tstart;
    trial_end_t   = site_lfp.trials(t).states([site_lfp.trials(t).states.id] == lfp_tfa_cfg.trialinfo.end_state).onset_t + lfp_tfa_cfg.trialinfo.ref_tend;
    site_lfp.trials(t).trialperiod = [trial_start_t, trial_end_t];
end

%% Time frequency spectrogram calculation
site_lfp = ecg_bna_compute_site_lfp_tfr( site_lfp, lfp_tfa_cfg );

% Noise rejection
site_lfp = lfp_tfa_reject_noisy_lfp_trials( site_lfp, lfp_tfa_cfg.noise );
end

