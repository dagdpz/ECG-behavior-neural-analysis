function [session_info, allsites_lfp]= ecg_bna_process_LFP( session_info, lfp_tfa_cfg )

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

close all;

load(session_info.Input_LFP{:}, 'sites');

% prepare results folder
results_fldr = fullfile(session_info.proc_results_fldr);
if ~exist(results_fldr, 'dir')
    mkdir(results_fldr);
end


% structure array to store lfp data for all sites
% to be used for cross power spectrum calculation
allsites_lfp = [];

for i = 1:length(sites)
    
    % struct to save data for a site
    site_lfp = struct();
    fprintf('=============================================================\n');
    fprintf('Processing site, %s\n', sites(i).site_ID);
    site_lfp.session = sites(i).site_ID(1:12);
    site_lfp.site_ID = sites(i).site_ID;
    site_lfp.target = sites(i).target;
    site_lfp.recorded_hemisphere = upper(sites(i).target(end));
    site_lfp.ref_hemisphere = lfp_tfa_cfg.ref_hemisphere;
    site_lfp.xpos = sites(i).grid_x;
    site_lfp.ypos = sites(i).grid_y;
    site_lfp.zpos = sites(i).electrode_depth;
    
    %% LS 2021: convert to ipsi/contra instead of fixed reference hemisphere
    %% remove non-completed trials later !
    sitetrials=[sites(i).trial];
    positions=[sitetrials.tar_pos]-[sitetrials.fix_pos];
    hemifields=num2cell(sign(real(positions)));
    fixations=num2cell([sitetrials.fix_pos]);
    positions=num2cell(positions);
    [sitetrials.position]=deal(positions{:});
    [sitetrials.hemifield]=deal(hemifields{:});
    [sitetrials.fixation]=deal(fixations{:});
    [sites(i).trial]=sitetrials;
    
    sites(i)=ph_LR_to_CI(lfp_tfa_cfg,sites(i));  %% convert... 
    
    %% now loop through each trial for this site
    for t = 1:length(sites(i).trial)
        % convert hand and space information into string labels (for some reason)
        hf=sites(i).trial(t).hemifield;
        
        % reach space
        if hf == -1
            reach_space = 'I';
        elseif hf == 1
            reach_space = 'C';
        else
            reach_space = 'N';
        end      
        
        rh = sites(i).trial(t).reach_hand; % 1 = left, 2 = right
        % reach hand
        if rh == 1
            reach_hand = 'I';
        elseif rh == 2
            reach_hand = 'C';
        else
            reach_hand = 'N';  % no hand labeling
        end
        hs_label=[reach_hand 'H ' reach_space 'S'];        
        
        site_lfp.trials(t).reach_hand  = reach_hand;
        site_lfp.trials(t).reach_space = reach_space;
        site_lfp.trials(t).hndspc_lbl  = hs_label;
              
        
        %% retrieve LFP data
        
        start_time = (sites(i).trial(t).TDT_LFPx_tStart); % trial start time
        fs = sites(i).trial(t).TDT_LFPx_SR; % sample rate
        LFP = sites(i).trial(t).LFP; % LFP data
        ts = (1/fs); % sample time
        nsamples = numel(LFP);
        end_time = start_time + (ts*(nsamples-1));
        timestamps = linspace(start_time, end_time, nsamples);        
        
        site_lfp.trials(t).lfp_data    = LFP;
        site_lfp.trials(t).time        = timestamps;
        site_lfp.trials(t).fsample     = fs;
        site_lfp.trials(t).tsample     = ts;

        %% need information about which trial it was(originally?)
        site_lfp.trials(t).n  = sites(i).trial(t).n;
        
        % save retrieved data into struct
                
        perturbation = sites(i).trial(t).perturbation; % 0 = control
        if isnan(perturbation)
            perturbation = 0;
        end        
        site_lfp.trials(t).perturbation  = perturbation;
        site_lfp.trials(t).completed    = sites(i).trial(t).completed;
        site_lfp.trials(t).success      = sites(i).trial(t).success;
        site_lfp.trials(t).type         = sites(i).trial(t).type;
        site_lfp.trials(t).effector     = sites(i).trial(t).effector;
        site_lfp.trials(t).run          = sites(i).trial(t).run;
        site_lfp.trials(t).block        = sites(i).trial(t).block;
        site_lfp.trials(t).choice_trial = sites(i).trial(t).choice;
                
        % flag to mark noisy trials, default False, filled in by lfp_tfa_reject_noisy_lfp.m
        site_lfp.trials(t).noisy = 0;
        
        % get state onset times and onset samples - test and delete
        site_lfp.trials(t).states = struct();
        
        for s = 1:length(sites(i).trial(t).states)
            % get state ID
            state_id = sites(i).trial(t).states(s);
            % get state onset time
            state_onset = sites(i).trial(t).states_onset(sites(i).trial(t).states == state_id);
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
    
    %%% Noise rejection - should this be included within processing check this? %%%
    %state_filt_lfp(i) = lfp_tfa_reject_noisy_lfp( state_lfp(i), lfp_tfa_cfg.noise );
    
    %% Time frequency spectrogram calculation
    site_lfp = ecg_bna_compute_site_lfp_tfr( site_lfp, lfp_tfa_cfg );
    %site_lfp = lfp_tfa_compute_site_tfr( site_lfp, lfp_tfa_cfg );
    
    % Noise rejection
    site_lfp = lfp_tfa_reject_noisy_lfp_trials( site_lfp, lfp_tfa_cfg.noise );
    
    % Baseline power calculation - this needs to move to somewhere entirely different (?)
    site_lfp = lfp_tfa_compute_site_baseline( site_lfp, session_info, lfp_tfa_cfg );
    
    allsites_lfp = [allsites_lfp, site_lfp];    
end

%% calculate cross power spectrum between sites and sync measure spectrogram - NOT WORKING at the moment
if any(strcmp(lfp_tfa_cfg.analyses, 'sync')) || any(strcmp(lfp_tfa_cfg.analyses, 'syncsp'))
    % prepare results folder
    results_fldr = fullfile(session_info.proc_results_fldr, 'crossspectrum');
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    % loop through each site
    for i = 1:length(allsites_lfp)-1
        site1_lfp = allsites_lfp(i);
        % pair a site
        for j = i+1:length(allsites_lfp)
            site2_lfp = allsites_lfp(j);
            fprintf('Computing cross power spectrum for site pair %s - %s\n', ...
                site1_lfp.site_ID, site2_lfp.site_ID);
            sitepair_crosspow = lfp_tfa_compute_sitepair_csd(site1_lfp, site2_lfp, lfp_tfa_cfg);
            % save data
            results_mat = fullfile(results_fldr, ['sites_crosspow_', sitepair_crosspow.sites{1} '-' sitepair_crosspow.sites{2} '.mat']);
            save(results_mat, 'sitepair_crosspow', '-v7.3');
            
        end
    end
end

%% calculate sync measure spectrum
%     if any(strcmp(lfp_tfa_cfg.analyses, 'syncspctrm'))
%         % loop through each site
%         for i = 1:length(allsites_lfp)-1
%             site1_lfp = allsites_lfp(i);
%             % pair a site
%             for j = i+1:length(allsites_lfp)
%                 site2_lfp = allsites_lfp(j);
%                 fprintf('Computing sync spectrum for site pair %s - %s\n', ...
%                     site1_lfp.site_ID, site2_lfp.site_ID);
%                 % compute ppc spectrum between sitepair
%                 % get the trial conditions for this session
%                 conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg, ...
%                     {session_info.Preinj_blocks, session_info.Postinj_blocks});
%                 sitepair_syncspctrm = lfp_tfa_sitepair_avg_syncspctrum(site1_lfp, site2_lfp, conditions, lfp_tfa_cfg);
%             end
%         end
%     end

end

