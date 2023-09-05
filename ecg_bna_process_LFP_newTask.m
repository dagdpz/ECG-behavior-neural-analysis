function [sessions_info, allsites_lfp]= ecg_bna_process_LFP_newTask( sessions_info, ecg_bna_cfg )

% ecg_bna_process_LFP_newTask - function to read in the site-wise LFP data for all
% sites recorded in a session, compute the LFP time frequency spectrogram,
% detect the noisy trials, and compute site-wise baseline power
%
% USAGE:
%	sessions_info = ecg_bna_process_LFP_newTask( sessions_info, ecg_bna_cfg )
%
% INPUTS:
%       sessions_info        - structure containing information about the
%       session to be processed, see settings/lfp_tfa_settings_example
%       Required fields:
%           Input               - path to the file which contains the
%                               trial-wise LFP data for all sites recorded
%                               in a session
%           proc_results_folder - folder where the results has to be
%           stored, see lfp_tfa_define_settings
%       ecg_bna_cfg         - structure containing configurations for
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
%		sessions_info            - same as input structure sessions_info
%
% REQUIRES: lfp_tfa_compute_site_tfr, lfp_tfa_reject_noisy_lfp_trials,
% lfp_tfa_compute_site_baseline, lfp_tfa_compute_sitepair_csd
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_define_settings,
% lfp_tfa_global_states, lfp_tfa_reject_noisy_lfp_trials,
% lfp_tfa_compute_site_baseline


% prepare results folder
results_fldr = fullfile(sessions_info.proc_results_fldr);
if ~exist(results_fldr, 'dir')
    mkdir(results_fldr);
end

d= dir(['Y:\Projects\Pulv_bodysignal\ephys\',ecg_bna_cfg.spikes_version filesep]);
sites_names = (strfind({d.name},[ 'sites_',sessions_info.session,'_']));
site_idx = find(~cellfun(@isempty,sites_names));
names = {d.name};
clear d;
% load(sessions_info.Input_LFP{:}, 'sites');
% for i = 1:length(site_idx)
%     sites(i)=load(['Y:\Projects\Pulv_bodysignal\ephys\',ecg_bna_cfg.spikes_version filesep,d(site_idx(i)).name],'sites');
% end
% adding the separate trial folders:
load(sessions_info.Input_trials, 'trials');

% structure array to store lfp data for all sites
% to be used for cross power spectrum calculation
allsites_lfp = [];

for i = 1:length(site_idx)
    
    load(['Y:\Projects\Pulv_bodysignal\ephys\',ecg_bna_cfg.spikes_version filesep,names{site_idx(i)}],'sites');
    % struct to save data for a site
    site_lfp = struct();
    fprintf('=============================================================\n');
    fprintf('Processing site, %s\n', sites(i).site_ID);
    site_lfp.session = sites(i).site_ID(1:12);
    site_lfp.site_ID = sites(i).site_ID;
    site_lfp.target = sites(i).target;
    site_lfp.recorded_hemisphere = upper(sites(i).target(end));
    %site_lfp.ref_hemisphere = ecg_bna_cfg.ref_hemisphere;
    site_lfp.xpos = sites(i).grid_x;
    site_lfp.ypos = sites(i).grid_y;
    site_lfp.zpos = sites(i).electrode_depth;
    
%     sitetrials=[sites(i).trial];
    sitetrials= ecg_bna_get_site_trials(sites(i),trials);
    
    positions=[sitetrials.tar_pos]-[sitetrials.fix_pos];
    hemifields=num2cell(sign(real(positions)));
    fixations=num2cell([sitetrials.fix_pos]);
    positions=num2cell(positions);
    [sitetrials.position]=deal(positions{:});
    [sitetrials.hemifield]=deal(hemifields{:});
    [sitetrials.fixation]=deal(fixations{:});
    [sites(i).trial]=sitetrials;
    
    sites(i)=ph_LR_to_CI(ecg_bna_cfg,sites(i));  %% convert... 
    
    %% now loop through each trial for this site to do what exactly?
    for t = 1:length(sites(i).trial)
%         % convert hand and space information into string labels (for some reason)
%         hf=sites(i).trial(t).hemifield;
%         
%         % reach space
%         if hf == -1
%             reach_space = 'I';
%         elseif hf == 1
%             reach_space = 'C';
%         else
%             reach_space = 'N';
%         end      
%         
%         rh = sites(i).trial(t).reach_hand; % 1 = left, 2 = right
%         % reach hand
%         if rh == 1
%             reach_hand = 'I';
%         elseif rh == 2
%             reach_hand = 'C';
%         else
%             reach_hand = 'N';  % no hand labeling
%         end
%         hs_label=[reach_hand 'H ' reach_space 'S'];        
%         
%         site_lfp.trials(t).reach_hand  = reach_hand;
%         site_lfp.trials(t).reach_space = reach_space;
%         site_lfp.trials(t).hndspc_lbl  = hs_label;
              
        
        %% retrieve LFP data
        
        start_time = (sites(i).trial(t).TDT_LFPx_tStart); % trial start time
        fs = sites(i).trial(t).TDT_LFPx_SR; % sample rate
%         LFP = sites(i).trial(t).LFP; % LFP data
        LFP = sites(i).LFP; % LFP data
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
        
        trial_start_t = site_lfp.trials(t).states([site_lfp.trials(t).states.id] == ecg_bna_cfg.trialinfo.start_state).onset_t + ecg_bna_cfg.trialinfo.ref_tstart;
        trial_end_t   = site_lfp.trials(t).states([site_lfp.trials(t).states.id] == ecg_bna_cfg.trialinfo.end_state).onset_t + ecg_bna_cfg.trialinfo.ref_tend;
        site_lfp.trials(t).trialperiod = [trial_start_t, trial_end_t];
    end
        
    %% Time frequency spectrogram calculation
    site_lfp = ecg_bna_compute_site_lfp_tfr( site_lfp, ecg_bna_cfg );
    
    % Noise rejection
    site_lfp = lfp_tfa_reject_noisy_lfp_trials( site_lfp, ecg_bna_cfg.noise );
    
    % Baseline power calculation - this needs to move to somewhere entirely different (?)
    % site_lfp = lfp_tfa_compute_site_baseline( site_lfp, sessions_info, ecg_bna_cfg );
    
    allsites_lfp = [allsites_lfp, site_lfp];    
end

%% calculate cross power spectrum between sites and sync measure spectrogram - PROBABLY NOT WORKING at the moment
if any(strcmp(ecg_bna_cfg.analyses, 'sync')) || any(strcmp(ecg_bna_cfg.analyses, 'syncsp'))
    % prepare results folder
    results_fldr = fullfile(sessions_info.proc_results_fldr, 'crossspectrum');
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
            sitepair_crosspow = lfp_tfa_compute_sitepair_csd(site1_lfp, site2_lfp, ecg_bna_cfg);
            % save data
            results_mat = fullfile(results_fldr, ['sites_crosspow_', sitepair_crosspow.sites{1} '-' sitepair_crosspow.sites{2} '.mat']);
            save(results_mat, 'sitepair_crosspow', '-v7.3');
            
        end
    end
end

%% calculate sync measure spectrum
%     if any(strcmp(ecg_bna_cfg.analyses, 'syncspctrm'))
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
%                 conditions = lfp_tfa_compare_conditions(ecg_bna_cfg, ...
%                     {sessions_info.Preinj_blocks, sessions_info.Postinj_blocks});
%                 sitepair_syncspctrm = lfp_tfa_sitepair_avg_syncspctrum(site1_lfp, site2_lfp, conditions, ecg_bna_cfg);
%             end
%         end
%     end

end

