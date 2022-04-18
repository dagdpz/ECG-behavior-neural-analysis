function session_info = ecg_bna_process_combined_LFP_ECG( session_info, ecg_bna_cfg )
% ecg_bna_process_combined_LFP_ECG - function which performs the follwing
% tasks for each site recorded in a session
% 1. read in the processed lfp data
% 2. store the trials information and lfp raw data
% 3. call to compute LFP time frequency spectrograms
% 4. call to detect noisy trials
% 5. call to get Rpeak data for each trial
% 6. call to get ECG raw data for each trial
% 7. store the trial info, LFP data and tfs, ECG raw data and Rpeaks
% compute the time frequency spectrogram for each trial
%
% USAGE:
%	session_info = ecg_bna_process_combined_LFP_ECG( session_info, ecg_bna_cfg )
%
% INPUTS:
%   session_info        - struct containing information about all sessions
%   to be processed
%   ecg_bna_cfg         - struct containing settings for analysis
% OUTPUTS:
%   session_info        - a copy of input struct session_info
%
% REQUIRES: ecg_bna_get_ECG_raw, ecg_bna_get_ECG_peaks,
% lfp_tfa_compute_site_tfr, lfp_tfa_reject_noisy_lfp_trials,
% lfp_tfa_compute_site_baseline,
%
% See also lfp_tfa_process_LFP, ecg_bna_read_preproc_ECG, ecg_bna_read_combined_ECG

close all;

if isfield(session_info, 'Input_LFP')
    if ~iscell(session_info.Input_LFP)
        session_info.Input_LFP = {session_info.Input_LFP};
    end
    combined_sites = cell(1, length(session_info.Input_LFP));
    for s = 1:length(combined_sites)
        % Read input LFP file
        load(session_info.Input_LFP{s}, 'sites');
        combined_sites{s} = sites;
    end
end

if isfield(session_info, 'Input_ECG')
    block_ECG = load(session_info.Input_ECG);
    %         if exist('out', 'var')
    %             block_ECG = out;
    %             clear out;
    %         end
end

% prepare results folder
results_fldr = fullfile(session_info.proc_lfp_fldr);
if ~exist(results_fldr, 'dir')
    mkdir(results_fldr);
end



% for future use
%     usable_sites_table = table;
%     if ~isempty(ecg_bna_cfg.sites_info)
%        usable_sites_table = ecg_bna_cfg.sites_info;
%     end
comp_trial = 0; % iterator for completed trials

% save data inside struct
% first loop through each site
for i = 1:(length(combined_sites{1}))
    
    % get info about site
    % for future use
    % find if this site's entry is available in usable_sites_table
    %         if isempty(usable_sites_table(strcmp(usable_sites_table.Site_ID, ...
    %                 sites(i).site_ID),:))
    %             continue;
    %         end
    fprintf('=============================================================\n');
    fprintf('Processing site, %s\n', sites(i).site_ID);
    % for future use
    % get 'Set' entry from usable_sites_table
    %         site_lfp.dataset = usable_sites_table(...
    %             strcmp(usable_sites_table.Site_ID, sites(i).site_ID), :).Set(1);

% struct to save data for a site - SCARY this was outside the loop, causing
% trials to accumulate and be taken over for next site (in case of less
% trials for next site!)

site_lfp = struct();

    site_lfp.site_ID = sites(i).site_ID;
    site_lfp.target = sites(i).target;
    site_lfp.recorded_hemisphere = upper(sites(i).target(end));
    site_lfp.xpos = sites(i).grid_x;
    site_lfp.ypos = sites(i).grid_y;
    site_lfp.zpos = sites(i).electrode_depth;
    site_lfp.session = sites(i).site_ID(1:12);
    site_lfp.ref_hemisphere = ecg_bna_cfg.ref_hemisphere;
    
    % loop through each input LFP file  --???
    for s = 1:length(combined_sites)
        sites = combined_sites{s};
        
        % iterator for completed trials
        if s == 1
            comp_trial = 0;
        end
        % now loop through each trial for this site
        for t = 1:length(sites(i).trial)
            
            type = sites(i).trial(t).type;
            effector = sites(i).trial(t).effector;
            run = sites(i).trial(t).run;
            block = sites(i).trial(t).block;
            dataset = sites(i).trial(t).dataset;
            completed = sites(i).trial(t).completed;
            success = sites(i).trial(t).success;
            % for future use
            % check if the block is usable
            %                 if isempty(usable_sites_table(strcmp(usable_sites_table.Site_ID, ...
            %                         sites(i).site_ID) && usable_sites_table.Block == block))
            %                     continue;
            %                 end
            choice_trial = sites(i).trial(t).choice;
            reach_hand = sites(i).trial(t).reach_hand; % 1 = left, 2 = right
            
            perturbation = nan;
            if isfield(session_info, 'Preinj_blocks') && ~isempty(session_info.Preinj_blocks) && ismember(block, session_info.Preinj_blocks)
                perturbation = 0;
            elseif exist('ses', 'var') && ...
                    block < ses.first_inj_block
                perturbation = 0;
            end
            if isnan(perturbation)
                if isfield(session_info, 'Postinj_blocks') && ~isempty(session_info.Postinj_blocks) && ismember(block, session_info.Postinj_blocks)
                    perturbation = 1;
                elseif exist('ses', 'var') && ...
                        block >= ses.first_inj_block
                    perturbation = 1;
                end
            end
            if isnan(perturbation)
                perturbation = 0;
            end
            
            tar_pos = sites(i).trial(t).tar_pos;
            fix_pos = sites(i).trial(t).fix_pos;
            
            % reach space
            if sign(real(tar_pos) - real(fix_pos)) == -1
                reach_space = 'L';
            elseif sign(real(tar_pos) - real(fix_pos)) == 1
                reach_space = 'R';
            else
                reach_space = 'N';
            end
            
            % reach hand
            if reach_hand == 1
                reach_hand = 'L';
            elseif reach_hand == 2
                reach_hand = 'R';
            else
                reach_hand = 'N';  % no hand labeling
            end
            
            
            % assign hand-space for the trial
            if strcmp(site_lfp.ref_hemisphere, reach_space)
                if strcmp(site_lfp.ref_hemisphere, reach_hand)
                    hs_label = 'IH IS';
                else
                    hs_label = 'CH IS';
                end
            else
                if strcmp(site_lfp.ref_hemisphere, reach_hand)
                    hs_label = 'IH CS';
                else
                    hs_label = 'CH CS';
                end
            end
            
            start_time = (sites(i).trial(t).TDT_LFPx_tStart); % trial start time
            fs = sites(i).trial(t).TDT_LFPx_SR; % sample rate
            ts = (1/fs); % sample time
            LFP = sites(i).trial(t).LFP; % LFP data
            nsamples = numel(LFP);
            end_time = start_time + (ts*(nsamples-1));
            timestamps = linspace(start_time, end_time, nsamples);
            trial_onset_time = sites(i).trial(t).trial_onset_time;
            
            % save retrieved data into struct
            comp_trial = comp_trial + 1;
            site_lfp.trials(comp_trial).completed = completed;
            site_lfp.trials(comp_trial).success = success;
            site_lfp.trials(comp_trial).type = type;
            site_lfp.trials(comp_trial).effector = effector;
            site_lfp.trials(comp_trial).run = run;
            site_lfp.trials(comp_trial).block = block;
            site_lfp.trials(comp_trial).dataset = dataset;
            site_lfp.trials(comp_trial).choice_trial = choice_trial;
            site_lfp.trials(comp_trial).time = timestamps;
            site_lfp.trials(comp_trial).time_from_rec_start = timestamps+sites(i).trial(t).TDT_LFPx_t0_from_rec_start; %% LS 20220328 THIS IS FROM REC ONSET
            site_lfp.trials(comp_trial).lfp_data = LFP;
            site_lfp.trials(comp_trial).fsample  = fs;
            site_lfp.trials(comp_trial).tsample = ts;
            site_lfp.trials(comp_trial).tstart = start_time;
            site_lfp.trials(comp_trial).reach_hand  = reach_hand;
            site_lfp.trials(comp_trial).reach_space  = reach_space;
            site_lfp.trials(comp_trial).hndspc_lbl  = hs_label;
            site_lfp.trials(comp_trial).perturbation  = perturbation;
            % flag to mark noisy trials
            site_lfp.trials(comp_trial).noisy = ~completed;
            
            % get state onset times and onset samples - test and delete
            site_lfp.trials(comp_trial).states = struct();
            
            for st = 1:length(sites(i).trial(t).states)
                % get state ID
                state_id = sites(i).trial(t).states(st);
                % get state onset time
                state_onset = sites(i).trial(t).states_onset(sites(i).trial(t).states == state_id);
                % get sample number of state onset time
                state_onset_sample = find(abs(timestamps - state_onset) == min(abs(timestamps - state_onset)));
                % save into struct
                site_lfp.trials(comp_trial).states(st).id = state_id;
                site_lfp.trials(comp_trial).states(st).onset_t  = state_onset;
                site_lfp.trials(comp_trial).states(st).onset_s  = state_onset_sample;
            end
            
            trial_start_t = site_lfp.trials(comp_trial).states([site_lfp.trials(comp_trial).states.id] == ecg_bna_cfg.trialinfo.start_state).onset_t + ecg_bna_cfg.trialinfo.ref_tstart;
            trial_end_t   = site_lfp.trials(comp_trial).states([site_lfp.trials(comp_trial).states.id] == ecg_bna_cfg.trialinfo.end_state).onset_t + ecg_bna_cfg.trialinfo.ref_tend;
            site_lfp.trials(comp_trial).trialperiod = [trial_start_t, trial_end_t];
            site_lfp.trials(comp_trial).trial_onset_time = trial_onset_time;
            
        end
        
        if s == length(combined_sites)
            
            % Get ECG raw data
            if isfield(session_info, 'Input_ECG_combined')
                site_lfp = ecg_bna_get_ECG_raw( site_lfp, session_info.Input_ECG_combined );
            end
            
            % Get ECG spikes
            if exist('block_ECG', 'var')
                site_lfp = ecg_bna_get_ECG_peaks( site_lfp, block_ECG );
            end
            
        end
    end
    
    
    %%% Noise rejection - should this be included within processing check this? %%%
    %state_filt_lfp(i) = ecg_bna_reject_noisy_lfp( state_lfp(i), ecg_bna_cfg.noise );
    
    %% Time frequency spectrogram calculation
    site_lfp = lfp_tfa_compute_site_tfr( site_lfp, ecg_bna_cfg );
    
    % Noise rejection
    if ecg_bna_cfg.noise.detect
        site_lfp = lfp_tfa_reject_noisy_lfp_trials( site_lfp, ecg_bna_cfg.noise );
    end
    
    % Baseline power calculation
    site_lfp = lfp_tfa_compute_site_baseline( site_lfp, session_info, ecg_bna_cfg );
    
    % save data
    results_mat = fullfile(results_fldr, ['site_lfp_pow_' site_lfp.site_ID '.mat']);
    save(results_mat, 'site_lfp', '-v7.3');
    
end

end

