function session_ecg = ecg_bna_read_preproc_ECG( session_info, plottrials )

% ecg_bna_read_preproc_ECG - function to read in the processed ECG and
% compute the time frequency spectrogram for each trial
%
% USAGE:
%	session_ecg = ecg_bna_read_preproc_ECG( session_info, plottrials )
% INPUTS: 
%   session_info    - struct containing relevant information about the
%   session to be analysed
%       Required fields:
%           Input_ECG_preproc   - path to the file containing preprocessed
%           ECG data containing trial information, ECG timeseries data and
%           timestamps for all trials in a session. For example,
%           Y:\Projects\PhysiologicalRecording\ephys\ECG_reaching\by_block_*.mat
%           Input_ECG           - path to the file containing ECG Rpeak data for a session. For example,
%           Y:\Projects\PhysiologicalRecording\Data\[Mpnkey]\[Date]\[Date]_ecg.mat
%           proc_ecg_fldr       - path to store the resulting variables
%           session             - name of the session
%       Optionl fields:
%           Preinj_blocks       - blocks to be considered as pre injection
%           Postinj_blocks      - blocks to be considred as post-injection
%
% OUTPUTS: 
%   session_ecg     - struct which stores the trial information, ECG data
%   and timestmps, and ECG Rpeak information for all available trials of a
%   session
%
% REQUIRES: ecg_bna_get_block_Rpeak_times
%
% See also ecg_bna_read_combined_ECG
    
    close all; 
    
    if nargin < 2
        plottrials = 0;
    end
    
    % struct to save data for a site
    session_ecg = struct();
        
    % functionality to read multiple by_block files
    session_info.Input_ECG_preproc = {session_info.Input_ECG_preproc};
    combined_Blocks = cell(1, length(session_info.Input_ECG_preproc));
    for s = 1:length(combined_Blocks)
        % Read input ECG file
        if ~exist(session_info.Input_ECG_preproc{s}, 'file')
            fprintf('No file found: %s\n', ...
                session_info.Input_ECG_preproc);
            return;
        else
            load(session_info.Input_ECG_preproc{s}, 'by_block');
            combined_Blocks{s} = by_block;
        end
        
    end
    
    
    % prepare results folder
    results_fldr = fullfile(session_info.proc_ecg_fldr);
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    if isfield(session_info, 'Input_ECG')
        if ~exist(session_info.Input_ECG, 'file')
            fprintf('No file found \n%s\n', ...
                session_info.Input_ECG);
            return;
        end
        load(session_info.Input_ECG);
        if exist('out', 'var')
            block_Rpeaks = out;
            clear out;
        end        
    end
          
    % for future use
%     usable_sites_table = table;
%     if ~isempty(ecg_bna_cfg.sites_info)
%        usable_sites_table = ecg_bna_cfg.sites_info;
%     end
    comp_trial = 0; % iterator for completed trials 
    
    % save data inside struct 
    % first loop through each site
    for b = 1:length(by_block)
        
        % get info about site
        % for future use
            % find if this site's entry is available in usable_sites_table
    %         if isempty(usable_sites_table(strcmp(usable_sites_table.Site_ID, ...
    %                 sites(i).site_ID),:))
    %             continue;
    %         end
            
        % check if ECG peaks exists for this block
        if length(block_Rpeaks) >= b
            block_Rpeak = block_Rpeaks(b);
            if isempty(block_Rpeak) || isempty(block_Rpeak.Rpeak_t)
                continue;
            end
        else
            continue;
        end

        % block ECG timestamps
        block_ecg_timestamps = [];
        trial_ecg_timestamps = [];           


        % for future use
        % get 'Set' entry from usable_sites_table
%         site_lfp.dataset = usable_sites_table(...
%             strcmp(usable_sites_table.Site_ID, sites(i).site_ID), :).Set(1);
        session_ecg.session = session_info.session;

        trial = by_block(b).trial;
        if isempty(trial)
            fprintf('No trials found for block, %g\n', b); 
            continue;
        end
        
        fprintf('=============================================================\n');
        fprintf('Reading ECG for block, %g\n', trial(1).block);
        fprintf('Number of trials %g\n', length(trial));

          

        %% get information common to all sites for a session

        % now loop through each trial for this site
        for t = 1:length(trial)

            block = trial(t).block;
            run = trial(t).run;
            type = trial(t).type;
            effector = trial(t).effector;
            completed = trial(t).completed;
            success = trial(t).success;
            choice_trial = trial(t).choice;
            fix_pos = trial(t).fix_pos;
            tar_pos = trial(t).tar_pos;
            sac_off = trial(t).sac_off;
            % for future use
            % check if the block is usable
%                 if isempty(usable_sites_table(strcmp(usable_sites_table.Site_ID, ...
%                         sites(i).site_ID) && usable_sites_table.Block == block))
%                     continue;
%                 end
            
            perturbation = nan;

            if isfield(session_info, 'Preinj_blocks') && ...
                ~isempty(session_info.Preinj_blocks) && ...
                ismember(block, session_info.Preinj_blocks)
                perturbation = 0;
%             elseif isfield(trial, 'perturbation') && ...
%                 ~isempty(trial(t).perturbation)
%                 perturbation = trial(t).perturbation;
            elseif exist('ses', 'var') && ...
                    isfield(ses, 'first_inj_block') && ...
                    ~isempty(ses.first_inj_block) && ...
                    block < ses.first_inj_block
                perturbation = 0;  
            end
            
            if isnan(perturbation)
                if isfield(session_info, 'Postinj_blocks') && ...
                    ~isempty(session_info.Postinj_blocks) && ...
                    ismember(block, session_info.Postinj_blocks)
                    perturbation = 1;
                elseif exist('ses', 'var') && ...
                    isfield(ses, 'first_inj_block') && ...
                    ~isempty(ses.first_inj_block) && ...
                    block >= ses.first_inj_block
                    perturbation = 1;
                elseif isfield(trial, 'perturbation') && ...
                    isempty(trial(t).perturbation)
                    perturbation = trial(t).perturbation;
                end
            end
            
            start_time = trial(t).TDT_ECG1_tStart; % trial start time
            fs = trial(t).TDT_ECG1_SR; % sample rate
            ts = (1/fs); % sample time
            ECG = trial(t).TDT_ECG1; % ecg data
            nsamples = numel(ECG);
            end_time = start_time + (ts*(nsamples-1));
            timestamps = linspace(start_time, end_time, nsamples);  
            if ~isempty(trial_ecg_timestamps)
                trial_ecg_timestamps = ...
                    timestamps - timestamps(1) + ts + trial_ecg_timestamps(end); 
            else
                trial_ecg_timestamps = timestamps - timestamps(1);
            end
%                         block_ecg_timestamps = [block_ecg_timestamps, ...
%                             trial_ecg_timestamps];
            % save retrieved data into struct
            comp_trial = comp_trial + 1;
            session_ecg.trials(comp_trial).completed = completed;
            session_ecg.trials(comp_trial).type = type;
            session_ecg.trials(comp_trial).effector = effector;
            session_ecg.trials(comp_trial).run = run;
            session_ecg.trials(comp_trial).block = block;
            session_ecg.trials(comp_trial).dataset = [];
            session_ecg.trials(comp_trial).choice_trial = choice_trial;
            session_ecg.trials(comp_trial).success = success;
            session_ecg.trials(comp_trial).reach_hand = 0;
            session_ecg.trials(comp_trial).reach_space = 0;
            session_ecg.trials(comp_trial).fix_pos = fix_pos;
            session_ecg.trials(comp_trial).eye_pos = tar_pos + sac_off;
            session_ecg.trials(comp_trial).hndspc_lbl  = [];
            session_ecg.trials(comp_trial).time = timestamps;
            session_ecg.trials(comp_trial).ecg_data = ECG;
            session_ecg.trials(comp_trial).fsample  = fs;
            session_ecg.trials(comp_trial).tsample = ts;
            session_ecg.trials(comp_trial).tstart = start_time;
            session_ecg.trials(comp_trial).trialperiod = ...
                [trial_ecg_timestamps(1) trial_ecg_timestamps(end)];
            session_ecg.trials(comp_trial).perturbation  = perturbation;
            % flag to mark noisy trials
            session_ecg.trials(comp_trial).noisy = ~completed;

            % get state onset times and onset samples - test and delete
            session_ecg.trials(comp_trial).states = struct();
            st_idx = 0;
            for st = 1:length(trial(t).states)
                % get state ID
                state_id = trial(t).states(st);
                st_idx = st_idx + 1;
                % get state onset time
                state_onset = trial(t).states_onset(trial(t).states == ...
                    state_id);
                % get sample number of state onset time
                state_onset_sample = find(abs(timestamps - state_onset(1)) == ...
                    min(abs(timestamps - state_onset(1))), 1);
                % save into struct
                session_ecg.trials(comp_trial).states(st_idx).id = state_id;
                session_ecg.trials(comp_trial).states(st_idx).onset_t  = state_onset(1);
                session_ecg.trials(comp_trial).states(st_idx).onset_s  = state_onset_sample;
            end
            %end


        end

        session_ecg = ecg_bna_get_block_Rpeak_times( session_ecg, block_Rpeak, block, plottrials, results_fldr );



    %%% Noise rejection - should this be included within processing check this? %%%
    %state_filt_lfp(i) = ecg_bna_reject_noisy_lfp( state_lfp(i), ecg_bna_cfg.noise );

        
    end
    
    results_mat = fullfile(results_fldr, ['session_ecg_' session_info.session '.mat']);
    save(results_mat, 'session_ecg', '-v7.3');
    
    
end

