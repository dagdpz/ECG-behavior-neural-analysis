function session_ecg = ecg_bna_read_preproc_ECG_simple( session_info, plottrials )

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
if ~iscell(session_info.Input_ECG_preproc)
    session_info.Input_ECG_preproc = {session_info.Input_ECG_preproc};
end
if ~exist(session_info.Input_ECG_preproc{1}, 'file')
    fprintf('No file found: %s\n', session_info.Input_ECG_preproc);
    return;
else
    combined_Blocks = cell(1, length(session_info.Input_ECG_preproc));
    for s = 1:length(combined_Blocks)
        % Read input ECG file
        load(session_info.Input_ECG_preproc{s}, 'by_block');

        % quick but hopefully soon no more necessary fix of same block being
        % split (should be fixed in spike analysis, reason is we have split
        % blocks in combined monkeypsych and TDT files)
        for b=1:numel(by_block)
        all_blocks(b)=by_block(b).trial(1).block;
        end
        multiple_blocks=fliplr(find([0 diff(all_blocks)==0])); % flip cause we wanna go backwards in the loop ahead
        for b=multiple_blocks
            for t=1:numel(by_block(b).trial)
                by_block(b).trial(t).n=by_block(b).trial(t).n+by_block(b-1).trial(end).n;
                by_block(b).trial(t).run_onset_time=by_block(b).trial(t).n+by_block(b-1).trial(end).run_onset_time;
                if isfield(by_block(b).trial(t),'TDT_ECG1_t0_from_rec_start')
                    by_block(b).trial(t).TDT_ECG1_t0_from_rec_start=by_block(b).trial(t).n+by_block(b-1).trial(end).TDT_ECG1_t0_from_rec_start;
                end
                if isfield(by_block(b).trial(t),'TDT_ECG4_t0_from_rec_start')
                    by_block(b).trial(t).TDT_ECG4_t0_from_rec_start=by_block(b).trial(t).n+by_block(b-1).trial(end).TDT_ECG4_t0_from_rec_start;
                end
                if isfield(by_block(b).trial(t),'TDT_LFPx_t0_from_rec_start')
                    by_block(b).trial(t).TDT_LFPx_t0_from_rec_start=by_block(b).trial(t).n+by_block(b-1).trial(end).TDT_LFPx_t0_from_rec_start;
                end
            end
            [by_block(b-1).trial]=[by_block(b-1).trial by_block(b).trial];
        end
        by_block(multiple_blocks)=[];
        
        combined_Blocks{s} = by_block;
    end
end

% prepare results folder
results_fldr = fullfile(session_info.proc_ecg_fldr);
if ~exist(results_fldr, 'dir')
    mkdir(results_fldr);
end


T = 0; % iterator for completed trials


% save data inside struct
% first loop through each file
for i_Files = 1:length(combined_Blocks)
    by_block_ECG = combined_Blocks{i_Files};
    for b = 1:length(by_block_ECG)
        % get info about site
        
        trial_ecg_timestamps = [];
        session_ecg.session = session_info.session;
        
        trial = by_block_ECG(b).trial;
        if isempty(trial)
            fprintf('No trials found for block, %g\n', b);
            continue;
        end
        
        fprintf('=============================================================\n');
        fprintf('Reading ECG for block %g\n', trial(1).block);
        fprintf('Number of trials %g\n', length(trial));
                
        % now loop through each trial for this site
        for t = 1:length(trial)
            
            block = trial(t).block;
            run = trial(t).run;
            type = trial(t).type;
            original_trial_n = trial(t).n;
            effector = trial(t).effector;
            completed = trial(t).completed;
            success = trial(t).success;
            choice_trial = trial(t).choice;
            fix_pos = trial(t).fix_pos;
            tar_pos = trial(t).tar_pos;
            sac_off = trial(t).sac_off;
            
            perturbation = nan;
            
            if isfield(session_info, 'Preinj_blocks') && ~isempty(session_info.Preinj_blocks) && ismember(block, session_info.Preinj_blocks)
                perturbation = 0;
                %             elseif isfield(trial, 'perturbation') && ...
                %                 ~isempty(trial(t).perturbation)
                %                 perturbation = trial(t).perturbation;
            elseif exist('ses', 'var') && (isempty(ses.first_inj_block) || block < ses.first_inj_block)
                perturbation = 0;
            end
            
            if isnan(perturbation) %% this part makes no sense to me
                if isfield(session_info, 'Postinj_blocks') && ~isempty(session_info.Postinj_blocks) && ismember(block, session_info.Postinj_blocks)
                    perturbation = 1;
                elseif exist('ses', 'var') && block >= ses.first_inj_block
                    perturbation = 1;
                elseif isfield(trial, 'perturbation') && ~isempty(trial(t).perturbation)
                    perturbation = sign(trial(t).perturbation);
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
                trial_ecg_timestamps = timestamps - timestamps(1) + ts + trial_ecg_timestamps(end);
            else
                trial_ecg_timestamps = timestamps - timestamps(1);
            end
            % save retrieved data into struct
            T = T + 1;
            %session_ecg.trials(T).time_from_rec_start = timestamps+trial(t).TDT_ECG1_t0_from_rec_start; %% LS 20220328 THIS IS FROM REC ONSET
            session_ecg.trials(T).t0_from_rec_start = trial(t).TDT_ECG1_t0_from_rec_start; %% LS 20220328 THIS IS FROM REC ONSET
            session_ecg.trials(T).completed = completed;
            session_ecg.trials(T).type = type;
            session_ecg.trials(T).effector = effector;
            session_ecg.trials(T).run = run;
            session_ecg.trials(T).block = block;
            session_ecg.trials(T).dataset = [];
            session_ecg.trials(T).choice_trial = choice_trial;
            session_ecg.trials(T).success = success;
            session_ecg.trials(T).reach_hand = 0;
            session_ecg.trials(T).reach_space = 0;
            session_ecg.trials(T).fix_pos = fix_pos;
            session_ecg.trials(T).eye_pos = tar_pos + sac_off;
            session_ecg.trials(T).hndspc_lbl  = [];
            %session_ecg.trials(T).time = timestamps;
            session_ecg.trials(T).ecg_data = ECG;
            session_ecg.trials(T).fsample  = fs;
            session_ecg.trials(T).tsample = ts;
            session_ecg.trials(T).tstart = start_time;
            session_ecg.trials(T).trialperiod = [trial_ecg_timestamps(1) trial_ecg_timestamps(end)];
            session_ecg.trials(T).perturbation  = perturbation;
            session_ecg.trials(T).n  = original_trial_n;
            
            % flag to mark noisy trials
            session_ecg.trials(T).noisy = ~completed;
            
            % get state onset times and onset samples - test and delete
            session_ecg.trials(T).states = struct();
            st_idx = 0;
            for st = 1:length(trial(t).states)
                % get state ID
                state_id = trial(t).states(st);
                st_idx = st_idx + 1;
                % get state onset time
                state_onset = trial(t).states_onset(trial(t).states == state_id);
                % get sample number of state onset time
                %state_onset_sample = find(abs(timestamps - state_onset(1)) == min(abs(timestamps - state_onset(1))), 1);
                % save into struct
                session_ecg.trials(T).states(st_idx).id = state_id;
                session_ecg.trials(T).states(st_idx).onset_t  = state_onset(1);
                %session_ecg.trials(T).states(st_idx).onset_s  = state_onset_sample;
            end
        end
        
        %session_ecg = ecg_bna_get_block_Rpeak_times( session_ecg, block_Rpeak, block, plottrials, results_fldr );
    end
    
    results_mat = fullfile(results_fldr, ['session_ecg_' session_info.session '.mat']);
    save(results_mat, 'session_ecg');
end

end

