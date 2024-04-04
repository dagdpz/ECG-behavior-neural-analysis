function [trial_starts_one_stream, trial_ends_one_stream,AT]=ecg_bna_get_trial_starts_and_ends(trcell,popcell,Rpeaks)
        % 0. Prepare data variables
        states_onset               = {trcell.states_onset};
        states                     = {trcell.states};
        rec_start                  = {trcell.TDT_ECG1_t0_from_rec_start}; % this is the relative time (First_trial_INI is included (?))
        block_nums                 = {trcell.block};
        arrival_times              = {popcell.arrival_times};
        state2_times               = cellfun(@(x,y) x(y == 2), states_onset, states, 'Uniformoutput', false); % trial starts = state 2
        state90_times              = cellfun(@(x,y) x(y == 90), states_onset, states, 'Uniformoutput', false); % trial ends = state 90
        % compute RR-intervals
%         valid_RRinterval_ends      = single([Rpeaks(b).(['RPEAK_ts' cfg.condition(c).Rpeak_field])]);
%         valid_RRinterval_starts    = single(valid_RRinterval_ends - [Rpeaks(b).(['RPEAK_dur' cfg.condition(c).Rpeak_field])]);
        % 0. figure out RR-intervals lying within trials
        
        
         % add trial onset time to each spike so its basically one stream again
         % also, make sure spikes aren't counted twice (because previous trial is appended in beginning;
         % (this will also remove overlapping spikes)
         % state time within trial(x) + trial onset relative to rec_start (y) + block separator (stored in triggers)
        
        trial_starts_one_stream    = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, state2_times,  rec_start, block_nums);
        trial_ends_one_stream      = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, state90_times, rec_start, block_nums);
        arrival_times              = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, arrival_times, rec_start, block_nums, 'Uniformoutput', false);
        
        AT=vertcat(arrival_times{:});
        AT(AT>trial_ends_one_stream(end))=[];
    end