function [One_stream]=ecg_bna_synchronize_with_trigger(Blockoffsets,trcell,popcell)
%ecg_bna_synchronize_with_trigger(trcell,popcell,Rpeaks)
%Bockoffsets=Rpeaks([Rpeaks.block] == z).offset
%Blockoffsets=Blockoffsets-Blockoffsets(1);
% 0. Prepare data variables
states_onset               = {trcell.states_onset};
states                     = {trcell.states};
rec_start                  = {trcell.TDT_LFPx_t0_from_rec_start}; % this is the relative time (First_trial_INI is NOT included (?))
block_nums                 = {trcell.block};
% compute RR-intervals
%         valid_RRinterval_ends      = single([Rpeaks(b).(['RPEAK_ts' cfg.condition(c).Rpeak_field])]);
%         valid_RRinterval_starts    = single(valid_RRinterval_ends - [Rpeaks(b).(['RPEAK_dur' cfg.condition(c).Rpeak_field])]);
% 0. figure out RR-intervals lying within trials


One_stream.state_onsets     = cellfun(@(x,y,z) x+y+Blockoffsets(z), states_onset,  rec_start, block_nums, 'Uniformoutput', false);

One_stream.state_onsets     =[One_stream.state_onsets{:}];
One_stream.states           =[states{:}];

One_stream.trial_starts     = One_stream.state_onsets(One_stream.states==2);
One_stream.trial_ends       = One_stream.state_onsets(One_stream.states==90);
%One_stream.trial_ends      = cellfun(@(x,y,z) x+y+Bockoffsets(z), state90_times, rec_start, block_nums);

%% is there a way, starts and ends can be of different length? hope not

if nargin>2
    arrival_times              = {popcell.arrival_times};
    arrival_times              = cellfun(@(x,y,z) x+y+Blockoffsets(z), arrival_times, rec_start, block_nums, 'Uniformoutput', false);
    
    One_stream.AT=vertcat(arrival_times{:});
    One_stream.AT(One_stream.AT>One_stream.trial_ends(end))=[];
end
end