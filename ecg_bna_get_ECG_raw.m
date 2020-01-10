function [ site_lfp ] = ecg_bna_get_ECG_raw( site_lfp, Input_ECG_folder )
% ecg_bna_get_ECG_raw - gets the raw ECG data for all trials of a
% session and combines it with LFP information in the incoming struct
%
% USAGE:
%	[ site_lfp ] = ecg_bna_get_ECG_raw( site_lfp, Input_ECG_folder ) 
%
% INPUTS:
%		site_lfp            - struct containing LFP data for all trials of 
%		a session from a single site
%       Input_ECG_folder    - path to the folder containing mat files which
%       store the combined LFP-ECG information for a session. 
%       This is usually of the form Y:\Data\[Monkey]\[Date]
%
% OUTPUTS:
%		site_lfp            - same as input structure but with additional
%		fields to store the ECG raw data for all trials in a block
% 
% See also ecg_bna_process_combined_LFP_ECG, ecg_bna_get_ECG_peaks 

% list mat files inside the given directory
block_ECG_files = dir(fullfile(Input_ECG_folder, '*.mat'));
ECG_filenames = {block_ECG_files.name};
ECG_filenames = cat(1, ECG_filenames{:});
ECG_blocks = ECG_filenames(:,end-5:end-4);
ECG_blocks = str2num(ECG_blocks);

%loop through each run
for b = (unique([site_lfp.trials.block]))
    fprintf('Reading ECG raw data for block %g\n', b);
    % load raw ECG for this block
    load(fullfile(Input_ECG_folder, ECG_filenames(ECG_blocks ==b,:)));
    if exist('trial', 'var')
        trials_ECG = trial;
        clear trial;
    end
    
    % trial ids for this block
    trials_idx = find([site_lfp.trials.block] == b);
    for t = 2:length(trials_idx) 
        % ECG sampling rate
        ECG_SR = trials_ECG(t).TDT_ECG1_samplingrate;
        ECG_ts = 1/ECG_SR;
        trial_ECG_raw = trials_ECG(t).TDT_ECG1;
        nsamples = length(trial_ECG_raw);
        start_time = site_lfp.trials(trials_idx(t)).tstart;
        end_time = start_time + ECG_ts * (nsamples-1);        
        ECG_timestamps = linspace(start_time, end_time, nsamples);
        
        % resample ECG data
        ECG_timeseries = timeseries(trial_ECG_raw, ECG_timestamps);
        trial_timestamps = site_lfp.trials(trials_idx(t)).time;
        resampled_ECG_timeseries = resample(ECG_timeseries, trial_timestamps);
        trial_ECG_raw = permute(resampled_ECG_timeseries.Data, [2 3 1]);
        
        % store ECG data
        site_lfp.trials(trials_idx(t)).ecg_data = trial_ECG_raw;
       
        
        
    end
    
end
