function ecg_bna_matchMPandTDT_times(YYYYMMDD)
% Example usage:
% ecg_bna_matchMPandTDT_times('20230518')
%
%
% Dependencies:
% ma1_check_timing_streams for dagdpz/ma1
%

% set up data folders based on session date
combined_folder = ['Y:\Data\Magnus_phys_combined_monkeypsych_TDT\' YYYYMMDD filesep];
bodySignals_folder = ['Y:\Data\BodySignals\ECG\Magnus\' YYYYMMDD filesep];

% set up combined and ECG data files based on session date
combined_files = dir([combined_folder '*.mat']);
ecg_file = dir([bodySignals_folder '*.mat']);

% find block numbers
block_id = cellfun(@(x) strfind(x, 'block_'), {combined_files.name}, 'UniformOutput', false);
blockNumbers = cellfun(@(x,y) str2double(x(y + (6:7))), {combined_files.name}, block_id);

% load ECG data
if length(ecg_file) == 1
    disp('Found one ECG mat file, we''re good! Proceeding...')
    ecg_filepath = [ecg_file.folder filesep ecg_file.name];
else
    error('Problem: more than 1 ECG file')
end

for blockNum = 1:length(combined_files)
    
    currFilePath = ...
        [combined_files(blockNum).folder filesep combined_files(blockNum).name];
    
    % trials are aligned to the 1st INI within this function
    trial = ma1_check_timing_streams(currFilePath, 0, 0, ecg_filepath, blockNumbers(blockNum));
    
    if length(trial(1).eye.tar) == 1
        disp(['Block ' num2str(blockNum) ' is a rest block. Skip, proceed with the next block'])
        continue
    end
    
    for trialNum = 1:length(trial)
    
        isVisibleTarget = nan(7,1);
        for iii = 1:7
            isVisibleTarget(iii) = ~isequal(trial(trialNum).eye.tar(iii).color_dim, [0 0 0]); % check if equals to background color
        end
        isVisibleTarget = sum(isVisibleTarget(2:7)); % check if one of the peripheral targets (2-7) isn't black; if this 1 - then peripheral target visible; if 0, then invisible
        
        if trial(trialNum).completed
            if trial(trialNum).rewarded && isVisibleTarget % rewarded, peripheral on - HIT
                trial(trialNum).SDT_trial_type = 1; % assign 1 to HIT trials
            elseif ~trial(trialNum).rewarded && isVisibleTarget % non-rewarded, peripheral on - MISS
                trial(trialNum).SDT_trial_type = 2; % assign 2 to MISS trials
            elseif ~trial(trialNum).rewarded && ~isVisibleTarget % non-rewarded, peripheral off - FALSE ALARM
                trial(trialNum).SDT_trial_type = 3; % assign 3 to FALSE ALARM
            elseif trial(trialNum).rewarded && ~isVisibleTarget % rewarded, peripheral off - CORRECT REJECTION
                trial(trialNum).SDT_trial_type = 4; % assign 4 to CORRECT REJECTION
            end
        else
            trial(trialNum).SDT_trial_type = 0; % don't analyze with STD
        end
    end
    
    save([combined_files(blockNum).name(1:end-4) '_lnv.mat'], 'trial')
    
end
