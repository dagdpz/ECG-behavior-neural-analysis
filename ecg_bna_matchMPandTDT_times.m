clear all, close all

% filenames
% 'Y:\Data\Magnus_phys_combined_monkeypsych_TDT\20220921\Magcombined2022-09-21_10_block_01.mat'
% 'Y:\Data\Magnus_phys_combined_monkeypsych_TDT\20230518\Magcombined2023-05-18_02_block_01.mat'
% 'Y:\Data\Magnus_phys_combined_monkeypsych_TDT\20230518\Magcombined2023-06-15_04_block_03.mat'
% 'Y:\Data\Magnus_phys_combined_monkeypsych_TDT\20230518\Magcombined2023-06-15_06_block_05.mat'

load('Y:\Data\Magnus_phys_combined_monkeypsych_TDT\20230518\Magcombined2023-05-18_02_block_01.mat')
load('Y:\Data\BodySignals\ECG\Magnus\20230518\20230518_ecg.mat')

for ii = 1:length(trial)
    
    trial(ii).TDT_state_onsets = trial(ii).TDT_state_onsets';
    trial(ii).TDT_state_onsets = trial(ii).TDT_state_onsets + trial(ii).states_onset(2);
    
    if ii < 11
       txt_MP = mat2cell(trial(ii).states, 1, ones(length(trial(ii).states),1));
       txt_MP = cellfun(@num2str, txt_MP, 'UniformOutput', false);
       txt_TDT = mat2cell(trial(ii).TDT_states, ones(length(trial(ii).TDT_states),1), 1);
       txt_TDT = cellfun(@num2str, txt_TDT, 'UniformOutput', false);
        
       figure,
       set(gcf, 'Position', [681 559 1197 420])
       stem(trial(ii).states_onset, ones(length(trial(ii).states_onset), 1), 'r')
       text(trial(ii).states_onset, ones(length(trial(ii).states_onset), 1)+0.05, txt_MP, 'Color','r')
       hold on
       stem(trial(ii).TDT_state_onsets, 1.2*ones(length(trial(ii).TDT_state_onsets), 1), 'b')
       text(trial(ii).TDT_state_onsets, 1.2*ones(length(trial(ii).TDT_state_onsets), 1)+0.05, txt_TDT, 'Color','b')
       title(['Magnus ' num2str(trial(ii).TDT_session) '; Trial ' num2str(trial(ii).trial_number(1))])
       ylim([0 1.3])
    end
    
    isVisibleTarget = nan(7,1);
    for iii = 1:7
        isVisibleTarget(iii) = ~isequal(trial(ii).eye.tar(iii).color_dim, [0 0 0]); % check if equals to background color
    end
    isVisibleTarget = sum(isVisibleTarget(2:7)); % check if one of the peripheral targets (2-7) isn't black; if this 1 - then peripheral target visible; if 0, then invisible
    
    if trial(ii).completed
        if trial(ii).rewarded && isVisibleTarget % rewarded, peripheral on - HIT
            trial(ii).SDT_trial_type = 1; % assign 1 to HIT trials
        elseif ~trial(ii).rewarded && isVisibleTarget % non-rewarded, peripheral on - MISS
            trial(ii).SDT_trial_type = 2; % assign 2 to MISS trials
        elseif ~trial(ii).rewarded && ~isVisibleTarget % non-rewarded, peripheral off - FALSE ALARM
            trial(ii).SDT_trial_type = 3; % assign 3 to FALSE ALARM
        elseif trial(ii).rewarded && ~isVisibleTarget % rewarded, peripheral off - CORRECT REJECTION
            trial(ii).SDT_trial_type = 4; % assign 4 to CORRECT REJECTION
        end
    else
        trial(ii).SDT_trial_type = 0; % don't analyze with STD
    end
    
end

TDT_state_onsets = [trial.TDT_state_onsets];
MP_state_onsets = [trial.states_onset];
Rpeak_times = out(1).Rpeak_t;

figure,
stem(TDT_state_onsets, 2*ones(length(TDT_state_onsets), 1), 'g')
hold on
stem(MP_state_onsets, ones(length(MP_state_onsets), 1), 'k')
stem(Rpeak_times, 0.5*ones(length(Rpeak_times), 1), 'm')

save('Magcombined2022-09-21_10_block_01.mat', 'trial', 'task', 'SETTINGS', 'First_trial_INI', 'preprocessing_settings')
