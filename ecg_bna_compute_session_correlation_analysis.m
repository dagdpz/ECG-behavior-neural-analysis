function ecg_bna_compute_session_correlation_analysis(trials,population,Rpeaks,cfg)
% Function to perform the correlation analysis between RR durations and
% corresponding FR within those cycles.
% This analysis had two important considerations. To be able to compare the
% results of shifting data related to one another between task and rest, we
% made spike and ECG data continuous and omitted invalid R-peaks of invalid
% spike periods. Took all trials, not only 

basepath_to_save=[cfg.SPK_root_results_fldr filesep 'correlation_analysis'];
if ~exist(basepath_to_save,'dir')
    mkdir(basepath_to_save);
end

%% load list of selected units - won't process all the crap
if cfg.spk.apply_exclusion_criteria
    
    unit_list = load([cfg.SPK_root_results_fldr filesep 'unit_lists' filesep cfg.spk.unit_list]);
    unitList = unique(unit_list.unit_ids);
    
    % figure out which units take from this session
    selected_this_session = ismember({population.unit_ID}, unitList);
    population = population(selected_this_session);
    
end

Rblocks=[Rpeaks.block];

for unitNum = 1:length(population)
    
    %% get to processing
    disp(['Processing unit ' num2str(unitNum) ' out of ' num2str(length(population))])
    
    pop=population(unitNum);
    
    T=ph_get_unit_trials(pop,trials);
    
    %% Make sure we only take overlapping blocks
    blocks_unit=unique([pop.block]);
    blocks=intersect(blocks_unit,Rblocks);
    b=ismember(Rblocks,blocks);
    
    %% preallocate 'data' structure
    data.unitId            = pop.unit_ID;
    data.target            = pop.target;
    data.channel           = pop.channel;
    data.unit              = pop.block_unit{2,1};
    data.quantSNR          = pop.avg_SNR;
    data.Single_rating     = pop.avg_single_rating;
    data.stability_rating  = pop.avg_stability;
    data.thresholds_microV = single([NaN; NaN; NaN; NaN]);
    data.FR                = single(mean(pop.FR_average));
    data.cc_lag_list       = cfg.spk.lag_list;
    data.criteria          = pop.criteria;
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        data.(L).meanHR_bpm                     = single(NaN);
        data.(L).medianHR_bpm                   = single(NaN);
        data.(L).stdHR_bpm                      = single(NaN);
        data.(L).SDNN_ms                        = single(NaN);
        data.(L).cc_PSTH_feature                = single(nan(1,4));
        data.(L).pp_PSTH_feature                = single(nan(1,4));
        data.(L).pperm_PSTH_feature             = single(nan(1,4)); % permuted p for the correlation between each feature and phase PSTH
        data.(L).timeRRstart                    = single(nan(1,1));
        data.(L).FRbyRR_Hz                      = single(nan(1,1));
        data.(L).cycleDurations_s               = single(nan(1,1));
        data.(L).pearson_r                      = single(nan(length(cfg.spk.lag_list), 1));
        data.(L).pearson_p                      = single(nan(length(cfg.spk.lag_list), 1));
        data.(L).permuted_p                     = single(nan(length(cfg.spk.lag_list), 1));
    end
    
    % find the corresponding WC file and load it
    chNum = data.channel;
    blkNum = unique([pop.block]);
    WCfile = ph_figure_out_waveclus_file_by_channel_and_blocks(chNum, blkNum, cfg.Input_WC);
    WC = load(WCfile, 'thr');
    if length(WC.thr) == 4
        data.thresholds_microV = 10^6 * WC.thr;
        data.thresholds_microV(3:4) = -1*data.thresholds_microV(3:4);
    end
    clear WC
    
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        
        % 1. take all the trials
        % 2. figure out spike and RR data for those trials
        % 3. figure out indices for invalid (wrong block, wrong condition, not completed)
        % 4. replace data belonging to wrong trials with NANs (both for spike and RR)
        % 5. replace data from RR intervals belonging to ITIs with NANs
        
        % What to do with inter-block intervals? - As one doesn't know how
        % many RR-intervals are there between blocks, let's replace those
        % with 12 RR slots - those will be nans and the corresponding FRs
        % as well
        
        %% get condition AND valid block trials only
        % create special trial configuration
        CT = ecg_bna_get_condition_trials(T, cfg.condition(c));
        tr=ismember([T.block],blocks) & CT;
        if sum(tr)<=1 || (~isfield(Rpeaks, 'RPEAK_ts_insp') && cfg.process_Rpeaks_inhalation_exhalation) || (~isfield(Rpeaks, 'RPEAK_ts_exp') && cfg.process_Rpeaks_inhalation_exhalation) % do calculations only if number of trials > 1
            continue
        end
        popcell=pop.trial(tr);
        trcell=T(tr);
        
        % 0. Prepare data variables
        states_onset               = {trcell.states_onset};
        states                     = {trcell.states};
        TDT_ECG1_t0_from_rec_start = {trcell.TDT_ECG1_t0_from_rec_start};
        block_nums                 = {trcell.block};
        state2_times               = cellfun(@(x,y) x(y == 2), states_onset, states, 'Uniformoutput', false); % trial starts = state 2
        state90_times              = cellfun(@(x,y) x(y == 90), states_onset, states, 'Uniformoutput', false); % trial ends = state 90
        % compute RR-intervals
        valid_RRinterval_ends      = single([Rpeaks(b).(['RPEAK_ts' cfg.condition(c).Rpeak_field])]);
        valid_RRinterval_starts    = single(valid_RRinterval_ends - [Rpeaks(b).(['RPEAK_dur' cfg.condition(c).Rpeak_field])]);
        % 0. figure out RR-intervals lying within trials
        trial_starts_one_stream    = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, state2_times, TDT_ECG1_t0_from_rec_start, block_nums);
        trial_ends_one_stream      = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, state90_times, TDT_ECG1_t0_from_rec_start, block_nums);
        
        RR_within_trial_idx = false(length(valid_RRinterval_starts),1);
        for RRnum = 1:length(RR_within_trial_idx)
            RR_within_trial_idx(RRnum) = ...
                any(valid_RRinterval_starts(RRnum)>trial_starts_one_stream & ...
                valid_RRinterval_ends(RRnum)<trial_ends_one_stream);
        end
        
        % 1. take arrival times and the corresponding waveforms
        AT = {popcell.arrival_times};
        WF = {popcell.waveforms};
        % 2. choose only those that happen after MP-state 1
        idx_after_state1 = cellfun(@(x,y) x>y, AT, state2_times, 'Uniformoutput', false);
        % 3. add TDT_ECG1_t0_from_rec_start and Rpeak block offset to spike times
        AT_one_stream_cell = cellfun(@(x,y,z,a) x(y)+z+Rpeaks([Rpeaks.block] == a).offset, AT, idx_after_state1, TDT_ECG1_t0_from_rec_start, block_nums, 'Uniformoutput', false);
        WF_one_stream_cell = cellfun(@(x,y) x(y,:), WF, idx_after_state1, 'Uniformoutput', false);
        % 4. merge all spike times and waveforms together
        AT_one_stream = cat(1, AT_one_stream_cell{:});
        WF_one_stream = cat(1, WF_one_stream_cell{:});
        if size(WF_one_stream,1) < 10
            continue
        end
        % 5. calculate heart cycle phase where individual spikes ended up
        [eventPhases, eventsTaken, cycleNums_withSpikes] = DAG_eventPhase(valid_RRinterval_starts, valid_RRinterval_ends, AT_one_stream); % all data
        
        % check the number of spikes left after computing phases
        if length(eventPhases) < 3 || length(valid_RRinterval_starts) < 3
            continue
        end
        
        % II. Compute unit firing rate per RR-interval
        data.(L).timeRRstart           = valid_RRinterval_starts;
        [FRbyRR_Hz, cycleDurations_s] = ...
            computeFRperCycle(valid_RRinterval_starts, valid_RRinterval_ends, AT_one_stream);
        FRbyRR_Hz(~RR_within_trial_idx)      = NaN;
        cycleDurations_s(~RR_within_trial_idx) = NaN;
        data.(L).FRbyRR_Hz = FRbyRR_Hz;
        data.(L).cycleDurations_s = cycleDurations_s;
        % compute correlation with different lag
        for lagNum = 1:length(cfg.spk.lag_list)
            % create data variables
            fr_hz = data.(L).FRbyRR_Hz;
            rr_s  = circshift(data.(L).cycleDurations_s, cfg.spk.lag_list(lagNum));
            
            % compute correlation coefficient
            [temp_r, temp_p] = corrcoef(fr_hz, rr_s, 'Rows', 'pairwise');
            data.(L).pearson_r(lagNum) = temp_r(2,1);
            data.(L).pearson_p(lagNum) = temp_p(2,1);
            data.(L).permuted_p(lagNum) = mult_comp_perm_corr(fr_hz, rr_s, cfg.spk.n_permutations, 0, 0.05, 'linear', 0);
            clear temp_r temp_p
%             xcorr
        end
    end
    save([basepath_to_save filesep data.unitId '_' data.target '__correlation.mat'], 'data', '-v7.3')
    clear data

end
end

function [FRbyRR_Hz, cycleDurations] = computeFRperCycle(interval_starts, interval_ends, eventTimes)
% make inputs vertical vectors
interval_starts = interval_starts';
interval_ends = interval_ends';
eventTimes = eventTimes';

% Calculate ECG cycle durations
cycleDurations = interval_ends - interval_starts;

% cycleNums = arrayfun(@(x) sum(x > interval_starts & x < interval_ends), eventTimes, 'UniformOutput', false);
spikeCounts = arrayfun(@(x,y) sum(eventTimes > x & eventTimes < y), interval_starts, interval_ends);
FRbyRR_Hz = spikeCounts ./ cycleDurations;

% figure,
% scatter(cycleDurations, FRbyRR_Hz)
end