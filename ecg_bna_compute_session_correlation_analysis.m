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
        state98_times              = cellfun(@(x,y) x(y == 98), states_onset, states, 'Uniformoutput', false); % trial ends = state 90
        % compute RR-intervals
        valid_RRinterval_ends      = single([Rpeaks(b).(['RPEAK_ts' cfg.condition(c).Rpeak_field])]);
        valid_RRinterval_starts    = single(valid_RRinterval_ends - [Rpeaks(b).(['RPEAK_dur' cfg.condition(c).Rpeak_field])]);
        % 0. figure out RR-intervals lying within trials
        trial_starts_one_stream    = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, state2_times, TDT_ECG1_t0_from_rec_start, block_nums);
        trial_ends_one_stream      = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, state98_times, TDT_ECG1_t0_from_rec_start, block_nums);
        
        RR_within_trial_idx = false(length(valid_RRinterval_starts),1);
        for RRnum = 1:length(RR_within_trial_idx)
            RR_within_trial_idx(RRnum) = ...
                any(valid_RRinterval_starts(RRnum)>trial_starts_one_stream & ...
                valid_RRinterval_ends(RRnum)<trial_ends_one_stream);
        end
        valid_RRinterval_ends      = valid_RRinterval_ends(RR_within_trial_idx); % get rid of RRs beyond the current set of trials
        valid_RRinterval_starts    = valid_RRinterval_starts(RR_within_trial_idx);
        
        % compute parameters of heart activity
        data.(L).meanHR_bpm   = mean(60 ./ (valid_RRinterval_ends - valid_RRinterval_starts));
        data.(L).medianHR_bpm = median(60 ./ (valid_RRinterval_ends - valid_RRinterval_starts));
        data.(L).stdHR_bpm    = std(60 ./ (valid_RRinterval_ends - valid_RRinterval_starts));
        data.(L).SDNN_ms      = std(1000 * (valid_RRinterval_ends - valid_RRinterval_starts));
        
        % implement median split to heart-cycle durations
        RRs     = valid_RRinterval_ends - valid_RRinterval_starts;
        M_IBI   = median(valid_RRinterval_ends - valid_RRinterval_starts);
        
        data.(L).IBI_median           = M_IBI;
        
        data.(L).lowIBI_timeRRstart   = valid_RRinterval_starts(RRs <= M_IBI);
        data.(L).lowIBI_timeRRend     = valid_RRinterval_ends(RRs <= M_IBI);
        data.(L).lowIBI_starts        = RRs(RRs <= M_IBI);
        data.(L).lowIBI_meanHR_bpm    = mean(60 ./ data.(L).lowIBI_starts);
        data.(L).lowIBI_medianHR_bpm  = median(60 ./ data.(L).lowIBI_starts);
        data.(L).lowIBI_stdHR_bpm     = std(60 ./ data.(L).lowIBI_starts);
        data.(L).lowIBI_SDNN_ms       = std(1000 * data.(L).lowIBI_starts);
        
        data.(L).highIBI_timeRRstart  = valid_RRinterval_starts(RRs > M_IBI);
        data.(L).highIBI_timeRRend    = valid_RRinterval_ends(RRs > M_IBI);
        data.(L).highIBI_starts       = RRs(RRs > M_IBI);
        data.(L).highIBI_meanHR_bpm   = mean(60 ./ data.(L).highIBI_starts);
        data.(L).highIBI_medianHR_bpm = median(60 ./ data.(L).highIBI_starts);
        data.(L).highIBI_stdHR_bpm    = std(60 ./ data.(L).highIBI_starts);
        data.(L).highIBI_SDNN_ms      = std(1000 * data.(L).highIBI_starts);
        
        
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
        [lowIBI_eventPhases, lowIBI_eventsTaken, lowIBI_cycleNums_withSpikes] = DAG_eventPhase(data.(L).lowIBI_timeRRstart, data.(L).lowIBI_timeRRend, AT_one_stream);
        [highIBI_eventPhases, highIBI_eventsTaken, highIBI_cycleNums_withSpikes] = DAG_eventPhase(data.(L).highIBI_timeRRstart, data.(L).highIBI_timeRRend, AT_one_stream);
        
        % check the number of spikes left after computing phases
        if length(eventPhases) < 3 || length(lowIBI_eventPhases) < 3 || length(highIBI_eventPhases) < 3 || ...
                length(valid_RRinterval_starts) < 3 || length(data.(L).lowIBI_timeRRstart) < 3 || length(data.(L).highIBI_timeRRstart) < 3
            continue
        end
        
        % 6. Put results for real data together
        data.(L).spike_phases_radians                 = eventPhases;
        data.(L).spike_phases_histogram               = hist(data.(L).spike_phases_radians, cfg.spk.phase_bin_centers); % compute overal phase histogram
        data.(L).spike_phases_histogram2              = hist3([eventPhases, cycleNums_withSpikes], 'ctrs', {cfg.spk.phase_bin_centers 1:length(valid_RRinterval_starts)});
        
        data.(L).lowIBI_spike_phases_radians          = lowIBI_eventPhases;
        data.(L).lowIBI_spike_phases_histogram        = hist(data.(L).lowIBI_spike_phases_radians, cfg.spk.phase_bin_centers);
        data.(L).lowIBI_spike_phases_histogram2       = hist3([lowIBI_eventPhases, lowIBI_cycleNums_withSpikes], 'ctrs', {cfg.spk.phase_bin_centers 1:length(data.(L).lowIBI_timeRRstart)});
        
        data.(L).highIBI_spike_phases_radians         = highIBI_eventPhases;
        data.(L).highIBI_spike_phases_histogram       = hist(data.(L).highIBI_spike_phases_radians, cfg.spk.phase_bin_centers);
        data.(L).highIBI_spike_phases_histogram2      = hist3([highIBI_eventPhases, highIBI_cycleNums_withSpikes], 'ctrs', {cfg.spk.phase_bin_centers 1:length(data.(L).highIBI_timeRRstart)});
        
        % Mosher's cosine fit the smoothed phase histogram
        [modIndex_phase_hist,phase_hist] = ...
            fitCardiacModulation(cfg.spk.phase_bin_centers, ...
            data.(L).spike_phases_histogram, {'PSTH'}, 0, [221]);
        
        data.(L).spike_phases_histogram_smoothed = phase_hist;
        data.(L).histogram_MI                    = modIndex_phase_hist(1);
        data.(L).histogram_p                     = modIndex_phase_hist(2);
        data.(L).histogram_phase                 = modIndex_phase_hist(3);
        data.(L).rsquared                        = modIndex_phase_hist(4);
        
        % LINEAR FITS
        % linear fit with 'fitlm' for all the data
        data.(L).linear              = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'linear');
        
        % linear fit for the lower IBI
        data.(L).lowIBI_linear       = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram2, length(data.(L).lowIBI_timeRRstart), 'linear');
        
        % linear fit for the higher IBI
        data.(L).highIBI_linear      = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).highIBI_spike_phases_histogram2, length(data.(L).highIBI_timeRRstart), 'linear');
        
        % COSINE FITS
        % cosine fit with 'fit' on all data points instead of smoothed mean
        data.(L).cosine              = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'cosine');
        
        % cosine fit for the lower IBI
        data.(L).lowIBI_cosine       = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram2, length(data.(L).lowIBI_timeRRstart), 'cosine');
        
        % cosine fit for the higher IBI
        data.(L).highIBI_cosine      = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).highIBI_spike_phases_histogram2, length(data.(L).highIBI_timeRRstart), 'cosine');
        
        % POSITIVE VON MISES FITS
        % fit all the data
        data.(L).vonMisesPos         = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'vonMises', 1);
        
        % pos. von Mises fits - low IBI
        data.(L).lowIBI_vonMisesPos  = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram2, length(data.(L).lowIBI_timeRRstart), 'vonMises', 1);
        
        % pos. von Mises fits - high IBI
        data.(L).highIBI_vonMisesPos = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).highIBI_spike_phases_histogram2, length(data.(L).highIBI_timeRRstart), 'vonMises', 1);
        
        % NEGATIVE VON MISES FITS
        % fit all the data
        data.(L).vonMisesNeg         = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'vonMises', -1);
        
        % neg. von Mises fits - low IBI
        data.(L).lowIBI_vonMisesNeg  = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).lowIBI_spike_phases_histogram2, length(data.(L).lowIBI_timeRRstart), 'vonMises', -1);
        
        % neg. von Mises fits - high IBI
        data.(L).highIBI_vonMisesNeg = ecg_bna_fit_neuronal_data(cfg, cfg.spk.phase_bin_centers, data.(L).highIBI_spike_phases_histogram2, length(data.(L).highIBI_timeRRstart), 'vonMises', -1);
        
        
        [~, ~, bin] = histcounts(data.(L).spike_phases_radians, cfg.spk.phase_bins);
        
        data.(L).waveforms_microvolts                 = 10^6 * WF_one_stream(eventsTaken,:);
        waveforms_upsampled                           = interpft(data.(L).waveforms_microvolts, 32*4, 2);
        data.(L).waveforms_upsampled_microvolts       = shift2peak(cfg.spk.wf_times_interp_ms, waveforms_upsampled);
        data.(L).waveforms_byBin_microvolts           = arrayfun(@(x) mean(data.(L).waveforms_upsampled_microvolts(bin == x,:),1), 1:cfg.spk.N_phase_bins, 'UniformOutput', false);
        data.(L).waveforms_byBin_microvolts           = cat(1,data.(L).waveforms_byBin_microvolts{:});
        % 7. Calculate spike features with Mosher's procedure
        tic
        sMetric = ...
            struct('extremAmp', cell(length(data.(L).spike_phases_radians),1), ...
            'widthHW', cell(length(data.(L).spike_phases_radians),1), ...
            'widthTP', cell(length(data.(L).spike_phases_radians),1), ...
            'repolTime', cell(length(data.(L).spike_phases_radians),1));
        parfor wfNum = 1:length(data.(L).spike_phases_radians)
            sMetric(wfNum)=spikeWaveMetrics(double(data.(L).waveforms_upsampled_microvolts(wfNum,:)), 37, cfg.spk.Fs*4, 0, [1 0 0 0]);
%             try %% spikeWaveMetrics missing! --> try catch is also VERY inefficient - more efficient is to skip excluded units
%                 sMetric(wfNum)=spikeWaveMetrics(double(data.(L).waveforms_upsampled_microvolts(wfNum,:)), 37, cfg.spk.Fs*4, 0, [1 1 0 0]); % 37 - index of th peak for updsampled data
%                 sMetric(wfNum).extremAmp   = sMetric(wfNum).extremAmp(1);
%                 sMetric(wfNum).widthHW     = sMetric(wfNum).widthHW(1);
%                 sMetric(wfNum).widthTP     = sMetric(wfNum).widthTP(1);
%                 sMetric(wfNum).repolTime   = sMetric(wfNum).repolTime(1);
%             catch ME
%                 sMetric(wfNum).extremAmp  = nan;
%                 sMetric(wfNum).widthHW    = nan;
%                 sMetric(wfNum).widthTP    = nan;
%                 sMetric(wfNum).repolTime  = nan;
%             end
        end
        toc
        % 8. put the resulting data together
        data.(L).AMP_microV               = [sMetric.extremAmp];
        data.(L).HW_ms                    = 10^3 * [sMetric.widthHW];
        data.(L).TPW_ms                   = 10^3 * [sMetric.widthTP];
        data.(L).REP_ms                   = 10^3 * [sMetric.repolTime];
        
        data.(L).AMP_microV_byBin         = arrayfun(@(x) nanmean(data.(L).AMP_microV(bin == x)), 1:cfg.spk.N_phase_bins); % mean by phase
        data.(L).HW_ms_byBin              = arrayfun(@(x) nanmean(data.(L).HW_ms(bin == x)), 1:cfg.spk.N_phase_bins);
        data.(L).TPW_ms_byBin             = arrayfun(@(x) nanmean(data.(L).TPW_ms(bin == x)), 1:cfg.spk.N_phase_bins);
        data.(L).REP_ms_byBin             = arrayfun(@(x) nanmean(data.(L).REP_ms(bin == x)), 1:cfg.spk.N_phase_bins);
        
        % 9. estimate significance with bootstrapping
        [data.(L).AMP_lowerPrctile_2_5, data.(L).AMP_upperPrctile_97_5, data.(L).AMP_reshuffled_avg] = ...
            compute_reshuffles(data.(L).AMP_microV, bin, cfg);
        [data.(L).HW_lowerPrctile_2_5, data.(L).HW_upperPrctile_97_5, data.(L).HW_reshuffled_avg] = ...
            compute_reshuffles(data.(L).HW_ms, bin, cfg);
        [data.(L).TPW_lowerPrctile_2_5, data.(L).TPW_upperPrctile_97_5, data.(L).TPW_reshuffled_avg] = ...
            compute_reshuffles(data.(L).TPW_ms, bin, cfg);
        [data.(L).REP_lowerPrctile_2_5, data.(L).REP_upperPrctile_97_5, data.(L).REP_reshuffled_avg] = ...
            compute_reshuffles(data.(L).REP_ms, bin, cfg);
        
        % 10. compute cosine fits and put those into the resulting
        % structure
        featureMatrix = ...
            [data.(L).AMP_microV_byBin; data.(L).HW_ms_byBin; ...
            data.(L).TPW_ms_byBin; data.(L).REP_ms_byBin];
        if size(featureMatrix,1) > 4
            featureMatrix = featureMatrix';
            if size(featureMatrix,1) ~= 4
                error('Dimensions of feature matrix aren''t suitable for the analysis')
            end
        end
        
        if length(data.(L).AMP_microV) - sum(isnan(data.(L).AMP_microV)) < 3
            continue
        end
        
        if sum(isnan(data.(L).AMP_microV)) == length(data.(L).AMP_microV) % if they're all nan
            data.(L).AMP_MI                   = nan(1,5);
        else
            coefs = cardioballistic_fit(data.(L).AMP_microV, eventPhases, cfg);
            data.(L).AMP_MI                   = coefs;
            clear coefs
        end
        
        if sum(isnan(data.(L).HW_ms)) == length(data.(L).HW_ms) % if they're all nan
            data.(L).HW_MI                    = nan(1,5);
        else
            coefs = cardioballistic_fit(data.(L).HW_ms, eventPhases, cfg);
            data.(L).HW_MI                    = coefs;
            clear coefs
        end
        
        if sum(isnan(data.(L).TPW_ms)) == length(data.(L).TPW_ms) % if they're all nan
            data.(L).TPW_MI                   = nan(1,5);
        else
            coefs = cardioballistic_fit(data.(L).TPW_ms, eventPhases, cfg);
            data.(L).TPW_MI                   = coefs;
            clear coefs
        end
        
        if sum(isnan(data.(L).REP_ms)) == length(data.(L).REP_ms) % if they're all nan
            data.(L).REP_MI                   = nan(1,5);
        else
            coefs = cardioballistic_fit(data.(L).REP_ms, eventPhases, cfg);
            data.(L).REP_MI                   = coefs;
            clear coefs
        end
        
        % store smoothed data for each measure
        if sum(isnan(data.(L).AMP_microV_byBin)) == length(data.(L).AMP_microV_byBin)
            data.(L).AMP_microV_byBin_smoothed    = nan(1,cfg.spk.N_phase_bins);
        else
            data.(L).AMP_microV_byBin_smoothed    = circ_smooth(data.(L).AMP_microV_byBin);
        end
        if sum(isnan(data.(L).HW_ms_byBin)) == length(data.(L).HW_ms_byBin)
            data.(L).HW_ms_byBin_smoothed         = nan(1,cfg.spk.N_phase_bins);
        else
            data.(L).HW_ms_byBin_smoothed         = circ_smooth(data.(L).HW_ms_byBin);
        end
        if sum(isnan(data.(L).TPW_ms_byBin)) == length(data.(L).TPW_ms_byBin)
            data.(L).TPW_ms_byBin_smoothed        = nan(1,cfg.spk.N_phase_bins);
        else
            data.(L).TPW_ms_byBin_smoothed        = circ_smooth(data.(L).TPW_ms_byBin);
        end
        if sum(isnan(data.(L).REP_ms_byBin)) == length(data.(L).REP_ms_byBin)
            data.(L).REP_ms_byBin_smoothed        = nan(1,cfg.spk.N_phase_bins);
        else
            data.(L).REP_ms_byBin_smoothed        = circ_smooth(data.(L).REP_ms_byBin);
        end
        
        %         data.(L).allCorr                      = allCorr;
        %         data.(L).allLinMod                    = allLinMod;
        
        [data.(L).AMP_max_consec_bins, data.(L).AMP_modulation_index] = ...
            significant_bins(data.(L).AMP_microV_byBin_smoothed', data.(L).AMP_lowerPrctile_2_5, data.(L).AMP_upperPrctile_97_5, data.(L).AMP_reshuffled_avg);
        [data.(L).HW_max_consec_bins, data.(L).HW_modulation_index] = ...
            significant_bins(data.(L).HW_ms_byBin_smoothed', data.(L).HW_lowerPrctile_2_5, data.(L).HW_upperPrctile_97_5, data.(L).HW_reshuffled_avg);
        [data.(L).TPW_max_consec_bins, data.(L).TPW_modulation_index] = ...
            significant_bins(data.(L).TPW_ms_byBin_smoothed', data.(L).TPW_lowerPrctile_2_5, data.(L).TPW_upperPrctile_97_5, data.(L).TPW_reshuffled_avg);
        [data.(L).REP_max_consec_bins, data.(L).REP_modulation_index] = ...
            significant_bins(data.(L).REP_ms_byBin_smoothed', data.(L).REP_lowerPrctile_2_5, data.(L).REP_upperPrctile_97_5, data.(L).REP_reshuffled_avg);
        
        % 11. estimate how much spike amplitude is far from the
        % thresholds
        % bin from the closest threshold first
        higher_threshold_detected = ...
            all(rmmissing(data.(L).AMP_microV) > data.thresholds_microV(1));
        lower_threshold_detected = ...
            all(rmmissing(data.(L).AMP_microV) > data.thresholds_microV(2));
        if higher_threshold_detected && lower_threshold_detected % both thresholds are lower than all the spike amplitudes
            % take the higher one
            data.(L).AMP_voltageBins   = data.thresholds_microV(1):500;
            data.(L).AMP_voltageBinned = histc(data.(L).AMP_microV, data.(L).AMP_voltageBins); % bin amplitudes
        elseif lower_threshold_detected && ~higher_threshold_detected
            % take the lower one
            data.(L).AMP_voltageBins   = data.thresholds_microV(2):500;
            data.(L).AMP_voltageBinned = histc(data.(L).AMP_microV, data.(L).AMP_voltageBins);
        elseif ~lower_threshold_detected && ~higher_threshold_detected
            error('It''s an unknown error: there are spikes below both thresholds')
        elseif ~lower_threshold_detected && higher_threshold_detected
            error('It''s an unknown error: there are spikes that are below only lower threshold')
        end
        
        first_nonzero_bin          = find(data.(L).AMP_voltageBinned,1,'first');
        
        data.(L).distance2thr = first_nonzero_bin; % basically, number of empty bins from the spike with the lowest amplitude and the closest threshold
        
        % 12. compute correlation between features phase dynamics and
        % spike dynamics
        [cc, pp] = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).AMP_microV_byBin_smoothed);
        data.(L).cc_PSTH_feature(1)    = cc(2,1); % for spike AMP
        data.(L).pp_PSTH_feature(1)    = pp(2,1);
        %         data.(L).pperm_PSTH_feature(1) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed, data.(L).AMP_microV_byBin_smoothed, cfg.spk.n_permutations, 0, 0.05, 'linear', 0);
        
        [cc, pp] = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).HW_ms_byBin_smoothed);
        data.(L).cc_PSTH_feature(2)    = cc(2,1); % for HW
        data.(L).pp_PSTH_feature(2)    = pp(2,1);
        %         data.(L).pperm_PSTH_feature(2) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed, data.(L).HW_ms_byBin_smoothed, cfg.spk.n_permutations, 0, 0.05, 'linear', 0);
        
        [cc, pp] = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).TPW_ms_byBin_smoothed);
        data.(L).cc_PSTH_feature(3)    = cc(2,1); % for spike TPW
        data.(L).pp_PSTH_feature(3)    = pp(2,1);
        %         data.(L).pperm_PSTH_feature(3) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed, data.(L).TPW_ms_byBin_smoothed, cfg.spk.n_permutations, 0, 0.05, 'linear', 0);
        
        [cc, pp] = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).REP_ms_byBin_smoothed);
        data.(L).cc_PSTH_feature(4)    = cc(2,1); % for spike REP
        data.(L).pp_PSTH_feature(4)    = pp(2,1);
        %         data.(L).pperm_PSTH_feature(4) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed, data.(L).REP_ms_byBin_smoothed, cfg.spk.n_permutations, 0, 0.05, 'linear', 0);
        
        % II. Compute unit firing rate per RR-interval
        data.(L).timeRRstart           = valid_RRinterval_starts;
        [data.(L).FRbyRR_Hz, ...
            data.(L).cycleDurations_s] = ...
            computeFRperCycle(valid_RRinterval_starts, valid_RRinterval_ends, AT_one_stream);
        % compute correlation with different lag
        for lagNum = 1:length(cfg.spk.lag_list)
            % create data variables
            fr_hz = data.(L).FRbyRR_Hz;
            rr_s  = circshift(data.(L).cycleDurations_s, cfg.spk.lag_list(lagNum));
            
            % compute correlation coefficient
            [temp_r, temp_p] = corrcoef(fr_hz, rr_s);
            data.(L).pearson_r(lagNum) = temp_r(2,1);
            data.(L).pearson_p(lagNum) = temp_p(2,1);
            data.(L).permuted_p(lagNum) = mult_comp_perm_corr(fr_hz, rr_s, cfg.spk.n_permutations, 0, 0.05, 'linear', 0);
            clear temp_r temp_p
            
        end
    end
    save([basepath_to_save filesep data.unitId '_' data.target '__correlation.mat'], 'data', '-v7.3')
    clear data

end

