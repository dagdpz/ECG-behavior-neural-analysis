function ecg_bna_compute_PEPH(trials,population,Triggers,Blockoffsets,cfg)

basepath_to_save=[cfg.SPK_root_results_fldr filesep 'cardioballistic'];
if ~exist(basepath_to_save,'dir')
    mkdir(basepath_to_save);
end

%% load list of selected units - won't process all the crap
%% How can we ever know this in advance??
% if cfg.spk.apply_exclusion_criteria
%     
%     unit_list = load([cfg.SPK_root_results_fldr filesep 'unit_lists' filesep cfg.spk.unit_list]);
%     unitList = unique(unit_list.unit_ids);
%     
%     % figure out which units take from this session
%     selected_this_session = ismember({population.unit_ID}, unitList);
%     population = population(selected_this_session);
%     
% end

% Rblocks=[Triggers.block];

for u = 1:numel(population)    
    pop=population(u);
    %% get to processing
    disp(['Processing ' pop.unit_ID ' ' num2str(u) '/' num2str(numel(population))]);
    T=ph_get_unit_trials(pop,trials);
    
%     %% Make sure we only take overlapping blocks
%     blocks_unit=unique([pop.block]);
%     blocks=intersect(blocks_unit,Rblocks);
%     b=ismember(Rblocks,blocks);
    
    %% preallocate 'data' structure
    data.unitId            = pop.unit_ID;
    data.target            = pop.target;
    data.quantSNR          = pop.avg_SNR;
    data.Single_rating     = pop.avg_single_rating;
    data.stability_rating  = pop.avg_stability;
    
%     data.channel           = pop.channel;
%     data.unit              = pop.block_unit{2,1};
%     data.thresholds_microV = single([NaN; NaN; NaN; NaN]);
%     data.FR                = single(mean(pop.FR_average));
%     data.criteria          = pop.criteria;


    % find the corresponding WC file and load it
    chNum = pop.channel;
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
        
        %% get condition AND valid block trials only - check, but this should really remove invalid blocks already
        tr = ecg_bna_get_condition_trials(T, cfg.condition(c));        
        data.(L).NrTrials=sum(tr);
        
%         CT = ecg_bna_get_condition_trials(T, cfg.condition(c));
%         tr=ismember([T.block],blocks) & CT;
        if sum(tr)<=1 
            continue
        end
        
        popcell=pop.trial(tr);
        trcell=T(tr);
        [O]=ecg_bna_synchronize_with_trigger(Blockoffsets,trcell,popcell);
        
        %% compute spike density as one continuous vector across all concatenated trials 
        % (hmmm there migth be a problem with interleaved trial types here)
        PSTH_time=O.trial_starts(1):cfg.spk.PSTH_binwidth:O.trial_ends(end);
        %[SD_stream, RAST]=ecg_bna_spike_density(O.AT,PSTH_time,cfg.spk);
%         
%         % 0. Prepare data variables
%         states_onset               = {trcell.states_onset};
%         states                     = {trcell.states};
%         TDT_ECG1_t0_from_rec_start = {trcell.TDT_ECG1_t0_from_rec_start};
%         block_nums                 = {trcell.block};
%         state2_times               = cellfun(@(x,y) x(y == 2), states_onset, states, 'Uniformoutput', false); % trial starts = state 2
%         state90_times              = cellfun(@(x,y) x(y == 90), states_onset, states, 'Uniformoutput', false); % trial ends = state 90
%         


        % 0. figure out RR-intervals lying within trials
        trial_starts_one_stream    = O.trial_starts;
        trial_ends_one_stream      = O.trial_ends;

           % E=events{e};
            E='R';
            
        % compute RR-intervals
        valid_RRinterval_ends      = Triggers.(E).ts;
        valid_RRinterval_starts    = Triggers.(E).ts  - Triggers.(E).intervals;
        RR_within_trial_idx=any(bsxfun(@ge,valid_RRinterval_starts',trial_starts_one_stream) & bsxfun(@le,valid_RRinterval_ends',trial_ends_one_stream),2)';
       
        
        
starts = valid_RRinterval_starts(RR_within_trial_idx);
ends = valid_RRinterval_ends(RR_within_trial_idx);
eventTimes = O.AT;

% Calculate ECG cycle durations
cycleDurations = interval_ends - interval_starts;

% Normalize event times to the respective ECG cycle durations
[events2include,cycleNums_withSpikes] = find(eventTimes > interval_starts' & ...
    eventTimes < interval_ends');
eventTimes = eventTimes(events2include); % include only events that landed within any interval
eventTimesNorm = (eventTimes - interval_starts(cycleNums_withSpikes)) ./ cycleDurations(cycleNums_withSpikes);
eventPhases = 2*pi*eventTimesNorm;
eventsTaken = single(find(events2include))';
        
        
        [eventPhases, eventsTaken, cycleNums_withSpikes] = ...
            DAG_eventPhase(valid_RRinterval_starts, valid_RRinterval_ends, AT_one_stream); % all data
        
        % compute shuffled RR-intervals
        shuffled_RRinterval_ends   = Triggers.(E).shuffled_ts;
        shuffled_RRinterval_starts = Triggers.(E).shuffled_ts  - Triggers.(E).shuffled_intervals;
        for s = 1:size(shuffled_RRinterval_starts,1)
            shuffledRR_within_trial_idx(s,:)=any(bsxfun(@ge,shuffled_RRinterval_starts(s,:)',trial_starts_one_stream) & ...
                bsxfun(@le,shuffled_RRinterval_ends(s,:)',trial_ends_one_stream),2)';
        end
        
        
%               @le             Less than or equal
%               @gt             Greater than
%               @ge             Greater than or equal
        
        
        % Create a logical array by vectorized comparison
         %RR_within_trial_idx = any(valid_RRinterval_starts' > trial_starts_one_stream & valid_RRinterval_ends' < trial_ends_one_stream, 2);
        
%         % Use logical indexing to filter valid RR intervals
%         valid_RRinterval_starts = valid_RRinterval_starts(RR_within_trial_idx);
%         valid_RRinterval_ends   = valid_RRinterval_ends(RR_within_trial_idx);
        
        % reshuffled - get rid of RRs beyond the current set of trials
      
        
        % compute parameters of heart activity - disregard trial
        % limitations here (?)
        data.(L).meanHR_bpm   = mean(60 ./ (valid_RRinterval_ends - valid_RRinterval_starts));
        data.(L).medianHR_bpm = median(60 ./ (valid_RRinterval_ends - valid_RRinterval_starts));
        data.(L).stdHR_bpm    = std(60 ./ (valid_RRinterval_ends - valid_RRinterval_starts));
        data.(L).SDNN_ms      = std(1000 * (valid_RRinterval_ends - valid_RRinterval_starts));
        
%         % implement median split to heart-cycle durations
%         RRs   = valid_RRinterval_ends - valid_RRinterval_starts;
%         M_IBI = median(valid_RRinterval_ends - valid_RRinterval_starts);
%         
%         data.(L).IBI_median = M_IBI;
        
        
%         % 1. take arrival times and the corresponding waveforms
%         AT = {popcell.arrival_times};
%         WF = {popcell.waveforms};
%         % 2. choose only those that happen after MP-state 1
%         idx_after_state1 = cellfun(@(x,y) x>y, AT, state2_times, 'Uniformoutput', false);
%         % 3. add TDT_ECG1_t0_from_rec_start and Rpeak block offset to spike times
%         AT_one_stream_cell = cellfun(@(x,y,z,a) x(y)+z+Triggers([Triggers.block] == a).offset, AT, idx_after_state1, TDT_ECG1_t0_from_rec_start, block_nums, 'Uniformoutput', false);
%         WF_one_stream_cell = cellfun(@(x,y) x(y,:), WF, idx_after_state1, 'Uniformoutput', false);
%         % 4. merge all spike times and waveforms together
%         AT_one_stream = cat(1, AT_one_stream_cell{:});
%         WF_one_stream = cat(1, WF_one_stream_cell{:});
%         if size(WF_one_stream,1) < 10
%             continue
%         end
        % 5. calculate heart cycle phase where individual spikes ended up
        
        
        % check the number of spikes left after computing phases
        if length(eventPhases) < 3 || length(lowIBI_eventPhases) < 3 || length(highIBI_eventPhases) < 3 || ...
                length(valid_RRinterval_starts) < 100 || length(data.(L).lowIBI_timeRRstart) < 3 || length(data.(L).highIBI_timeRRstart) < 3
            continue
        end
        
        % for reshuffled data - calculate heart cycle phase where individual spikes ended up
        shuffled_SDF         = nan(cfg.time.n_permutations,cfg.phase.N_phase_bins);
        shuffled_SDF_lowIBI  = nan(cfg.time.n_permutations,cfg.phase.N_phase_bins);
        shuffled_SDF_highIBI = nan(cfg.time.n_permutations,cfg.phase.N_phase_bins);
        
        tic
%         medianIBI = data.(L).IBI_median;
        phase_bin_centers = cfg.phase.phase_bin_centers;
        parfor s = 1:length(shuffled_RRinterval_starts)
            % all data
            [shuffled_eventPhases, ~, shuffled_cycleNums_withSpikes] = ...
                DAG_eventPhase(shuffled_RRinterval_starts{s}, shuffled_RRinterval_ends{s}, AT_one_stream);
            
            % all data
            shuffled_histogram = ...
                hist3([shuffled_eventPhases, shuffled_cycleNums_withSpikes], ...
                'ctrs', {phase_bin_centers 1:length(shuffled_RRinterval_starts{s})});
            
            shuffled_SDF(s,:) = mean(average_smooth_data(shuffled_histogram,cfg));
        end
        toc
        
        % 6. Put results for real data together
        data.(L).spike_phases_radians    = eventPhases;
        data.(L).spike_phases_histogram  = hist(data.(L).spike_phases_radians, cfg.phase.phase_bin_centers) / length(valid_RRinterval_starts) / length(cfg.phase.phase_bin_centers); % compute overal phase histogram
        data.(L).spike_phases_histogram2 = hist3([eventPhases, cycleNums_withSpikes], 'ctrs', {cfg.phase.phase_bin_centers 1:length(valid_RRinterval_starts)});
       
        % all data
        real_data = average_smooth_data(data.(L).spike_phases_histogram2,cfg);
        SD        = ecg_bna_do_statistics(real_data,shuffled_SDF,1:cfg.phase.N_phase_bins);
        
        data.(L).SD          = SD.SD_mean;
        data.(L).SD_STD      = SD.SD_STD;
        data.(L).SD_SEM      = SD.SD_SEM ;
        data.(L).SDP         = SD.SDPmean ;
        data.(L).SDPCL       = SD.SDPconf(1,:) ;
        data.(L).SDPCu       = SD.SDPconf(2,:) ;
        data.(L).sig_all     = SD.sig_all;
        data.(L).sig         = SD.sig;
        data.(L).sig_FR_diff = SD.sig_FR_diff;
        data.(L).sig_time    = SD.sig_time;
        data.(L).sig_n_bins  = SD.sig_n_bins;
        data.(L).sig_sign    = SD.sig_sign;
        data.(L).SDsubstractedSDP            = data.(L).SD - data.(L).SDP; % spikes/s, difference between mean and jittered data
        data.(L).SDsubstractedSDP_normalized = data.(L).SDsubstractedSDP ./ data.(L).SDP *100; % percent signal change
        data.(L).FR_ModIndex_SubtrSDP        = max(data.(L).SDsubstractedSDP) - min(data.(L).SDsubstractedSDP); % difference between max and min FR
        data.(L).FR_ModIndex_PcS             = max(data.(L).SDsubstractedSDP_normalized) - min(data.(L).SDsubstractedSDP_normalized); % difference between max and min % signal change
        
       
        % Mosher's cosine fit the smoothed phase histogram
        [modIndex_phase_hist,phase_hist] = ...
            fitCardiacModulation(cfg.phase.phase_bin_centers, ...
            data.(L).spike_phases_histogram, {'PSTH'}, 0, [221]);
        
        data.(L).spike_phases_histogram_smoothed = phase_hist;
        data.(L).histogram_MI                    = modIndex_phase_hist(1);
        data.(L).histogram_p                     = modIndex_phase_hist(2);
        data.(L).histogram_phase                 = modIndex_phase_hist(3);
        data.(L).rsquared                        = modIndex_phase_hist(4);
        
        %% LINEAR FITS
        % linear fit with 'fitlm' for all the data
%         data.(L).linear                 = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'linear');
        
        % linear fit with 'fitlm' for smoothed data
%         data.(L).linear_smoothed        = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram_smoothed, length(valid_RRinterval_starts), 'linear');
     
        %% COSINE FITS
        % cosine fit with 'fit' on all data points instead of smoothed mean
%         data.(L).cosine              = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'cosine');
        
        % cosine fit with 'fit' for smoothed data
%         data.(L).cosine_smoothed     = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram_smoothed, length(valid_RRinterval_starts), 'cosine');
     
        %% POSITIVE VON MISES FITS
        % fit all the data
%         data.(L).vonMisesPos          = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'vonMises', 1);
        
        % fit all the data - smoothed
%         data.(L).vonMisesPos_smoothed = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram_smoothed, length(valid_RRinterval_starts), 'vonMises', 1);
    
        %% NEGATIVE VON MISES FITS
        % fit all the data
%         data.(L).vonMisesNeg          = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram2, length(valid_RRinterval_starts), 'vonMises', -1);
        
        % fit all the data - smoothed
%         data.(L).vonMisesNeg_smoothed = ecg_bna_fit_neuronal_data(cfg, cfg.phase.phase_bin_centers, data.(L).spike_phases_histogram_smoothed, length(valid_RRinterval_starts), 'vonMises', -1);
      
        [~, ~, bin] = histcounts(data.(L).spike_phases_radians, cfg.phase.phase_bins);
        
        data.(L).waveforms_microvolts                 = 10^6 * WF_one_stream(eventsTaken,:);
        waveforms_upsampled                           = interpft(data.(L).waveforms_microvolts, 32*4, 2);
        data.(L).waveforms_upsampled_microvolts       = shift2peak(cfg.phase.wf_times_interp_ms, waveforms_upsampled);
        data.(L).waveforms_byBin_microvolts           = arrayfun(@(x) mean(data.(L).waveforms_upsampled_microvolts(bin == x,:),1), 1:cfg.phase.N_phase_bins, 'UniformOutput', false);
        data.(L).waveforms_byBin_microvolts           = cat(1,data.(L).waveforms_byBin_microvolts{:});
        % 7. Calculate spike features with Mosher's procedure
        tic
        sMetric = ...
            struct('extremAmp', cell(length(data.(L).spike_phases_radians),1), ...
            'widthHW', cell(length(data.(L).spike_phases_radians),1), ...
            'widthTP', cell(length(data.(L).spike_phases_radians),1), ...
            'repolTime', cell(length(data.(L).spike_phases_radians),1));
        parfor wfNum = 1:length(data.(L).spike_phases_radians)
            sMetric(wfNum)=spikeWaveMetrics(double(data.(L).waveforms_upsampled_microvolts(wfNum,:)), 37, cfg.phase.Fs*4, 0, [1 0 0 0]); % 37 - index of th peak for updsampled data
        end
        toc
        % 8. put the resulting data together
        data.(L).AMP_microV               = [sMetric.extremAmp];
        data.(L).HW_ms                    = 10^3 * [sMetric.widthHW];
        data.(L).TPW_ms                   = 10^3 * [sMetric.widthTP];
        data.(L).REP_ms                   = 10^3 * [sMetric.repolTime];
        
        data.(L).AMP_microV_byBin         = arrayfun(@(x) nanmean(data.(L).AMP_microV(bin == x)), 1:cfg.phase.N_phase_bins); % mean by phase
        data.(L).HW_ms_byBin              = arrayfun(@(x) nanmean(data.(L).HW_ms(bin == x)), 1:cfg.phase.N_phase_bins);
        data.(L).TPW_ms_byBin             = arrayfun(@(x) nanmean(data.(L).TPW_ms(bin == x)), 1:cfg.phase.N_phase_bins);
        data.(L).REP_ms_byBin             = arrayfun(@(x) nanmean(data.(L).REP_ms(bin == x)), 1:cfg.phase.N_phase_bins);
        
        % 9. estimate significance with bootstrapping
        [data.(L).AMP_lowerPrctile_2_5, data.(L).AMP_upperPrctile_97_5, data.(L).AMP_reshuffled_avg] = ...
            compute_reshuffles(data.(L).AMP_microV, bin, cfg.phase);
        [data.(L).HW_lowerPrctile_2_5, data.(L).HW_upperPrctile_97_5, data.(L).HW_reshuffled_avg] = ...
            compute_reshuffles(data.(L).HW_ms, bin, cfg.phase);
        [data.(L).TPW_lowerPrctile_2_5, data.(L).TPW_upperPrctile_97_5, data.(L).TPW_reshuffled_avg] = ...
            compute_reshuffles(data.(L).TPW_ms, bin, cfg.phase);
        [data.(L).REP_lowerPrctile_2_5, data.(L).REP_upperPrctile_97_5, data.(L).REP_reshuffled_avg] = ...
            compute_reshuffles(data.(L).REP_ms, bin, cfg.phase);
        
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
        
        % use Mosher's cosine fitting procedure with binned data
        modIndex        = fitCardiacModulation(cfg.phase.phase_bin_centers, featureMatrix, {'AMP', 'HW', 'TPW', 'REP'}, 0, 0);
        data.(L).AMP_MI_binned = modIndex(1,:);
        data.(L).HW_MI_binned  = modIndex(2,:);
        data.(L).TPW_MI_binned = modIndex(3,:);
        data.(L).REP_MI_binned = modIndex(4,:);
        
        % store smoothed data for each measure
        if sum(isnan(data.(L).AMP_microV_byBin)) == length(data.(L).AMP_microV_byBin)
            data.(L).AMP_microV_byBin_smoothed    = nan(1,cfg.phase.N_phase_bins);
        else
            data.(L).AMP_microV_byBin_smoothed    = DAG_circ_smooth(data.(L).AMP_microV_byBin);
        end
        if sum(isnan(data.(L).HW_ms_byBin)) == length(data.(L).HW_ms_byBin)
            data.(L).HW_ms_byBin_smoothed         = nan(1,cfg.phase.N_phase_bins);
        else
            data.(L).HW_ms_byBin_smoothed         = DAG_circ_smooth(data.(L).HW_ms_byBin);
        end
        if sum(isnan(data.(L).TPW_ms_byBin)) == length(data.(L).TPW_ms_byBin)
            data.(L).TPW_ms_byBin_smoothed        = nan(1,cfg.phase.N_phase_bins);
        else
            data.(L).TPW_ms_byBin_smoothed        = DAG_circ_smooth(data.(L).TPW_ms_byBin);
        end
        if sum(isnan(data.(L).REP_ms_byBin)) == length(data.(L).REP_ms_byBin)
            data.(L).REP_ms_byBin_smoothed        = nan(1,cfg.phase.N_phase_bins);
        else
            data.(L).REP_ms_byBin_smoothed        = DAG_circ_smooth(data.(L).REP_ms_byBin);
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
        [cc, pp]                       = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).AMP_microV_byBin_smoothed);
        data.(L).cc_PSTH_feature(1)    = cc(2,1); % for spike AMP
        data.(L).pp_PSTH_feature(1)    = pp(2,1);
        data.(L).pperm_PSTH_feature(1) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed', data.(L).AMP_microV_byBin_smoothed, cfg.phase.n_shuffles, cfg.phase.tail, cfg.phase.alpha_level, cfg.phase.stat, cfg.phase.reports, cfg.phase.seed_state);
%         [c,p] = corrcoef(data.(L).spike_phases_radians, data.(L).AMP_microV);
        
        [cc, pp]                       = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).HW_ms_byBin_smoothed);
        data.(L).cc_PSTH_feature(2)    = cc(2,1); % for HW
        data.(L).pp_PSTH_feature(2)    = pp(2,1);
%         data.(L).pperm_PSTH_feature(2) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed', data.(L).HW_ms_byBin_smoothed, cfg.phase.correlation.n_permutations, 0, 0.05, 'linear', 0);
        
        [cc, pp]                       = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).TPW_ms_byBin_smoothed);
        data.(L).cc_PSTH_feature(3)    = cc(2,1); % for spike TPW
        data.(L).pp_PSTH_feature(3)    = pp(2,1);
%         data.(L).pperm_PSTH_feature(3) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed', data.(L).TPW_ms_byBin_smoothed, cfg.phase.correlation.n_permutations, 0, 0.05, 'linear', 0);
        
        [cc, pp]                       = corrcoef(data.(L).spike_phases_histogram_smoothed, data.(L).REP_ms_byBin_smoothed);
        data.(L).cc_PSTH_feature(4)    = cc(2,1); % for spike REP
        data.(L).pp_PSTH_feature(4)    = pp(2,1);
%         data.(L).pperm_PSTH_feature(4) = mult_comp_perm_corr(data.(L).spike_phases_histogram_smoothed', data.(L).REP_ms_byBin_smoothed, cfg.phase.correlation.n_permutations, 0, 0.05, 'linear', 0);
        
    end
    save([cfg.cardioballistic_folder filesep data.unitId '_' data.target '__spikes_ECGphase.mat'], 'data', '-v7.3')
    clear data
end
end

function real_data = average_smooth_data(input, cfg)
Kernel = normpdf(-5*cfg.time.gaussian_kernel:cfg.time.PSTH_binwidth:5*cfg.time.gaussian_kernel,0,cfg.time.gaussian_kernel); % build kernel
tmp = [ [zeros(1,80); input(:,1:end-1)'] ...
    input' ...
    [input(:,2:end)'; zeros(1,80)] ];
conv_data = conv2(tmp,Kernel,'same');
real_data = conv_data(:,cfg.phase.N_phase_bins+1:end-cfg.phase.N_phase_bins);
% average_data = mean(conv_data);
% std_data = std(conv_data);
% ste_data = sterr(conv_data);
end

function output = cardioballistic_fit(feature_data, eventPhases, cfg)
scaled_feature = feature_data / nanmean(feature_data);

% compute starting point for scaling factor
a1 = (max(scaled_feature) - min(scaled_feature))/2;

% compute starting point for phase
[y,dropnan]  = rmmissing(scaled_feature);
ph           = eventPhases(~dropnan);
b1 = mod(circ_mean(ph, y'), 2*pi); % add modulo by 2pi as circ_mean and circ_median can return negative output even having input within 0-2pi

% starting point for intercept
c1 = nanmean(scaled_feature);

startPoint = [a1 b1 c1];

[fittedmdl,gof,~,lin_mdl] = cosine_fit(eventPhases, scaled_feature', cfg, startPoint);

coefs = coeffvalues(fittedmdl);
coefs(2) = mod(coefs(2),2*pi);

% 5 coefficient related to cosine fitting
% - the modulation index, the slope of the cosine function
% - p-value of the modulation index
% - phase of modulation
% - R-squared
% - intercept from linear model
output = [coefs(1) lin_mdl.Coefficients.pValue(2) coefs(2) gof.rsquare lin_mdl.Coefficients.Estimate(1)];

end

function [max_consec_bins, feature_modulation_index] = significant_bins(average_real, lowerPercentile_2_5, upperPercentile_97_5, average_reshuffled)
% figure out significant differences of spike feature dynamics
sig_above = average_real > upperPercentile_97_5;
sig_below = average_real < lowerPercentile_2_5;

consec_above = diff([0 find(diff(sig_above)) numel(sig_above)]); % number of repetitions of the same element
consec_below =  diff([0 find(diff(sig_below)) numel(sig_below)]);

if sig_above(1) % if starts with zero
    clust_above = consec_above(1:2:end); % then take lengths of non-zero clusters
else
    clust_above = consec_above(2:2:end);
end

if sig_below(1) % if starts with zero
    clust_below = consec_below(1:2:end);
else
    clust_below = consec_below(2:2:end);
end

max_consec_bins = max([clust_above clust_below]); % max number of consecutive significant bins

if isempty(max_consec_bins)
    max_consec_bins = 0;
end

feature_modulation_index = ...
    (max(average_real) - min(average_real)) / mean(average_reshuffled, 'omitnan');

end

function [fittedmdl,gof,yfit,lin_mdl] = cosine_fit(x, y, cfg, startPoint_cos)
% drop nans
[y,dropnans] = rmmissing(y);
x            = x(~dropnans);

% do the non-linear fit
[fittedmdl,gof] = fit(double(x),y,cfg.fit.cos_mod,'StartPoint', startPoint_cos, 'Lower', cfg.fit.cos_lower, 'Upper', cfg.fit.cos_upper);

coefs = coeffvalues(fittedmdl); % get model coefficients
coefs(2) = mod(coefs(2),2*pi);

yfit_all = cfg.fit.cos_mod(coefs(1), coefs(2), coefs(3), x);
yfit     = cfg.fit.cos_mod(coefs(1), coefs(2), coefs(3), cfg.phase.phase_bin_centers);

% employ a linear fit to get a p-value vs. fitting with a
% constant model
lin_mdl = fitlm(y(:), yfit_all);

end

function spikes_realigned_microV = shift2peak(wf_times_interp_ms, waveforms_upsampled)
peak_idx = discretize(wf_times_interp_ms, [0.32 0.48]) == 1;
[~, idx] = max(abs(waveforms_upsampled(:,peak_idx)), [], 2); % search for max only around the trough
idx = idx+find(peak_idx, 1, 'first')-1; % return to indices of the interpolated data
spikes_realigned_microV = bsxfun(@(x,y) circshift(x,y,1), waveforms_upsampled', 37 - idx'); % shift to peak (idx 37 for the interpolated data), in microvolts
spikes_realigned_microV = spikes_realigned_microV';
end

function [lowerPercentile_2_5, upperPercentile_97_5, average_reshuffled] = compute_reshuffles(data, bin, cfg)
%% compute reshuffles
rng(0)
[~, reshuffled_spike_order] = sort(rand(cfg.n_permutations, length(data)), 2); % get random order of elements
data_reshuffled      = data(reshuffled_spike_order);
data_reshuffled      = arrayfun(@(x) nanmean(data_reshuffled(:, bin == x),2), 1:cfg.N_phase_bins, 'UniformOutput', false); % mean by phase
data_reshuffled      = cat(2, data_reshuffled{:});
average_reshuffled = mean(data_reshuffled, 1, 'omitnan');
lowerPercentile_2_5    = prctile(data_reshuffled, 2.5, 1);
upperPercentile_97_5  = prctile(data_reshuffled, 97.5, 1);
end

% function Y_microV = upsample_spikes(waveforms, wf_times, wf_times_interp)
% %% Luba's interpolation
% y = waveforms;
% y2 = interp1(wf_times, y', wf_times_interp, 'cubic');
% peak_idx = discretize(wf_times_interp, [0.32 0.48]) == 1;
% [~, idx] = max(abs(y2(peak_idx,:))); % search for max only around the trough
% idx = idx+find(peak_idx, 1, 'first')-1; % return to indices of the interpolated data
% 
% Y_microV = bsxfun(@(x,y) 10^6 * circshift(x,y), y2, 37 - idx); % shift to peak (idx 37 for the interpolated data), in microvolts
% end

% function AMP = compute_AMP(waveforms, spike_phases, wf_times, wf_times_interp, phase_bins, nReshuffles)
% %% phase bins
% [AMP.counts,~,AMP.bin] = histcounts(spike_phases, phase_bins);
% AMP.bin = single(AMP.bin);
% 
% %% Luba's interpolation
% y = waveforms;
% y2 = interp1(wf_times, y', wf_times_interp, 'cubic');
% peak_idx = discretize(wf_times_interp, [0.32 0.48]) == 1;
% [~, idx] = max(abs(y2(peak_idx,:))); % search for max only around the trough
% idx = idx+find(peak_idx, 1, 'first')-1; % return to indices of the interpolated data
% 
% AMP.waveforms_microV = bsxfun(@(x,y) 10^6 * circshift(x,y), y2, 37 - idx); % shift to peak (idx 37 for the interpolated data), in microvolts
% 
% AMP.real_microV                   = AMP.waveforms_microV(37,:); % real peak amplitude
% AMP.abs_microV                    = abs(AMP.real_microV); % absolute peak amplitude
% AMP.sign_peak                     = sign(mean(AMP.real_microV)); % signum at the peak
% AMP.mean_by_phase                 = arrayfun(@(x) mean(AMP.abs_microV(AMP.bin == x)), unique(AMP.bin)); % mean by phase
% AMP.mean_by_phase_smoothed        = circ_smooth(AMP.mean_by_phase, 16); % pi/128 * 16 --> pi/8
% AMP.std_by_phase                  = arrayfun(@(x) std(AMP.abs_microV(AMP.bin == x)), unique(AMP.bin)); % standard deviation by phase
% 
% %% compute reshuffles
% [~, reshuffled_spike_order] = sort(rand(nReshuffles, size(waveforms,1)), 2); % get random order of elements
% AMP_reshuffled = AMP.abs_microV(reshuffled_spike_order);
% AMP_reshuffled      = arrayfun(@(x) mean(AMP_reshuffled(:, AMP.bin == x),2), unique(AMP.bin), 'UniformOutput', false); % mean by phase
% AMP_reshuffled      = cat(2, AMP_reshuffled{:});
% AMP.lowerPrctile_2_5    = prctile(AMP_reshuffled, 2.5, 1);
% AMP.uppperPrctile_97_5  = prctile(AMP_reshuffled, 97.5, 1);
% 
% % average waveforms by phase and then figure out spike parameters
% AMP.WF_by_phase = arrayfun(@(x) mean(AMP.waveforms_microV(:,AMP.bin == x),2), unique(AMP.bin), 'UniformOutput', false);
% AMP.WF_by_phase = cat(2, AMP.WF_by_phase{:});
% AMP.mean_by_bin = AMP.WF_by_phase(37,:)';
% end



























%         % reshuffle based on the computed median
%         %% lowIBI - set up settings
%         cfg.time.IBI       = 1;
%         cfg.time.IBI_low   = 1;
%         cfg.time.IBI_high  = 0;
%         cfg.time.IBI_thrsh = data.(L).IBI_median* ones(1, length(cfg.condition));
%         Rpeaks_lowIBI      = ecg_bna_compute_session_shuffled_Rpeaks(sessions_info,cfg.time);
%         for XXX=1:numel(Rpeaks_lowIBI)
%             Rpeaks_lowIBI(XXX).RPEAK_ts=Rpeaks_lowIBI(XXX).RPEAK_ts-Rpeaks_lowIBI(XXX).offset+Triggers(XXX).offset;
%             Rpeaks_lowIBI(XXX).shuffled_ts=Rpeaks_lowIBI(XXX).shuffled_ts-Rpeaks_lowIBI(XXX).offset+Triggers(XXX).offset;
%             Rpeaks_lowIBI(XXX).offset=Triggers(XXX).offset;
%         end
        
%         %% compute within trial lowIBI RR intervals for real and shuffled data
%         % compute RR-intervals
%         lowIBI_valid_RRinterval_ends      = single([Rpeaks_lowIBI(b).(['RPEAK_ts' cfg.condition(c).Rpeak_field])]);
%         lowIBI_valid_RRinterval_starts    = single(lowIBI_valid_RRinterval_ends - [Rpeaks_lowIBI(b).(['RPEAK_dur' cfg.condition(c).Rpeak_field])]);
%         % compute shuffled RR-intervals
%         lowIBI_shuffled_RRinterval_ends   = single([Rpeaks_lowIBI(b).(['shuffled_ts' cfg.condition(c).Rpeak_field])]);
%         lowIBI_shuffled_RRinterval_starts = single(lowIBI_shuffled_RRinterval_ends - [Rpeaks_lowIBI(b).(['shuffled_dur' cfg.condition(c).Rpeak_field])]);
%         
%         lowIBI_shuffled_RRinterval_ends   = mat2cell(lowIBI_shuffled_RRinterval_ends, ones(size(lowIBI_shuffled_RRinterval_ends,1),1), size(lowIBI_shuffled_RRinterval_ends,2));
%         lowIBI_shuffled_RRinterval_starts = mat2cell(lowIBI_shuffled_RRinterval_starts, ones(size(lowIBI_shuffled_RRinterval_starts,1),1), size(lowIBI_shuffled_RRinterval_starts,2));
%         
% %         % 0. figure out RR-intervals lying within trials
% %         lowIBI_trial_starts_one_stream    = cellfun(@(x,y,z) x+y+Rpeaks_lowIBI([Rpeaks_lowIBI.block] == z).offset, state2_times, TDT_ECG1_t0_from_rec_start, block_nums);
% %         lowIBI_trial_ends_one_stream      = cellfun(@(x,y,z) x+y+Rpeaks_lowIBI([Rpeaks_lowIBI.block] == z).offset, state90_times, TDT_ECG1_t0_from_rec_start, block_nums);
%         
%         % Create a logical array by vectorized comparison
%         lowIBI_RR_within_trial_idx = any(lowIBI_valid_RRinterval_starts' > trial_starts_one_stream & ...
%             lowIBI_valid_RRinterval_ends' < trial_ends_one_stream, 2);
%         
%         % Use logical indexing to filter valid RR intervals
%         lowIBI_valid_RRinterval_starts = lowIBI_valid_RRinterval_starts(lowIBI_RR_within_trial_idx);
%         lowIBI_valid_RRinterval_ends   = lowIBI_valid_RRinterval_ends(lowIBI_RR_within_trial_idx);
%         
%         % reshuffled - get rid of RRs beyond the current set of trials
%         for shuffNum = 1:length(lowIBI_shuffled_RRinterval_starts)
%             lowIBI_shuffledRR_within_trial_idx = any(lowIBI_shuffled_RRinterval_starts{shuffNum}' > trial_starts_one_stream & ...
%                 lowIBI_shuffled_RRinterval_ends{shuffNum}' < trial_ends_one_stream, 2);
%             
%             lowIBI_shuffled_RRinterval_ends{shuffNum}   = lowIBI_shuffled_RRinterval_ends{shuffNum}(lowIBI_shuffledRR_within_trial_idx); 
%             lowIBI_shuffled_RRinterval_starts{shuffNum} = lowIBI_shuffled_RRinterval_starts{shuffNum}(lowIBI_shuffledRR_within_trial_idx);
%         end
%         
%         lowIBI_RRs   = lowIBI_valid_RRinterval_ends - lowIBI_valid_RRinterval_starts;
%         
%         %% highIBI - set up settings
%         cfg.time.IBI       = 1;
%         cfg.time.IBI_low   = 0;
%         cfg.time.IBI_high  = 1;
%         cfg.time.IBI_thrsh = data.(L).IBI_median* ones(1, length(cfg.condition));
%         Rpeaks_highIBI      = ecg_bna_compute_session_shuffled_Rpeaks(sessions_info,cfg.time);
%         for XXX=1:numel(Rpeaks_highIBI)
%             Rpeaks_highIBI(XXX).RPEAK_ts=Rpeaks_highIBI(XXX).RPEAK_ts-Rpeaks_highIBI(XXX).offset+Triggers(XXX).offset;
%             Rpeaks_highIBI(XXX).shuffled_ts=Rpeaks_highIBI(XXX).shuffled_ts-Rpeaks_highIBI(XXX).offset+Triggers(XXX).offset;
%             Rpeaks_highIBI(XXX).offset=Triggers(XXX).offset;
%         end
%         
%         % !!! IMPORTANT !!! erase settings related to median split
%         % reshuffles and re-initiate them for the next unit
%         cfg.time = rmfield(cfg.time, {'IBI', 'IBI_low', 'IBI_high', 'IBI_thrsh'});
%         
%         %% compute within trial highIBI RR intervals for real and shuffled data
%         % compute RR-intervals
%         highIBI_valid_RRinterval_ends      = single([Rpeaks_highIBI(b).(['RPEAK_ts' cfg.condition(c).Rpeak_field])]);
%         highIBI_valid_RRinterval_starts    = single(highIBI_valid_RRinterval_ends - [Rpeaks_highIBI(b).(['RPEAK_dur' cfg.condition(c).Rpeak_field])]);
%         % compute shuffled RR-intervals
%         highIBI_shuffled_RRinterval_ends   = single([Rpeaks_highIBI(b).(['shuffled_ts' cfg.condition(c).Rpeak_field])]);
%         highIBI_shuffled_RRinterval_starts = single(highIBI_shuffled_RRinterval_ends - [Rpeaks_highIBI(b).(['shuffled_dur' cfg.condition(c).Rpeak_field])]);
%         
%         highIBI_shuffled_RRinterval_ends   = mat2cell(highIBI_shuffled_RRinterval_ends, ones(size(highIBI_shuffled_RRinterval_ends,1),1), size(highIBI_shuffled_RRinterval_ends,2));
%         highIBI_shuffled_RRinterval_starts = mat2cell(highIBI_shuffled_RRinterval_starts, ones(size(highIBI_shuffled_RRinterval_starts,1),1), size(highIBI_shuffled_RRinterval_starts,2));
%         
% %         % 0. figure out RR-intervals lying within trials
% %         trial_starts_one_stream    = cellfun(@(x,y,z) x+y+Rpeaks_highIBI([Rpeaks_highIBI.block] == z).offset, state2_times, TDT_ECG1_t0_from_rec_start, block_nums);
% %         trial_ends_one_stream      = cellfun(@(x,y,z) x+y+Rpeaks_highIBI([Rpeaks_highIBI.block] == z).offset, state90_times, TDT_ECG1_t0_from_rec_start, block_nums);
%         
%         % Create a logical array by vectorized comparison
%         highIBI_RR_within_trial_idx = any(highIBI_valid_RRinterval_starts' > trial_starts_one_stream & ...
%             highIBI_valid_RRinterval_ends' < trial_ends_one_stream, 2);
%         
%         % Use logical indexing to filter valid RR intervals
%         highIBI_valid_RRinterval_starts = highIBI_valid_RRinterval_starts(highIBI_RR_within_trial_idx);
%         highIBI_valid_RRinterval_ends   = highIBI_valid_RRinterval_ends(highIBI_RR_within_trial_idx);
%         
%         % reshuffled - get rid of RRs beyond the current set of trials
%         for shuffNum = 1:length(highIBI_shuffled_RRinterval_starts)
%             highIBI_shuffledRR_within_trial_idx = any(highIBI_shuffled_RRinterval_starts{shuffNum}' > trial_starts_one_stream & ...
%                 highIBI_shuffled_RRinterval_ends{shuffNum}' < trial_ends_one_stream, 2);
%             
%             highIBI_shuffled_RRinterval_ends{shuffNum}   = highIBI_shuffled_RRinterval_ends{shuffNum}(highIBI_shuffledRR_within_trial_idx); 
%             highIBI_shuffled_RRinterval_starts{shuffNum} = highIBI_shuffled_RRinterval_starts{shuffNum}(highIBI_shuffledRR_within_trial_idx);
%         end
%         
%         highIBI_RRs   = highIBI_valid_RRinterval_ends - highIBI_valid_RRinterval_starts;
%         
% %         lowIBIids  = realPSTHs.RDs{1} < Output.(L).IBI_median;
% %         highIBIids = realPSTHs.RDs{1} > Output.(L).IBI_median;
%         
%         data.(L).lowIBI_timeRRstart   = lowIBI_valid_RRinterval_starts;
%         data.(L).lowIBI_timeRRend     = lowIBI_valid_RRinterval_ends;
%         data.(L).lowIBI_starts        = lowIBI_RRs;
%         data.(L).lowIBI_meanHR_bpm    = mean(60 ./ data.(L).lowIBI_starts);
%         data.(L).lowIBI_medianHR_bpm  = median(60 ./ data.(L).lowIBI_starts);
%         data.(L).lowIBI_stdHR_bpm     = std(60 ./ data.(L).lowIBI_starts);
%         data.(L).lowIBI_SDNN_ms       = std(1000 * data.(L).lowIBI_starts);
%         
%         data.(L).highIBI_timeRRstart  = highIBI_valid_RRinterval_starts;
%         data.(L).highIBI_timeRRend    = highIBI_valid_RRinterval_ends;
%         data.(L).highIBI_starts       = highIBI_RRs;
%         data.(L).highIBI_meanHR_bpm   = mean(60 ./ data.(L).highIBI_starts);
%         data.(L).highIBI_medianHR_bpm = median(60 ./ data.(L).highIBI_starts);
%         data.(L).highIBI_stdHR_bpm    = std(60 ./ data.(L).highIBI_starts);
%         data.(L).highIBI_SDNN_ms      = std(1000 * data.(L).highIBI_starts);


function [during_trial_index,bins]=get_valid_PSTH_indexes(PSTH_time,trial_onsets,trial_ends,cfg)
%% define which parts of the continous PSTH are during a trial
trial_onset_samples=ceil((trial_onsets-PSTH_time(1))/cfg.spk.PSTH_binwidth);
trial_ends_samples=floor((trial_ends-PSTH_time(1))/cfg.spk.PSTH_binwidth);
trial_onset_samples=trial_onset_samples+1;
trial_ends_samples=trial_ends_samples+1;

bins_before=round(cfg.window(1)/cfg.spk.PSTH_binwidth);
bins_after=round(cfg.window(2)/cfg.spk.PSTH_binwidth);
bins=bins_before:bins_after;

during_trial_index=false(size(PSTH_time));
drop_samples = trial_onset_samples < 1 | trial_ends_samples < 1; % in very rare cases samples have negative values, drop those (??)
trial_onset_samples = trial_onset_samples(~drop_samples);
trial_ends_samples = trial_ends_samples(~drop_samples);
for t=1:numel(trial_onset_samples)
    during_trial_index(trial_onset_samples(t)-(bins_before-1):trial_ends_samples(t)-bins_after)=true;
end
during_trial_index(1:-bins_before)=false;
during_trial_index(end-bins_after:end)=false;
%
% %% remove samples that would land outside
% rpeaks2exclude = ...
%     RPEAK_samples>=numel(SD)-bins_after | RPEAK_samples<=-bins_before;
%
% RPEAK_ts(rpeaks2exclude)       = NaN;
% RPEAK_samples(rpeaks2exclude)  = NaN;
% RPEAK_dur(rpeaks2exclude)      = NaN;
end

function out=compute_PSTH(trig,RAST,SD,PSTH_time,during_trial_index,bins,mode,cfg)

switch mode
    case 'real'
        ts=trig.ts;
        iv=trig.intervals;
    case 'shuffled'
        ts=trig.shuffled_ts;
        iv=trig.shuffled_intervals;
end

ts_samples=round((ts-PSTH_time(1))/cfg.spk.PSTH_binwidth);

%% preallocate "out" structure
n=size(ts_samples,1);
%m=numel(bins);

            out.raster = false(size(bins));
% out.mean      = nan([n m]);
% out.std       = nan([n m]);
% out.SEM       = nan([n m]);
% out.n_events  = NaN([n 1]);
% out.RTs       = {};
% out.RDs       = {};

%% loop through rows of RPEAK_samples: 1 row for real, nReshuffles rows of reshuffled data
for p=1:n
    % %% remove samples that would land outside
    valid = ts_samples(p,:)<=numel(SD) & ts_samples(p,:)>0;
    
    %RT=ts(p,valid);      % take time-stamps
    TS=ts_samples(p,valid); % take samples
    IV=iv(p,valid);     % take durations
    
    within_trial_idx = during_trial_index(TS);
    
    %RT=RT(within_trial_idx);
    TS=TS(within_trial_idx);
    IV=IV(within_trial_idx);
    if isempty(TS)
            PSTH = NaN(size(bins));
    else
        idx_by_event = bsxfun(@plus,TS',bins); % bsxfun produces a RPEAK-(nBins-1) matrix with samples taken from SDF for each R-peak
        if n == 1
            out.raster = RAST(idx_by_event);
        end
        PSTH = SD(idx_by_event);
    end
    
    out.mean(p,:)=mean(PSTH, 1);
    out.std(p,:) = std(PSTH, [], 1);
    out.SEM(p,:)=sterr(PSTH,1);
    out.n_events(p)=numel(TS);
    %out.RTs{p}=RT;
    out.intervals{p}=IV;
end
end

