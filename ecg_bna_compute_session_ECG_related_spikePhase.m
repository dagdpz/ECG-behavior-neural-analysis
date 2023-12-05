function ecg_bna_compute_session_ECG_related_spikePhase(trials,population,Rpeaks,cfg)

basepath_to_save=[cfg.SPK_root_results_fldr filesep 'cardioballistic'];
if ~exist(basepath_to_save,'dir')
    mkdir(basepath_to_save);
end

%% load list of selected units - rather go through all units really
% unit_list = load([session_info.SPK_fldr '\unitInfo_after_exclusion_stableTaskAndRest.mat']);
% unitList = unique(unit_list.unit_ids_after_exclusion);

%% settings for this analysis - subject to be moved to the cfg
Fs = 2.441406250000000e+04; % sampling frequency of BB signal, Hz
wf_times_ms = 1000 * (1/Fs:1/Fs:32/Fs); % in ms
wf_times_interp_ms = 1000 * (1/4/Fs:1/4/Fs:32/Fs); % in ms
peak_id = 10; % sample number of the trough in the spike waveform
phase_bins          = linspace(0, 2*pi, cfg.spk.N_phase_bins+1);
phase_bin_centers   = 2*pi/cfg.spk.N_phase_bins:2*pi/cfg.spk.N_phase_bins:2*pi;
lag_list = [-11 -7 -3 0 3 7 11];


Rblocks=[Rpeaks.block];

%% really would rather not exclude here as discussed
% figure out which units take from this session
% selected_this_session = ismember({population.unit_ID}, unitList);
% population = population(selected_this_session);

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
    data.unitId                                            = pop.unit_ID;
    data.target                                            = pop.target;
    data.channel                                           = pop.channel;
    data.unit                                              = pop.block_unit{2,1};
    data.quantSNR                                          = pop.avg_SNR;
    data.Single_rating                                     = pop.avg_single_rating;
    data.stability_rating                                  = pop.avg_stability;
    data.thresholds_microV                                 = single([0; 0; 0; 0]);
    data.FR                                                = single(mean(pop.FR_average));
    data.cc_lag_list                                       = lag_list;
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        data.(L).spike_phases_radians           = single(NaN);
        data.(L).waveforms_microvolts           = single(nan(1,32));
        data.(L).waveforms_upsampled_microvolts = single(nan(1,128));
        data.(L).waveforms_byBin_microvolts     = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).AMP_microV                     = single(NaN);
        data.(L).HW_ms                          = single(NaN);
        data.(L).TPW_ms                         = single(NaN);
        data.(L).REP_ms                         = single(NaN);
        data.(L).AMP_microV_byBin               = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).HW_ms_byBin                    = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).TPW_ms_byBin                   = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).REP_ms_byBin                   = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).AMP_reshuffled_avg             = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).AMP_lowerPrctile_2_5           = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).AMP_upperPrctile_97_5          = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).HW_reshuffled_avg              = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).HW_lowerPrctile_2_5            = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).HW_upperPrctile_97_5           = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).TPW_reshuffled_avg             = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).TPW_lowerPrctile_2_5           = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).TPW_upperPrctile_97_5          = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).REP_reshuffled_avg             = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).REP_lowerPrctile_2_5           = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).REP_upperPrctile_97_5          = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).AMP_MI                         = single(nan(1,5));
        data.(L).HW_MI                          = single(nan(1,5));
        data.(L).TPW_MI                         = single(nan(1,5));
        data.(L).REP_MI                         = single(nan(1,5));
        data.(L).AMP_microV_byBin_smoothed      = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).HW_ms_byBin_smoothed           = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).TPW_ms_byBin_smoothed          = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).REP_ms_byBin_smoothed          = single(nan(1,cfg.spk.N_phase_bins));
        data.(L).allCorr                        = single(nan(1,6));
        data.(L).allLinMod                      = single(nan(6,2));
        data.(L).AMP_max_consec_bins            = single(nan(1,1));
        data.(L).AMP_modulation_index           = single(nan(1,1));
        data.(L).HW_max_consec_bins             = single(nan(1,1));
        data.(L).HW_modulation_index            = single(nan(1,1));
        data.(L).TPW_max_consec_bins            = single(nan(1,1));
        data.(L).TPW_modulation_index           = single(nan(1,1));
        data.(L).REP_max_consec_bins            = single(nan(1,1));
        data.(L).REP_modulation_index           = single(nan(1,1));
        data.(L).FRbyRR_Hz                      = single(nan(1,1));
        data.(L).cycleDurations_s               = single(nan(1,1));
        data.(L).pearson_r                      = single(nan(length(lag_list), 1));
        data.(L).pearson_p                      = single(nan(length(lag_list), 1));
    end
    
    % find the corresponding WC file and load it
    chNum = data.channel;
    blkNum = unique([pop.block]);
    WCfile = ph_figure_out_waveclus_file_by_channel_and_blocks(chNum, blkNum, cfg.Input_WC);
    WC = load(WCfile);
    if length(WC.thr) == 4
        data.thresholds_microV = 10^6 * WC.thr;
        data.thresholds_microV(3:4) = -1*data.thresholds_microV(3:4);
    end
    clear WC
    
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        
        %% get condition AND valid block trials only
        CT = ecg_bna_get_condition_trials(T, cfg.condition(c));
        tr=ismember([T.block],blocks) & CT;
        if sum(tr)>1 % do calculations only if number of trials > 1
            popcell=pop.trial(tr);
            trcell=T(tr);
            
            % 0. Prepare data variables
            states_onset               = {trcell.states_onset};
            states                     = {trcell.states};
            TDT_ECG1_t0_from_rec_start = {trcell.TDT_ECG1_t0_from_rec_start};
            block_nums                 = {trcell.block};
            state1_times               = cellfun(@(x,y) x(y == 1), states_onset, states, 'Uniformoutput', false); % trial starts
            state98_times              = cellfun(@(x,y) x(y == 98), states_onset, states, 'Uniformoutput', false); % trial ends
            % compute RR-intervals
            valid_RRinterval_ends      = single([Rpeaks(b).RPEAK_ts]);
            valid_RRinterval_starts    = single(valid_RRinterval_ends - [Rpeaks(b).RPEAK_dur]);
            % 0. figure out RR-intervals lying within trials
            trial_starts_one_stream    = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, state1_times, TDT_ECG1_t0_from_rec_start, block_nums);
            trial_ends_one_stream      = cellfun(@(x,y,z) x+y+Rpeaks([Rpeaks.block] == z).offset, state98_times, TDT_ECG1_t0_from_rec_start, block_nums);
            
            RR_within_trial_idx = false(length(valid_RRinterval_starts),1);
            for RRnum = 1:length(RR_within_trial_idx)
                RR_within_trial_idx(RRnum) = ...
                    any(valid_RRinterval_starts(RRnum)>trial_starts_one_stream & ...
                    valid_RRinterval_ends(RRnum)<trial_ends_one_stream);
            end
            valid_RRinterval_ends      = valid_RRinterval_ends(RR_within_trial_idx); % get rid of RRs beyond the current set of trials
            valid_RRinterval_starts    = valid_RRinterval_starts(RR_within_trial_idx);
            % 1. take arrival times and the corresponding waveforms
            AT = {popcell.arrival_times};
            WF = {popcell.waveforms};
            % 2. choose only those that happen after MP-state 1
            idx_after_state1 = cellfun(@(x,y) x>y, AT, state1_times, 'Uniformoutput', false);
            % 3. add TDT_ECG1_t0_from_rec_start and Rpeak block offset to spike times
            AT_one_stream_cell = cellfun(@(x,y,z,a) x(y)+z+Rpeaks([Rpeaks.block] == a).offset, AT, idx_after_state1, TDT_ECG1_t0_from_rec_start, block_nums, 'Uniformoutput', false);
            WF_one_stream_cell = cellfun(@(x,y) x(y,:), WF, idx_after_state1, 'Uniformoutput', false);
            % 4. merge all spike times and waveforms together
            AT_one_stream = cat(1, AT_one_stream_cell{:});
            WF_one_stream = cat(1, WF_one_stream_cell{:});
            % 5. calculate heart cycle phase where individual spikes ended up
            [eventPhases, eventsTaken] = DAG_eventPhase(valid_RRinterval_starts, valid_RRinterval_ends, AT_one_stream);
            % 6. Put results for real data together
            data.(L).spike_phases_radians                 = eventPhases;
            
            [~,~,bin] = histcounts(data.(L).spike_phases_radians, phase_bins);
            
            data.(L).waveforms_microvolts                 = 10^6 * WF_one_stream(eventsTaken,:);
            waveforms_upsampled                       = interpft(data.(L).waveforms_microvolts, 32*4, 2);
            data.(L).waveforms_upsampled_microvolts       = shift2peak(wf_times_interp_ms, waveforms_upsampled);
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
                try %% spikeWaveMetrics missing! --> try catch is also VERY inefficient
                    sMetric(wfNum)=spikeWaveMetrics(double(data.(L).waveforms_upsampled_microvolts(wfNum,:)), 37, Fs*4, 0); % 37 - index of th peak for updsampled data
                    sMetric(wfNum).extremAmp   = sMetric(wfNum).extremAmp(1);
                    sMetric(wfNum).widthHW     = sMetric(wfNum).widthHW(1);
                    sMetric(wfNum).widthTP     = sMetric(wfNum).widthTP(1);
                    sMetric(wfNum).repolTime   = sMetric(wfNum).repolTime(1);
                catch ME
                    sMetric(wfNum).extremAmp  = nan;
                    sMetric(wfNum).widthHW    = nan;
                    sMetric(wfNum).widthTP    = nan;
                    sMetric(wfNum).repolTime  = nan;
                end
            end
            toc
            data.(L).AMP_microV               = [sMetric.extremAmp];
            data.(L).HW_ms                    = 10^3 * [sMetric.widthHW];
            data.(L).TPW_ms                   = 10^3 * [sMetric.widthTP];
            data.(L).REP_ms                   = 10^3 * [sMetric.repolTime];
            
            data.(L).AMP_microV_byBin         = arrayfun(@(x) nanmean(data.(L).AMP_microV(bin == x)), 1:cfg.spk.N_phase_bins); % mean by phase
            data.(L).HW_ms_byBin              = arrayfun(@(x) nanmean(data.(L).HW_ms(bin == x)), 1:cfg.spk.N_phase_bins);
            data.(L).TPW_ms_byBin             = arrayfun(@(x) nanmean(data.(L).TPW_ms(bin == x)), 1:cfg.spk.N_phase_bins);
            data.(L).REP_ms_byBin             = arrayfun(@(x) nanmean(data.(L).REP_ms(bin == x)), 1:cfg.spk.N_phase_bins);
            
            [data.(L).AMP_lowerPrctile_2_5, data.(L).AMP_upperPrctile_97_5, data.(L).AMP_reshuffled_avg] = ...
                compute_reshuffles(data.(L).AMP_microV, bin, cfg);
            [data.(L).HW_lowerPrctile_2_5, data.(L).HW_upperPrctile_97_5, data.(L).HW_reshuffled_avg] = ...
                compute_reshuffles(data.(L).HW_ms, bin, cfg);
            [data.(L).TPW_lowerPrctile_2_5, data.(L).TPW_upperPrctile_97_5, data.(L).TPW_reshuffled_avg] = ...
                compute_reshuffles(data.(L).TPW_ms, bin, cfg);
            [data.(L).REP_lowerPrctile_2_5, data.(L).REP_upperPrctile_97_5, data.(L).REP_reshuffled_avg] = ...
                compute_reshuffles(data.(L).REP_ms, bin, cfg);
            
            featureMatrix = ...
                [data.(L).AMP_microV_byBin; data.(L).HW_ms_byBin; ...
                data.(L).TPW_ms_byBin; data.(L).REP_ms_byBin];
            if size(featureMatrix,1) > 4
                featureMatrix = featureMatrix';
                if size(featureMatrix,1) ~= 4
                    error('Dimensions of feature matrix aren''t suitable for the analysis')
                end
            end
            
            [modIndex,removeNoise,allCorr,allLinMod] = ...
                fitCardiacModulation(phase_bins(1:end-1), ...
                featureMatrix, {'AMP', 'HW', 'TPW', 'REP'}, 0, [221 222 223 224]);
            % why does this one not exist yet ??
            modIndex(:,5)=nanmean(featureMatrix,2);
            % 4 coefficient related to cosine fitting
            % - the modulation index, the slope of the cosine function
            % - p-value of the modulation index
            % - phase of modulation
            % - p-value of the modulation index ?? (mdl.Rsquared.ordinary)
            data.(L).AMP_MI                   = modIndex(1, :);
            data.(L).HW_MI                    = modIndex(2, :);
            data.(L).TPW_MI                   = modIndex(3, :);
            data.(L).REP_MI                   = modIndex(4, :);
            
            % store smoothed data for each measure
            data.(L).AMP_microV_byBin_smoothed    = removeNoise(1,:);
            data.(L).HW_ms_byBin_smoothed         = removeNoise(2,:);
            data.(L).TPW_ms_byBin_smoothed        = removeNoise(3,:);
            data.(L).REP_ms_byBin_smoothed        = removeNoise(4,:);
            
            data.(L).allCorr                      = allCorr;
            data.(L).allLinMod                    = allLinMod;
            
            [data.(L).AMP_max_consec_bins, data.(L).AMP_modulation_index] = ...
                significant_bins(data.(L).AMP_microV_byBin_smoothed, data.(L).AMP_lowerPrctile_2_5, data.(L).AMP_upperPrctile_97_5, data.(L).AMP_reshuffled_avg);
            [data.(L).HW_max_consec_bins, data.(L).HW_modulation_index] = ...
                significant_bins(data.(L).HW_ms_byBin_smoothed, data.(L).HW_lowerPrctile_2_5, data.(L).HW_upperPrctile_97_5, data.(L).HW_reshuffled_avg);
            [data.(L).TPW_max_consec_bins, data.(L).TPW_modulation_index] = ...
                significant_bins(data.(L).TPW_ms_byBin_smoothed, data.(L).TPW_lowerPrctile_2_5, data.(L).TPW_upperPrctile_97_5, data.(L).TPW_reshuffled_avg);
            [data.(L).REP_max_consec_bins, data.(L).REP_modulation_index] = ...
                significant_bins(data.(L).REP_ms_byBin_smoothed, data.(L).REP_lowerPrctile_2_5, data.(L).REP_upperPrctile_97_5, data.(L).REP_reshuffled_avg);
            
            % II. Compute unit firing rate per RR-interval
            [data.(L).FRbyRR_Hz, ...
                data.(L).cycleDurations_s] = ...
                computeFRperCycle(valid_RRinterval_starts, valid_RRinterval_ends, AT_one_stream);
            % compute correlation with different lag
            for lagNum = 1:length(lag_list)
                [temp_r, temp_p] = corrcoef(data.(L).FRbyRR_Hz, circshift(data.(L).cycleDurations_s, lag_list(lagNum)));
                data.(L).pearson_r(lagNum) = temp_r(2,1);
                data.(L).pearson_p(lagNum) = temp_p(2,1);
            end
        end
    end
    save([basepath_to_save filesep data.unitId '_' data.target '__spikes_ECGphase.mat'], 'data', '-v7.3')
    clear data
end
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

function spikes_realigned_microV = shift2peak(wf_times_interp_ms, waveforms_upsampled)
peak_idx = discretize(wf_times_interp_ms, [0.32 0.48]) == 1;
[~, idx] = max(abs(waveforms_upsampled(:,peak_idx)), [], 2); % search for max only around the trough
idx = idx+find(peak_idx, 1, 'first')-1; % return to indices of the interpolated data
spikes_realigned_microV = bsxfun(@(x,y) circshift(x,y,1), waveforms_upsampled', 37 - idx'); % shift to peak (idx 37 for the interpolated data), in microvolts
spikes_realigned_microV = spikes_realigned_microV';
end

function [FRbyRR_Hz, cycleDurations] = computeFRperCycle(interval_starts, interval_ends, eventTimes)
% make inputs vertical vectors
if size(interval_starts, 2) > size(interval_starts, 1)
    interval_starts = interval_starts';
end
if size(interval_ends, 2) > size(interval_ends, 1)
    interval_ends = interval_ends';
end
if size(eventTimes, 2) > size(eventTimes, 1)
    eventTimes = eventTimes';
end

% Calculate ECG cycle durations
cycleDurations = interval_ends - interval_starts;

% cycleNums = arrayfun(@(x) sum(x > interval_starts & x < interval_ends), eventTimes, 'UniformOutput', false);
spikeCounts = arrayfun(@(x,y) sum(eventTimes > x & eventTimes < y), interval_starts, interval_ends);
FRbyRR_Hz = spikeCounts ./ cycleDurations;

% figure,
% scatter(cycleDurations, FRbyRR_Hz)
end

function [lowerPercentile_2_5, upperPercentile_97_5, average_reshuffled] = compute_reshuffles(data, bin, cfg)
%% compute reshuffles
rng(0)
[~, reshuffled_spike_order] = sort(rand(cfg.spk.n_permutations, length(data)), 2); % get random order of elements
data_reshuffled      = data(reshuffled_spike_order);
data_reshuffled      = arrayfun(@(x) nanmean(data_reshuffled(:, bin == x),2), 1:cfg.spk.N_phase_bins, 'UniformOutput', false); % mean by phase
data_reshuffled      = cat(2, data_reshuffled{:});
average_reshuffled = mean(data_reshuffled, 1, 'omitnan');
lowerPercentile_2_5    = prctile(data_reshuffled, 2.5, 1);
upperPercentile_97_5  = prctile(data_reshuffled, 97.5, 1);
end

function Y_microV = upsample_spikes(waveforms, wf_times, wf_times_interp)
%% Luba's interpolation
y = waveforms;
y2 = interp1(wf_times, y', wf_times_interp, 'cubic');
peak_idx = discretize(wf_times_interp, [0.32 0.48]) == 1;
[~, idx] = max(abs(y2(peak_idx,:))); % search for max only around the trough
idx = idx+find(peak_idx, 1, 'first')-1; % return to indices of the interpolated data

Y_microV = bsxfun(@(x,y) 10^6 * circshift(x,y), y2, 37 - idx); % shift to peak (idx 37 for the interpolated data), in microvolts
end

function AMP = compute_AMP(waveforms, spike_phases, wf_times, wf_times_interp, phase_bins, nReshuffles)
%% phase bins
[AMP.counts,~,AMP.bin] = histcounts(spike_phases, phase_bins);
AMP.bin = single(AMP.bin);

%% Luba's interpolation
y = waveforms;
y2 = interp1(wf_times, y', wf_times_interp, 'cubic');
peak_idx = discretize(wf_times_interp, [0.32 0.48]) == 1;
[~, idx] = max(abs(y2(peak_idx,:))); % search for max only around the trough
idx = idx+find(peak_idx, 1, 'first')-1; % return to indices of the interpolated data

AMP.waveforms_microV = bsxfun(@(x,y) 10^6 * circshift(x,y), y2, 37 - idx); % shift to peak (idx 37 for the interpolated data), in microvolts

AMP.real_microV                   = AMP.waveforms_microV(37,:); % real peak amplitude
AMP.abs_microV                    = abs(AMP.real_microV); % absolute peak amplitude
AMP.sign_peak                     = sign(mean(AMP.real_microV)); % signum at the peak
AMP.mean_by_phase                 = arrayfun(@(x) mean(AMP.abs_microV(AMP.bin == x)), unique(AMP.bin)); % mean by phase
AMP.mean_by_phase_smoothed        = circ_smooth(AMP.mean_by_phase, 16); % pi/128 * 16 --> pi/8
AMP.std_by_phase                  = arrayfun(@(x) std(AMP.abs_microV(AMP.bin == x)), unique(AMP.bin)); % standard deviation by phase

%% compute reshuffles
[~, reshuffled_spike_order] = sort(rand(nReshuffles, size(waveforms,1)), 2); % get random order of elements
AMP_reshuffled = AMP.abs_microV(reshuffled_spike_order);
AMP_reshuffled      = arrayfun(@(x) mean(AMP_reshuffled(:, AMP.bin == x),2), unique(AMP.bin), 'UniformOutput', false); % mean by phase
AMP_reshuffled      = cat(2, AMP_reshuffled{:});
AMP.lowerPrctile_2_5    = prctile(AMP_reshuffled, 2.5, 1);
AMP.uppperPrctile_97_5  = prctile(AMP_reshuffled, 97.5, 1);

% average waveforms by phase and then figure out spike parameters
AMP.WF_by_phase = arrayfun(@(x) mean(AMP.waveforms_microV(:,AMP.bin == x),2), unique(AMP.bin), 'UniformOutput', false);
AMP.WF_by_phase = cat(2, AMP.WF_by_phase{:});
AMP.mean_by_bin = AMP.WF_by_phase(37,:)';
end

function out = circ_smooth(input, smooth_win)

A = repmat(input, 3, 1);
A_smoothed = smooth(A, smooth_win);

out = A_smoothed(length(input)+1:end-length(input));

end
