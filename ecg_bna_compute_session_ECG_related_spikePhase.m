function ecg_bna_compute_session_ECG_related_spikePhase(session_info,Rpeaks,ecg_bna_cfg)

basepath_to_save=[session_info.SPK_fldr filesep 'cardioballistic'];
if ~exist(basepath_to_save,'dir')
    mkdir(basepath_to_save);
end

%% settings for this analysis - subject to be moved to the ecg_bna_cfg
nReshuffles = ecg_bna_cfg.n_permutations;
Fs = 2.441406250000000e+04; % sampling frequency of BB signal, Hz
wf_times_ms = 1000 * (1/Fs:1/Fs:32/Fs); % in ms
wf_times_interp_ms = 1000 * (1/4/Fs:1/4/Fs:32/Fs); % in ms
peak_id = 10; % sample number of the trough in the spike waveform
phase_bins          = linspace(0, 2*pi, 64);
% phase_bin_centers   = phase_bins + phase_bins(2)/2;
WCfolder    = ['Y:\Data\Sortcodes\Magnus_phys\' session_info.Date '\WC\'];

condition_labels={'Rest','Task'};

load(session_info.Input_spikes, 'population'); % 'population' structure comes from here
load(session_info.Input_trials, 'trials'); % 'trials' structure comes from here

Rblocks=[Rpeaks.block];

% one needs a population and a by_block file to make it work
for unitNum = 1:length(population)
    
    disp(['Processing unit ' num2str(unitNum) ' out of ' num2str(length(population))])
    
    pop=population(unitNum);
    
    T=ph_get_unit_trials(pop,trials);
    
    T_acc=[T.accepted] & [T.completed];
    T=T(T_acc);
    pop.trial=pop.trial(T_acc);
    
    %% Make sure we only take overlapping blocks
    blocks_unit=unique([pop.block]);
    blocks=intersect(blocks_unit,Rblocks);
    b=ismember(Rblocks,blocks);
    
    data.unitId                               = pop.unit_ID;
    data.channel                              = pop.channel;
    data.target                               = pop.target;
    data.SNR                                  = pop.avg_SNR;
    data.stability                            = pop.avg_stability;
    data.single_rating                        = pop.avg_single_rating;
    data.FR                                   = mean([pop.FR_average]);
    
    for tasktype = 1:2
        tType = condition_labels{tasktype};
        
        %% take trials from the necessary blocks and current tasktype
        tr=ismember([T.block],blocks) & [T.type]==tasktype;
        popcell=pop.trial(tr);
        trcell=T(tr);
        
        % find the corresponding WC file and load it
        chNum = ['00' num2str(population(unitNum).channel)];
        chNum = chNum(end-2:end);
        WCfile = dir([WCfolder filesep 'dataspikes_ch' chNum '*.mat']);
        if length(WCfile) == 1
            WC = load([WCfile.folder filesep WCfile.name]);
        end
        if length(WC.thr) == 4
            data.(tType).thresholds_microV = 10^6 * WC.thr;
            data.(tType).thresholds_microV(3:4) = -1*data.(tType).thresholds_microV(3:4);
        end
        clear WC
        % 0. Prepare data variables
        states_onset               = {trcell.states_onset};
        states                     = {trcell.states};
        TDT_ECG1_t0_from_rec_start = {trcell.TDT_ECG1_t0_from_rec_start};
        block_nums                 = {trcell.block};
        % compute RR-intervals
        valid_RRinterval_ends =             single([Rpeaks(b).RPEAK_ts]);
        valid_RRinterval_starts =           single(valid_RRinterval_ends - [Rpeaks(b).RPEAK_dur]);
        % 1. take arrival times and the corresponding waveforms
        AT = {popcell.arrival_times};
        WF = {popcell.waveforms};
        % 2. choose only those that happen after MP-state 1
        state1_times = cellfun(@(x,y) x(y == 1), states_onset, states, 'Uniformoutput', false);
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
        data.(tType).spike_phases_radians                 = eventPhases;
        
        [~,~,bin] = histcounts(data.(tType).spike_phases_radians, phase_bins);
        
        data.(tType).waveforms_microvolts                 = 10^6 * WF_one_stream(eventsTaken,:);
        waveforms_upsampled                       = interpft(data.(tType).waveforms_microvolts, 32*4, 2);
        data.(tType).waveforms_upsampled_microvolts       = shift2peak(wf_times_interp_ms, waveforms_upsampled);
        data.(tType).waveforms_byBin_microvolts           = arrayfun(@(x) mean(data.(tType).waveforms_upsampled_microvolts(bin == x,:),1), unique(bin), 'UniformOutput', false);
        data.(tType).waveforms_byBin_microvolts           = cat(1,data.(tType).waveforms_byBin_microvolts{:});
        % 7. Calculate spike features with Mosher's procedure
        tic
        sMetric = ...
            struct('extremAmp', cell(length(data.(tType).spike_phases_radians),1), 'widthHW', cell(length(data.(tType).spike_phases_radians),1), 'widthTP', cell(length(data.(tType).spike_phases_radians),1), 'repolTime', cell(length(data.(tType).spike_phases_radians),1));
        for wfNum = 1:length(data.(tType).spike_phases_radians)
            try
                sMetric(wfNum)=spikeWaveMetrics(data.(tType).waveforms_upsampled_microvolts(wfNum,:), 37, Fs*4, 0); % 37 - index of th peak for updsampled data
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
        data.(tType).AMP_microV               = [sMetric.extremAmp];
        data.(tType).HW_ms                    = 10^3 * [sMetric.widthHW];
        data.(tType).TPW_ms                   = 10^3 * [sMetric.widthTP];
        data.(tType).REP_ms                   = 10^3 * [sMetric.repolTime];
        
        data.(tType).AMP_microV_byBin         = arrayfun(@(x) nanmean(data.(tType).AMP_microV(bin == x)), unique(bin)); % mean by phase
        data.(tType).HW_ms_byBin              = arrayfun(@(x) nanmean(data.(tType).HW_ms(bin == x)), unique(bin));
        data.(tType).TPW_ms_byBin             = arrayfun(@(x) nanmean(data.(tType).TPW_ms(bin == x)), unique(bin));
        data.(tType).REP_ms_byBin             = arrayfun(@(x) nanmean(data.(tType).REP_ms(bin == x)), unique(bin));
        
        [data.(tType).AMP_lowerPrctile_2_5, data.(tType).AMP_upperPrctile_97_5] = ...
            compute_reshuffles(data.(tType).AMP_microV, bin, nReshuffles);
        [data.(tType).HW_lowerPrctile_2_5, data.(tType).HW_upperPrctile_97_5] = ...
            compute_reshuffles(data.(tType).HW_ms, bin, nReshuffles);
        [data.(tType).TPW_lowerPrctile_2_5, data.(tType).TPW_upperPrctile_97_5] = ...
            compute_reshuffles(data.(tType).TPW_ms, bin, nReshuffles);
        [data.(tType).REP_lowerPrctile_2_5, data.(tType).REP_upperPrctile_97_5] = ...
            compute_reshuffles(data.(tType).REP_ms, bin, nReshuffles);
        
        [modIndex,removeNoise,allCorr,allLinMod] = ...
            fitCardiacModulation(phase_bins(1:end-1), ...
            [data.(tType).AMP_microV_byBin data.(tType).HW_ms_byBin data.(tType).TPW_ms_byBin data.(tType).REP_ms_byBin]', ...
            {'AMP', 'HW', 'TPW', 'REP'}, 0, [221 222 223 224]);
        
        % 4 coefficient related to cosine fitting
        % - the modulation index, the slope of the cosine function
        % - p-value of the modulation index
        % - phase of modulation
        % - p-value of the modulation index ?? (mdl.Rsquared.ordinary)
        data.(tType).AMP_MI                   = modIndex(1, :);
        data.(tType).HW_MI                    = modIndex(2, :);
        data.(tType).TPW_MI                   = modIndex(3, :);
        data.(tType).REP_MI                   = modIndex(4, :);
        
        % store smoothed data for each measure
        data.(tType).AMP_microV_byBin_smoothed    = removeNoise(1,:);
        data.(tType).HW_ms_byBin_smoothed         = removeNoise(2,:);
        data.(tType).TPW_ms_byBin_smoothed        = removeNoise(3,:);
        data.(tType).REP_ms_byBin_smoothed        = removeNoise(4,:);
        
        data.(tType).allCorr                  = allCorr;
        data.(tType).allLinMod                = allLinMod;
        
        % II. Compute unit firing rate per RR-interval
        [data.(tType).FRbyRR_Hz, ...
            data.(tType).cycleDurations_s] = ...
            computeFRperCycle(valid_RRinterval_starts, valid_RRinterval_ends, AT_one_stream);
        [temp_r, temp_p] = corrcoef(data.(tType).FRbyRR_Hz, data.(tType).cycleDurations_s);
        data.(tType).pearson_r = temp_r(2,1);
        data.(tType).pearson_p = temp_p(2,1);
        % III. Put all results together
        
    end
	save([basepath_to_save filesep data.unitId '_spikes_ECGphase.mat'], 'data', '-v7.3')
	clear data
end
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

function [lowerPercentile_2_5, upperPercentile_97_5] = compute_reshuffles(data, bin, nReshuffles)
%% compute reshuffles
rng(0)
[~, reshuffled_spike_order] = sort(rand(nReshuffles, length(data)), 2); % get random order of elements
data_reshuffled      = data(reshuffled_spike_order);
data_reshuffled      = arrayfun(@(x) nanmean(data_reshuffled(:, bin == x),2), unique(bin), 'UniformOutput', false); % mean by phase
data_reshuffled      = cat(2, data_reshuffled{:});
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
