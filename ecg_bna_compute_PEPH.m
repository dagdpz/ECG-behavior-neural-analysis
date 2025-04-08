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

Nbins      = cfg.phase.N_phase_bins;
phase_bin_edges         = linspace(0, 2*pi, Nbins+1);
phase_bins  = pi/Nbins : 2*pi/Nbins : 2*pi-pi/Nbins;

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
    
    %% do we need those??
    data.channel           = pop.channel;
    data.unit              = pop.block_unit{2,1};
    data.criteria          = pop.criteria;
    data.FR                = single(mean(pop.FR_average));
    
    
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
        
        if sum(tr)<=1
            continue
        end
        
        popcell=pop.trial(tr);
        trcell=T(tr);
        [O]=ecg_bna_synchronize_with_trigger(Blockoffsets,trcell,popcell);
        
        % 0. figure out RR-intervals lying within trials
        trial_starts    = O.trial_starts;
        trial_ends      = O.trial_ends;
        
        % E=events{e};
        E='R';
        
        % compute RR-intervals
        interval_ends      = Triggers.(E).ts;
        interval_starts    = Triggers.(E).ts  - Triggers.(E).intervals;
        trig_in_trial=any(bsxfun(@ge,interval_starts',trial_starts) & bsxfun(@le,interval_ends',trial_ends),2)';
        
        starts = interval_starts(trig_in_trial)';
        ends = interval_ends(trig_in_trial)';
        durations=ends-starts;
        spiketimes = O.AT;
        amps = 10^6 *O.amps;
        
        
        % compute parameters of heart activity - disregard trial limitations here (?)
        data.(L).meanHR_bpm   = mean(60./durations);
        data.(L).medianHR_bpm = median(60./durations);
        data.(L).stdHR_bpm    = std(60./durations);
        data.(L).SDNN_ms      = std(1000 *durations);
        data.(L).NrTrigs      = sum(trig_in_trial);
        
        % Normalize event times to the respective ECG cycle durations
        [events2include,cycleNums_withSpikes] = find(bsxfun(@ge,spiketimes,starts') & bsxfun(@le,spiketimes,ends'));
        ts = spiketimes(events2include); % include only events that landed within any interval
        amps = amps(events2include); % include only events that landed within any interval
        spikephases = 2*pi*(ts - starts(cycleNums_withSpikes)) ./ durations(cycleNums_withSpikes);
        %eventsTaken = single(find(events2include))';
        real_histogram = hist3([spikephases, cycleNums_withSpikes], 'ctrs', {phase_bins 1:length(starts)});
        real_SDF = real_histogram.*repmat(numel(phase_bins)./durations',numel(phase_bins),1);
        
        
        replace_nans=0;
        [spikesperbin, ~, bin] = histcounts(spikephases, phase_bin_edges);
        [~,sort_idx]=sort(-amps);
        minspikes=min(spikesperbin);
        amps_s=amps(sort_idx);
        bin_s=bin(sort_idx);
        if minspikes>=1
            data.(L).AMP_byBin             = arrayfun(@(x) nanmean(amps_s(find(bin_s == x,minspikes))), 1:Nbins); % mean by phase
            data.(L).AMP_byBin_smoothed    = DAG_circ_smooth(data.(L).AMP_byBin); % don't think this has to be smoothed (?)
        else
            data.(L).AMP_byBin             = arrayfun(@(x) nanmean(amps(bin == x)), 1:Nbins); % mean by phase
            %data.(L).AMP_byBin             = nan(1,Nbins); % set all NaNs since its not meaningful
            data.(L).AMP_byBin_smoothed    = nan(1,Nbins); % don't think this has to be smoothed (?)
        end
        
        data.(L).AMP_byBin_min         = arrayfun(@(x) max([min(amps(bin == x)),data.thresholds_microV(2)]), 1:Nbins); % mean by phase
        data.(L).AMP_byBin_max         = arrayfun(@(x) max([amps(bin == x);data.thresholds_microV(2)]), 1:Nbins); % mean by phase
        
        
        data.(L).AMP_byBin(isnan(data.(L).AMP_byBin))=replace_nans;
        
        
        % 6. Put results for real data together
        %data.(L).spike_phases_radians    = spikephases;
        %data.(L).spike_phases_histogram  = hist(data.(L).spike_phases_radians, cfg.phase.phase_bin_centers) / length(interval_starts) / length(cfg.phase.phase_bin_centers); % compute overal phase histogram
        %data.(L).spike_phases_histogram = real_histogram;
        %data.(L).spike_phases_histogram_smoothed = hist3([spikephases, cycleNums_withSpikes], 'ctrs', {cfg.phase.phase_bin_centers 1:length(interval_starts)});
        
        
        % compute shuffled intervals
        shuffled_ends   = Triggers.(E).shuffled_ts;
        shuffled_starts = Triggers.(E).shuffled_ts  - Triggers.(E).shuffled_intervals;
        % for reshuffled data - calculate heart cycle phase where individual spikes ended up
        shuffled_SDF         = nan(size(shuffled_starts,1),Nbins);
        for s = 1:size(shuffled_starts,1)
            shuffled_within_trial_idx=any(bsxfun(@ge,shuffled_starts(s,:)',trial_starts) & bsxfun(@le,shuffled_ends(s,:)',trial_ends),2)';
            starts = shuffled_starts(s,shuffled_within_trial_idx)';
            ends = shuffled_ends(s,shuffled_within_trial_idx)';
            durations=ends-starts;
            [events2include,cycleNums_withSpikes] = find(bsxfun(@ge,spiketimes,starts') & bsxfun(@le,spiketimes,ends'));
            ts = spiketimes(events2include); % include only events that landed within any interval
            spikephases = 2*pi*(ts - starts(cycleNums_withSpikes))./durations(cycleNums_withSpikes);
            %eventsTaken = single(find(events2include))';
            % all data
            shuffled_histogram = hist3([spikephases, cycleNums_withSpikes], 'ctrs', {phase_bins 1:length(starts)});
            %shuffled_SDF(s,:) = mean(average_smooth_data(shuffled_histogram,cfg));
            shuffled_SDF(s,:) = sum(shuffled_histogram.*repmat(numel(phase_bins)./durations',numel(phase_bins),1),2)/numel(durations);
        end
        SD        = ecg_bna_do_statistics(real_SDF',shuffled_SDF,1:Nbins);
        
        % 12. compute correlation between features phase dynamics and
        % spike dynamics
        [cc, pp]             = corrcoef(SD.SD_mean-SD.SDPmean, data.(L).AMP_byBin);
        data.(L).FR_amp_cc   = cc(2,1);
        data.(L).FR_amp_p    = pp(2,1);
        
        
        % all data
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
        data.(L).SDsubtractedSDP            = data.(L).SD - data.(L).SDP; % spikes/s, difference between mean and jittered data
        data.(L).SDsubtractedSDP_normalized = data.(L).SDsubtractedSDP ./ data.(L).SDP *100; % percent signal change
        data.(L).FR_ModIndex_SubtrSDP       = max(data.(L).SDsubtractedSDP) - min(data.(L).SDsubtractedSDP); % difference between max and min FR
        data.(L).FR_ModIndex_PcS            = max(data.(L).SDsubtractedSDP_normalized) - min(data.(L).SDsubtractedSDP_normalized); % difference between max and min % signal change
        
        
        % 11. estimate how much spike amplitude is far from the
        % thresholds
        % bin from the closest threshold first
        higher_threshold_detected = all(amps > data.thresholds_microV(1));
        data.(L).AMP_voltageBinned = hist(amps, 1:250); % bin amplitudes
        data.(L).AMP_all_high = higher_threshold_detected; % bin amplitudes
        
        
        %         if higher_threshold_detected % both thresholds are lower than all the spike amplitudes
        %             % take the higher one - this indicates to me that
        %             data.(L).AMP_voltageBins   = data.thresholds_microV(1):500;
        %             data.(L).AMP_voltageBinned = histc(amps, data.(L).AMP_voltageBins); % bin amplitudes
        %         else
        %             data.(L).AMP_voltageBins   = data.thresholds_microV(2):500;
        %             data.(L).AMP_voltageBinned = histc(amps, data.(L).AMP_voltageBins); % bin amplitudes
        %         end
        
        first_nonzero_bin     = find(data.(L).AMP_voltageBinned,1,'first');
        data.(L).distance2thr = first_nonzero_bin-data.thresholds_microV(2); % basically, number of empty bins from the spike with the lowest amplitude and the closest threshold
        
        
        %         %% dont like this condition at all
        %         % check the number of spikes left after computing phases
        %         if length(spikephases) < 3 || length(lowIBI_eventPhases) < 3 || length(highIBI_eventPhases) < 3 || ...
        %                 length(interval_starts) < 100 || length(data.(L).lowIBI_timeRRstart) < 3 || length(data.(L).highIBI_timeRRstart) < 3
        %             continue
        %         end
        %
        
        %% IF we are fitting, we want to use shuffled data too!
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
        
        
    end
    save([cfg.cardioballistic_folder filesep data.unitId '_' data.target '__spikes_ECGphase.mat'], 'data')
    clear data
end
end

function real_data = average_smooth_data(input, cfg)
Kernel = normpdf(-5*cfg.time.gaussian_kernel:cfg.time.PSTH_binwidth:5*cfg.time.gaussian_kernel,0,cfg.time.gaussian_kernel); % build kernel
tmp = [ [zeros(1,80); input(:,1:end-1)'] ...
    input' ...
    [input(:,2:end)'; zeros(1,80)] ];
conv_data = conv2(tmp,Kernel,'same');
real_data = conv_data(:,Nbins+1:end-Nbins);
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


