function ecg_bna_plots_per_session( sites_data, site_conditions,cfg, varargin )


%lfp_tfa_plot_hs_tuned_tfr_multiple_img  - Plots the LFP time frequency spectrogram
%averages for different hand-space conditions to be compared
%
% USAGE:
%   lfp_tfa_plot_hs_tuned_tfr_multiple_img( avg_tfr, cfg, plottitle, results_file )
%   lfp_tfa_plot_hs_tuned_tfr_multiple_img( avg_tfr, cfg, plottitle, results_file, cm )
%   lfp_tfa_plot_hs_tuned_tfr_multiple_img( avg_tfr, cfg, plottitle, results_file, cm, plot_significant )
%
%
% INPUTS:
%       avg_tfr         - average LFP time frequency response for different
%       hand-space conditions to be compared
%		cfg     - struct containing the required settings
%           Required Fields: see settings/lfp_tfa_settings_example
%               1. baseline_method             - method used for baseline
%               normalization
%               2. compare.reach_hands          - hand labels to compare
%               3. compare.reach_spaces         - space labels to compare
%       plottitle       - title for the plot
%       results_file    - path to filename to store the resulting image
%       varargin        - colormap to be used (default = 'jet', can be any
%                       standard colormap additionally supported is 'bluewhitered')
%                       - flag to indicate if only significant difference
%                       bins should be plotted
%
% REQUIRES:	bluewhitered, export_fig
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_plot_site_average_tfr,
% lfp_tfa_avg_tfr_across_sessions, lfp_tfa_avg_tfr_across_sites,
% bluewhitered, colormap, lfp_tfa_compute_difference_condition_tfr
%
% Author(s):	S.Nair, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2019-02-15:	Created function (Sarath Nair)
% 2019-03-05:	First Revision
% ...
% $Revision: 1.0 $  $Date: 2019-03-05 17:18:00 $

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%

plot_significant = 0;
if nargin > 4
    plot_significant = varargin{2};
end

% colorbar title
if strcmp(cfg.baseline_method, 'zscore')
    cbtitle = 'Z-score';
    imscale = [-1, 1];
elseif strcmp(cfg.baseline_method, 'division')
    cbtitle = 'P / \mu';
    imscale = [0, 2];
elseif strcmp(cfg.baseline_method, 'subtraction')
    cbtitle = 'P - \mu';
    imscale = [0 1e-8];
elseif strcmp(cfg.baseline_method, 'relchange')
    cbtitle = '(P - \mu) / \mu';
    imscale = [-1, 1];
end

% number of subplots required
nhandlabels = length(cfg.compare.reach_hands);
nspacelabels = length(cfg.compare.reach_spaces);


%% loop through conditions
for cn= 1:numel(sites_data.condition)    
    if isempty(fieldnames(sites_data.condition(cn))) || isempty(fieldnames(sites_data.condition(cn).state_hs)) % &&~isempty([sites_data(i).condition(cn).state_hs.lfp])
        continue
    end    
    if site_conditions(cn).perturbation == 0
        injection = 'Pre';
    elseif site_conditions(cn).perturbation == 1
        injection = 'Post';
    else
        injection = 'Any';
    end
    plottitle = ['LFP (' injection '): Site ' sites_data.site_ID ', Target ' sites_data.target '(ref_' cfg.ref_hemisphere '), '  site_conditions(cn).label];
    if site_conditions(cn).choice == 0
        plottitle = [plottitle 'Instructed trials'];
    elseif site_conditions(cn).choice == 1
        plottitle = [plottitle 'Choice trials'];
    end
    
    % create 3 figures with handles
    plot_names={'LFP_TFR','LFP_PHA','LFP_PBP','LFP_Evoked'};
    results_folders={[cfg.sites_lfp_fldr filesep 'Rpeak_evoked_TFS'],[cfg.sites_lfp_fldr filesep 'Rpeak_evoked_PHA'],[cfg.sites_lfp_fldr filesep 'Rpeak_evoked_PHABP'],[cfg.sites_lfp_fldr filesep 'Rpeak_evoked_LFP']};
    for f=1:numel(plot_names)
        h(f) = figure;
        set(h, 'position', [100, 100,900, 675]);
        hold on
    end
    
    avg_tfr=sites_data.condition(cn).state_hs;
    states_valid=[];
    % loop through handspace
    for hs = 1:size(avg_tfr, 2)
        freq =  cat(2, avg_tfr(:, hs).freq);
        if isempty(cat(3, freq.powspctrm)) % this is a strange break condition to be honest
            continue;
        end
        % concatenate tfs for different state windows for plotting
        concat_states_tfs.powspctrm = [];
        concat_states_tfs.phasespctrm = [];
        concat_states_tfs.state_time = [];
        concat_states_tfs.freq = avg_tfr(1, hs).freq.freq; %% this is actually in the settings...
        concat_states_tfs.label = avg_tfr(1, hs).hs_label; %% this is in the settings too        
        
        state_info = struct();
        for st = 1:size(avg_tfr, 1)
            % state timing information
            % state onset sample number
            state_info(st).onset_s = find(avg_tfr(st, hs).freq.time <= 0, 1, 'last');
            % state onset time
            state_info(st).onset_t = 0;
            % start start sample
            state_info(st).start_s = 1;
            % state start time
            state_info(st).start_t = avg_tfr(st, hs).freq.time(1);
            % start finish sample
            state_info(st).finish_s = length(avg_tfr(st, hs).freq.time);
            % start end sample
            state_info(st).finish_t = avg_tfr(st, hs).freq.time(end);
            
            % state onset, start and finish samples for further states offset from previous state window
            if st > 1
                state_info(st).start_s  = length(concat_states_tfs.state_time) + state_info(st).start_s;
                state_info(st).finish_s = length(concat_states_tfs.state_time) + state_info(st).finish_s;
                state_info(st).onset_s  = length(concat_states_tfs.state_time) + state_info(st).onset_s;
            end
            % concatenate across states with a NaN separation in between
            concat_states_tfs.powspctrm     = cat(3, concat_states_tfs.powspctrm,       avg_tfr(st, hs).freq.powspctrm,   nan(size(avg_tfr(st, hs).freq.powspctrm, 1), length(concat_states_tfs.freq), 100/25));
            concat_states_tfs.phasespctrm   = cat(3, concat_states_tfs.phasespctrm,     avg_tfr(st, hs).freq.phasespctrm, nan(size(avg_tfr(st, hs).freq.phasespctrm, 1), length(concat_states_tfs.freq), 100/25));
            concat_states_tfs.state_time    = [concat_states_tfs.state_time, avg_tfr(st, hs).freq.time, nan(1, 100/25)];
            
            % somehow needed for (not) labelling not existing alignments
            if ~all(isnan(avg_tfr(st, hs).freq.time))
                states_valid=[states_valid st];
            end
            
            %% i think this part is for differences plots... not sure what to do here with
            avg_tfr(st, hs).freq.powspctrm(isnan(avg_tfr(st, hs).freq.powspctrm))=0;
            avg_tfr(st, hs).freq.phasespctrm(isnan(avg_tfr(st, hs).freq.phasespctrm))=0;
            avg_tfr(st, hs).freq.phasesBP(isnan(avg_tfr(st, hs).freq.phasesBP))=0;
            state_lfp_powspctrm = nanmean(avg_tfr(st, hs).freq.powspctrm, 1);
            state_lfp_phasespctrm = nanmean(avg_tfr(st, hs).freq.phasespctrm, 1);
            state_lfp_phasesBP = nanmean(avg_tfr(st, hs).freq.phasesBP, 1);
            state_lfp_phasesBP = squeeze(state_lfp_phasesBP);
            if plot_significant && isfield(avg_tfr(st, hs).freq, 'stat_test') && ~isempty(avg_tfr(st, hs).freq.stat_test.h)
                avg_tfr(st, hs).freq.stat_test.h(isnan(avg_tfr(st, hs).freq.stat_test.h))=0;
                state_lfp_powspctrm = state_lfp_powspctrm .* avg_tfr(st, hs).freq.stat_test.h;
            end
            toplot={state_lfp_powspctrm,state_lfp_phasespctrm};
            for figr=1:2 % frequency spectra
                figure(h(figr));
                subplot(nhandlabels, nspacelabels, hs);
                hold on;
                imagesc(...
                    linspace(state_info(st).start_s, state_info(st).finish_s, ...
                    state_info(st).finish_s - state_info(st).start_s + 1), ...
                    linspace(1, length(avg_tfr(st, hs).freq.freq), ...
                    length(avg_tfr(st, hs).freq.freq)), ...
                    squeeze(toplot{figr}) , imscale);
                
                % horizontal lines to separate frequency bands
                fbandstart = [2, 4, 8, 12, 18, 32, 80];
                fbandstart_idx = zeros(size(fbandstart));
                for f = fbandstart
                    f_idx = find(abs(avg_tfr(st, hs).freq.freq - f) == min(abs(avg_tfr(st, hs).freq.freq - f)), 1, 'first');
                    line([state_info(st).start_s state_info(st).finish_s], [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                    fbandstart_idx(fbandstart == f) = f_idx;
                end
            end
            
            figure(h(4));
            subplot(nhandlabels, nspacelabels, hs);
            hold on;
            plot(avg_tfr.time, avg_tfr.mean) %, 'Color', colors(i,:));
            
            %             if ploterr
            %                 for i = 1:size(evoked_lfp(1, hs).mean, 1)
            %                     plot(evoked_lfp(1, hs).time, evoked_lfp(1, hs).mean(i, :) + evoked_lfp(1, hs).std(i, :), [colors(mod(i, length(colors))) ':']);
            %                     plot(evoked_lfp(1, hs).time, evoked_lfp(1, hs).mean(i, :) - evoked_lfp(1, hs).std(i, :), [colors(mod(i, length(colors))) ':']);
            %                 end
            %             end
            
            %             if isfield(evoked_lfp(1, hs), 'shuffled_mean') && isfield(evoked_lfp(1, hs), 'shuffled_std')
            %                 for i = 1:size(evoked_lfp(1, hs).shuffled_mean, 1)
%             plot(avg_tfr.time, avg_tfr.shuffled_mean) %, 'Color', colors(i),'linestyle','-.');
%             plot(avg_tfr.time, avg_tfr.shuffled_mean + avg_tfr.shuffled_std) %, [colors(mod(i, length(colors))) ':']);
%             plot(avg_tfr.time, avg_tfr.shuffled_mean - avg_tfr.shuffled_std) %, [colors(mod(i, length(colors))) ':']);
            lineprops={};
            shadedErrorBar(avg_tfr.time, avg_tfr.shuffled_mean,avg_tfr.shuffled_std,lineprops,1);
            %                 end
            %             end
            figure(h(3));
            subplot(nhandlabels, nspacelabels, hs);
            plot(avg_tfr.freq.time, state_lfp_phasesBP)
            legend(strcat(num2str(round(cfg.tfr.frequency_bands(:,1))), '-',num2str(round(cfg.tfr.frequency_bands(:,2))), ' Hz'));
            hold on;
            
            
        end
        %concat_states_tfs.time = 1:1:size(concat_states_tfs.powspctrm, 3);
        
        %% format figure
        state_onsets = find(concat_states_tfs.state_time == 0);
        states_names={avg_tfr(states_valid, hs).state_name};
        state_samples = sort([state_info.start_s, state_info.onset_s, state_info.finish_s]);
        
        subplottitle = concat_states_tfs.label{1};
        if isfield(avg_tfr(1, hs), 'nsessions')
            subplottitle = [subplottitle ' (nsessions = ' num2str(avg_tfr(1, hs).nsessions) ')'];
        elseif isfield(avg_tfr(1, hs), 'nsites')
            subplottitle = [subplottitle ' (nsites = ' num2str(avg_tfr(1, hs).nsites) ')'];
        elseif isfield(avg_tfr(1, hs), 'ntrials') && ~isempty(avg_tfr(1, hs).ntrials)
            subplottitle = [subplottitle ' (ntrials = ' num2str(avg_tfr(1, hs).ntrials) ')'];
        end
        
        
        for figr=1:2 % frequency spectra
            figure(h(figr));
            subplot(nhandlabels, nspacelabels, hs);
            %subplot(nhandlabels, nspacelabels, hs)
            %imagesc(concat_states_tfs.time, [1:numel(concat_states_tfs.freq)], squeeze(concat_states_tfs.powspctrm), [-1 1]);
            %axis xy, 
            cb = colorbar;
            set(get(cb,'title'),'string', cbtitle, 'fontsize',8);
            set(gca,'TickDir','out')
            % log y axis ticks
            %set(gca, 'ytick', ([1:8:numel(concat_states_tfs.freq)]));
            set(gca, 'ytick', (fbandstart_idx));
            set(gca, 'yticklabel', fbandstart);
            %round(concat_states_tfs.freq([1:8:numel(concat_states_tfs.freq)])));
            for so = state_onsets
                line([so so], ylim, 'color', 'k');
                if isfield(avg_tfr(state_onsets == so, hs), 'state_name') && ~isempty(states_names(state_onsets == so))
                    state_name = states_names{state_onsets == so};
                    text(so+1, 10, state_name, 'fontsize', 8);
                end
            end
            
            % add 0.5 at end since the time value is the center of the bin
            % add 0 at beginning to make x-axis visible
            set(gca, 'ylim', [0 numel(avg_tfr(st, hs).freq.freq) + 0.5]);
            % mark state onsets
            state_ticks=round(concat_states_tfs.state_time(state_samples), 1);
            set(gca,'xtick',state_samples(~isnan(state_ticks)))
            set(gca,'xticklabels', state_ticks(~isnan(state_ticks)), 'fontsize', 8)
            set(gca, 'xticklabelrotation', 45)
            % add 0.5 since the time value is the center of the bin
            % add 0 at the beginning to make the y-axis visible
            set(gca, 'xlim', [0 state_samples(end) + 0.5]);
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            
            title(subplottitle);
            
            %change aspect ratio if only 2 conditions
%             if size(avg_tfr, 2) < 3
%                 set(gca,'DataAspectRatio', [1 0.6 1]);
%             end
        end
        
        % evoked LFP
        figure(h(3));
        hold on
        subplot(nhandlabels, nspacelabels, hs);
        line([0 0], ylim, 'color', 'k');
        title(subplottitle);
        
        xlabel('Time(s)');
        %ylabel(yaxislabel);
    end
    
    for figr=1:2
        figure(h(figr));
        cm = colormap('jet'); 
        if nargin > 3
            cm = colormap(varargin{1});
            colorbar;
        end
        colormap(cm);
    end
    
    for figr=1:numel(h)
        fldr=results_folders{figr};
        [~,PRTS]=fileparts(fldr);
        if ~exist(fldr,'dir')
           mkdir(cfg.sites_lfp_fldr,PRTS);
        end
        results_file = fullfile(fldr, [plot_names '_' sites_data.site_ID '_' site_conditions(cn).label]);   
        mtit(plottitle)
        export_fig(h(figr), results_file, '-pdf');
    end
    close all
end
end

