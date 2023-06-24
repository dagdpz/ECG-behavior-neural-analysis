function ecg_bna_plots_per_session( data, con_info,cfg, varargin )


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
    imscale = [0 1e-8];
end

% number of subplots required
nhandlabels = length(cfg.compare.reach_hands);
nspacelabels = length(cfg.compare.reach_spaces);


%% loop through conditions
for cn= 1:numel(data.condition)    
    if isempty(fieldnames(data.condition(cn))) || isempty(fieldnames(data.condition(cn).state_hs)) % &&~isempty([sites_data(i).condition(cn).state_hs.lfp])
        continue
    end    
    if con_info(cn).perturbation == 0
        injection = 'Pre';
    elseif con_info(cn).perturbation == 1
        injection = 'Post';
    else
        injection = 'Any';
    end
    plottitle = ['LFP (' injection '): Site ' data.site_ID ', Target ' data.target '(ref_' cfg.ref_hemisphere '), '  con_info(cn).label];
    if con_info(cn).choice == 0
        plottitle = [plottitle 'Instructed trials'];
    elseif con_info(cn).choice == 1
        plottitle = [plottitle 'Choice trials'];
    end
    
    % create figures with handles
    plot_names={'POW','ITPC','ITPC_BP','LFP_Evoked'};
    results_folders={[cfg.sites_lfp_fldr filesep 'Rpeak_evoked_TFS'],[cfg.sites_lfp_fldr filesep 'Rpeak_evoked_PHA'],[cfg.sites_lfp_fldr filesep 'Rpeak_evoked_PHABP'],[cfg.sites_lfp_fldr filesep 'Rpeak_evoked_LFP']};
    for f=1:numel(plot_names)
        h(f) = figure('Name',plot_names{f});
        %set(h(f), 'position', [100, 100,900, 675]);
        hold on
    end
    
    con_data=data.condition(cn).state_hs;
    states_valid=[];
    collim{1}=[];
    collim{2}=[];
    % loop through handspace
    for hs = 1:size(con_data, 2)
        %freq =  cat(2, con_data(:, hs).freq);
        if isempty(cat(3, con_data(:, hs).pow.mean)) % this is a strange break condition to be honest
            continue;
        end
        % concatenate tfs for different state windows for plotting
        concat.pow = [];
        concat.itpc = [];
        concat.itpcbp = [];
        concat.lfp = [];
        concat.tfr_time = [];
        concat.lfp_time = [];
        
        concat.freq = con_data(1, hs).freq; %% this is actually in the settings...
        concat.label = con_data(1, hs).hs_label; %% this is in the settings too        
        
        state_info = struct();
        for st = 1:size(con_data, 1)
            % state timing information
            % state onset sample number
            con=con_data(st, hs);
            state_info(st).onset_s = find(con.tfr_time <= 0, 1, 'last');
            % state onset time
            state_info(st).onset_t = 0;
            % start start sample
            state_info(st).start_s = 1;
            % state start time
            state_info(st).start_t = con.tfr_time(1);
            % start finish sample
            state_info(st).finish_s = length(con.tfr_time);
            % start end sample
            state_info(st).finish_t = con.tfr_time(end);
            
            % state onset, start and finish samples for further states offset from previous state window
            if st > 1
                state_info(st).start_s  = length(concat.state_time) + state_info(st).start_s;
                state_info(st).finish_s = length(concat.state_time) + state_info(st).finish_s;
                state_info(st).onset_s  = length(concat.state_time) + state_info(st).onset_s;
            end
            
            
            con.pow.mean(isnan(con.pow.mean))=0;
            con.itpc.mean(isnan(con.itpc.mean))=0;
            con.itpcbp.mean(isnan(con.itpcbp.mean))=0;
%             phasesBP = nanmean(con.itpcbp.mean, 1);
%             phasesBP = squeeze(phasesBP);
            %% i think this part is for differences plots... not sure what to do here with
            if plot_significant && isfield(con, 'stat_test') && ~isempty(con.stat_test.h)
                con.stat_test.h(isnan(con.stat_test.h))=0;
                powspctrm = powspctrm .* con.stat_test.h;
            end
            
            % concatenate across states with a NaN separation in between
            concat.pow      = cat(3, concat.pow,   con.pow.mean,   nan(size(con.pow.mean, 1),   size(con.pow.mean, 2),   100/25));
            concat.itpc     = cat(3, concat.itpc,  con.itpc.mean, nan(size(con.itpc.mean, 1), size(con.itpc.mean, 2), 100/25));
            concat.itpcbp   = cat(3, concat.itpcbp,con.itpcbp.mean,    nan(size(con.itpcbp.mean, 1),    size(con.itpcbp.mean, 2),    100));
            concat.lfp      = cat(2, concat.lfp,   con.lfp.mean,     nan(size(con.lfp.mean, 1),        100));
            concat.tfr_time = [concat.tfr_time, con.tfr_time, nan(1, 100/25)];
            concat.lfp_time = [concat.lfp_time, con.time,     nan(1, 100)];
            
            % somehow needed for (not) labelling not existing alignments
            if ~all(isnan(con.tfr_time))
                states_valid=[states_valid st];
            end
            
%             powspctrm = nanmean(con.pow.mean, 1);
%             phasespctrm = nanmean(con.itpc.mean, 1);
            
            
        end
        
        %% plot
        
        state_onsets = find(concat.tfr_time == 0);
        states_names={con_data(states_valid, hs).state_name};
        state_samples = sort([state_info.start_s, state_info.onset_s, state_info.finish_s]);
        
        subplottitle = concat.label{1};
        if isfield(con_data(1, hs), 'nsessions')
            subplottitle = [subplottitle ' (nsessions = ' num2str(con_data(1, hs).nsessions) ')'];
        elseif isfield(con_data(1, hs), 'nsites')
            subplottitle = [subplottitle ' (nsites = ' num2str(con_data(1, hs).nsites) ')'];
        elseif isfield(con_data(1, hs), 'ntrials') && ~isempty(con_data(1, hs).ntrials)
            subplottitle = [subplottitle ' (ntrials = ' num2str(con_data(1, hs).ntrials) ')'];
        end
        
        toplot={concat.pow,concat.itpc};
        for figr=1:2 % frequency spectra
            figure(h(figr));
            sph{figr}(hs)=subplot(nhandlabels, nspacelabels, hs);
%             imagesc(1:size(toplot{figr},3), con.freq, squeeze(toplot{figr}));
            image(1:size(toplot{figr},3), 1:numel(con.freq), squeeze(toplot{figr}),'CDataMapping','scaled');
            hold on;

            nonnan=toplot{figr};nonnan(isnan(nonnan))=[];
            collim{figr}=[min([collim{figr}(:); nonnan(:)]) max([collim{figr}(:); nonnan(:)])];

            % horizontal lines to separate frequency bands
            fbandstart = [2, 4, 8, 12, 18, 32, 80];
            fbandstart_idx = zeros(size(fbandstart));
            for f = fbandstart
                f_idx = find(abs(con.freq - f) == min(abs(con.freq - f)), 1, 'first');
                line([state_info(st).start_s state_info(st).finish_s], [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                fbandstart_idx(fbandstart == f) = f_idx;
            end
            
            set(gca,'TickDir','out')
            % log y axis ticks
            %set(gca, 'ytick', ([1:8:numel(concat_states_tfs.freq)]));
            set(gca, 'ytick', fbandstart_idx);
            set(gca, 'yticklabel', fbandstart);
            % add 0.5 at end since the time value is the center of the bin
            % add 0 at beginning to make x-axis visible
            set(gca, 'ylim', [0.5,numel(con.freq) + 0.5]);
            %round(concat_states_tfs.freq([1:8:numel(concat_states_tfs.freq)])));
            for so = state_onsets
                line([so so], ylim, 'color', 'k');
                if isfield(con_data(state_onsets == so, hs), 'state_name') && ~isempty(states_names(state_onsets == so))
                    state_name = states_names{state_onsets == so};
                    text(so+1, 10, state_name, 'fontsize', 8);
                end
            end
            
            % mark state onsets
            state_ticks=round(concat.tfr_time(state_samples), 1);
            set(gca,'xtick',state_samples(~isnan(state_ticks)))
            set(gca,'xticklabels', state_ticks(~isnan(state_ticks)), 'fontsize', 8)
            set(gca, 'xticklabelrotation', 45)
            % add 0.5 since the time value is the center of the bin
            % add 0 at the beginning to make the y-axis visible
            set(gca, 'xlim', [0 state_samples(end) + 0.5]);
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            
            title(subplottitle);
            
            
        end
        
        figure(h(4));
        subplot(nhandlabels, nspacelabels, hs);
        hold on;
        plot(con.time, con.lfp.mean) %, 'Color', colors(i,:));
        lineprops={};
        shadedErrorBar(con.time, con.lfp_shuff.mean,con.lfp_shuff.std,lineprops,1);
        line([0 0], ylim, 'color', 'k');
        title(subplottitle);        
        xlabel('Time(s)');
        
        figure(h(3));
        subplot(nhandlabels, nspacelabels, hs);
        plot(repmat(concat.lfp_time,size(concat.itpcbp,2),1)', squeeze(concat.itpcbp)')
        legend(strcat(num2str(round(cfg.tfr.frequency_bands(:,1))), '-',num2str(round(cfg.tfr.frequency_bands(:,2))), ' Hz'));
        hold on;
        
        
        %concat_states_tfs.time = 1:1:size(concat_states_tfs.powspctrm, 3);
        
        
        
%         
%         for figr=1:2 % frequency spectra
%             figure(h(figr));
%             subplot(nhandlabels, nspacelabels, hs);
%             %subplot(nhandlabels, nspacelabels, hs)
%             %imagesc(concat_states_tfs.time, [1:numel(concat_states_tfs.freq)], squeeze(concat_states_tfs.powspctrm), [-1 1]);
%             %axis xy, 
%             
%             %change aspect ratio if only 2 conditions
% %             if size(avg_tfr, 2) < 3
% %                 set(gca,'DataAspectRatio', [1 0.6 1]);
% %             end
%         end
%         
%         % evoked LFP
%         figure(h(3));
%         hold on
%         subplot(nhandlabels, nspacelabels, hs);
%         %ylabel(yaxislabel);
    end
    
    %% format spectra colors
    for figr=1:2
        figure(h(figr));
        
        %set(gcf,'CLim',collim{figr})
        for hs = 1:size(con_data, 2)
            subplot(sph{figr}(hs));
            set(gca,'CLim',collim{figr})
            %clim(collim{figr});
        end
        
        cm = colormap('jet');
        if nargin > 3
            cm = colormap(varargin{1});
        end
            cb = colorbar;
            set(get(cb,'title'),'string', cbtitle, 'fontsize',8);
        colormap(cm);
    end
    
    for figr=1:numel(h)
        fldr=results_folders{figr};
        [~,PRTS]=fileparts(fldr);
        if ~exist(fldr,'dir')
           mkdir(cfg.sites_lfp_fldr,PRTS);
        end
        results_file = fullfile(fldr, [plot_names '_' data.site_ID '_' con_info(cn).label]);   
        mtit(plottitle)
        export_fig(h(figr), results_file, '-pdf');
    end
    close all
end
end

