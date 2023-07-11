function ecg_bna_plots_per_session( data, con_info,cfg, PlotMethod,varargin )


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
    elseif strcmp(baseline_method, 'none')
    cbtitle = 'not-Normalized';
    imscale = [-1, 1];
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
    plottitle = [data.site_ID ', Target ' data.target ' (' injection '): ' con_info(cn).label];
    
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
        if isempty(cat(3, con_data(:, hs).(PlotMethod).pow.mean)) % this is a strange break condition to be honest
            continue;
        end
        % concatenate tfs for different state windows for plotting
        concat.pow = [];
        concat.itpc = [];
        concat.itpcbp = [];
        concat.lfp = [];
        concat.pow_sgnf = [];
        concat.itpc_sgnf = [];
        concat.itpcbp_sgnf = [];
        concat.lfp_sgnf = [];
        concat.tfr_time = [];
        concat.lfp_time = [];
        
        concat.freq = con_data(1, hs).(PlotMethod).freq; %% this is actually in the settings...
        concat.label = con_data(1, hs).hs_label; %% this is in the settings too   
        
        state_info = struct();
        for st = 1:size(con_data, 1)
            % state timing information
            % state onset sample number
            con=con_data(st, hs).(PlotMethod);
            con.pow_sgnf=con_data(st, hs).significance.pow;
            con.itpc_sgnf=con_data(st, hs).significance.itpc;
            con.itpcbp_sgnf=con_data(st, hs).significance.itpcbp;
            con.lfp_sgnf=con_data(st, hs).significance.lfp;
            con.shuffled=con_data(1, hs).shuffled;
            
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
%             %% i think this part is for differences plots... not sure what to do here with
%             if plot_significant && isfield(con, 'stat_test') && ~isempty(con.stat_test.h)
%                 con.stat_test.h(isnan(con.stat_test.h))=0;
%                 powspctrm = powspctrm .* con.stat_test.h;
%             end
            
            % concatenate across states with a NaN separation in between
            concat.pow      = cat(3, concat.pow,   con.pow.mean,   nan(size(con.pow.mean, 1),   size(con.pow.mean, 2),   100/25));
            concat.itpc     = cat(3, concat.itpc,  con.itpc.mean, nan(size(con.itpc.mean, 1), size(con.itpc.mean, 2), 100/25));
            concat.itpcbp   = cat(3, concat.itpcbp,con.itpcbp.mean,  nan(size(con.itpcbp.mean, 1),    size(con.itpcbp.mean, 2),    100));
            concat.lfp      = cat(2, concat.lfp,   con.lfp.mean,     nan(size(con.lfp.mean, 1),        100));
            concat.tfr_time = [concat.tfr_time, con.tfr_time, nan(1, 100/25)];
            concat.lfp_time = [concat.lfp_time, con.time,     nan(1, 100)];
            concat.pow_sgnf = cat(3, concat.pow_sgnf,   con.pow_sgnf,   nan(size(con.pow_sgnf, 1),   size(con.pow_sgnf, 2),   100/25));
            concat.itpc_sgnf = cat(3, concat.itpc_sgnf,   con.itpc_sgnf,   nan(size(con.itpc_sgnf, 1),   size(con.itpc_sgnf, 2),   100/25));
            concat.itpcbp_sgnf = cat(3, concat.itpcbp_sgnf,   con.itpcbp_sgnf,   nan(size(con.itpcbp_sgnf, 1),   size(con.itpcbp_sgnf, 2),   100));
            concat.lfp_sgnf = cat(2, concat.lfp_sgnf,   con.lfp_sgnf,     nan(size(con.lfp_sgnf, 1),        100));
                       
            % somehow needed for (not) labelling not existing alignments
            if ~all(isnan(con.tfr_time))
                states_valid=[states_valid st];
            end
        end
        
        %% plot
        
        state_onsets = find(concat.tfr_time == 0);
        states_names={con_data(states_valid, hs).(PlotMethod).state_name};
        state_samples = sort([state_info.start_s, state_info.onset_s, state_info.finish_s]);
        
        subplottitle =concat.label{1};
        if isfield(con_data(1, hs), 'nsessions')
            subplottitle = [subplottitle ' (nsessions = ' num2str(con_data(1, hs).nsessions) ')'];
        elseif isfield(con_data(1, hs), 'nsites')
            subplottitle = [subplottitle ' (nsites = ' num2str(con_data(1, hs).nsites) ')'];
        elseif isfield(con_data(1, hs), 'ntrials') && ~isempty(con_data(1, hs).ntrials)
            subplottitle = [subplottitle ' (ntrials = ' num2str(con_data(1, hs).ntrials) ')'];
        end
        

        %% POW and ITPC
        toplot={concat.pow,concat.itpc};
        sigplot = {concat.pow_sgnf, concat.itpc_sgnf};
        for figr=1:2 % frequency spectra
            figure(h(figr));
            sph{figr}(hs)=subplot(nhandlabels, nspacelabels, hs);
            image(1:size(toplot{figr},3), 1:numel(concat.freq), squeeze(toplot{figr}),'CDataMapping','scaled');
            hold on;
            if strcmp(PlotMethod,'real')
                % calculate the significance here:
                significance = double(squeeze(sigplot{figr})); %+repmat(randi([0 1], 60,19),[1,11]); checking the plot results
                significance(significance==0)=NaN;
                % plotting a contour around the significant parts:
                contour(1:size(toplot{figr},3),1:numel(con.freq),significance,1,'linecolor','k')
            end
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
            
            clear significance;
            
        end
        
        %% Evoked LFP
        % Smoothing of the evoked LFP here:
        jnk = [];
        win = 1:cfg.smoothWin; win=win-(numel(win)+1)/2;
        %gauss=normpdf(win,0,numel(win)/6);
        for k=1:size(concat.lfp,1)
            jnk(k,:)=conv(concat.lfp(k,:), normpdf(win,0,numel(win)/6),'same'); %./max(conv(ones(100,1), gausswin(win)));
        end
        con_lfp_mean_smooth = jnk;%(:,win:size(con.lfp.mean,2)-win); % cutting the zero padding part of the conv, from begin and end of the results
        
        figure(h(4));
        subplot(nhandlabels, nspacelabels, hs);
        hold on;
%         plot(con.time(win:size(con.lfp.mean,2)-win), con_lfp_mean_smooth) %, 'Color', colors(i,:));
        plot(concat.lfp_time', squeeze(con_lfp_mean_smooth)')
        line([0 0], ylim, 'color', 'k');
        title(subplottitle,'Fontsize',10);
        xlabel('Time(s)');
        
        if strcmp(PlotMethod,'real')
            lineprops={};
            shadedErrorBar(con.time, con.shuffled.lfp.mean,con.shuffled.lfp.std,lineprops,1);
            
            clear significance
            %t = con.time(win:size(con.lfp.mean,2)-win);
            ylm = get(gca,'Ylim');
            significance = double(squeeze(concat.lfp_sgnf));
            significance(significance==0)=NaN;
            significance=significance.*ylm(1);
            % adding the signifiance horizontal lines:
            plot(repmat(concat.lfp_time,size(concat.itpcbp,2),1)', significance','linewidth',3);
        end
        %% Bandpassed ITPC        
        % Smoothing of the itpcbp here:
        jnk = [];
        win = cfg.smoothWin;
        for m = 1: size(concat.itpcbp,1)
            for k=1:size(concat.itpcbp,2)
                jnk(m,k,:)= conv(squeeze(concat.itpcbp(m,k,:)), gausswin(win),'same');%./max(conv(ones(100,1), gausswin(win)));
            end
        end
        con_itpcbp_smooth = jnk;%(:,:,win:size(concat.itpcbp,3)-win); % cutting the zero padding part of the conv, from begin and end of the results
        % see above, why now we don't need to cut
        
        figure(h(3));
        subplot(nhandlabels, nspacelabels, hs);
        set(gca,'ColorOrder',jet(size(concat.itpcbp,2)));               %% change color order to something nicer
        plot(repmat(concat.lfp_time,size(concat.itpcbp,2),1)', squeeze(con_itpcbp_smooth)')
        hold on;
        
        if strcmp(PlotMethod,'real')
            % adding the signifiance horizontal lines:
            clear significance
            %t = con.time(win:size(con.lfp.mean,2)-win);
            ylm = get(gca,'Ylim');
            %ytck = get(gca,'Ytick');
            stp = (ylm(2)-ylm(1))/100;
            significance = double(squeeze(concat.itpcbp_sgnf));         % i needed to create concat.itpcbp_sgnf, it basically appends Nans for a (potential) separator with a second alignment
            significance(significance==0)=NaN;                          % replacing zeros with Nans means once we plot, lines will be discontinoous there
            multiplicator=(1:size(significance,1))*-1*stp;              % multiplicator basically defines position of significance line
            significance=significance.*repmat(multiplicator',1,size(significance,2));
            plot(repmat(concat.lfp_time,size(concat.itpcbp,2),1)', significance','linewidth',3);
        end
        legend(strcat(num2str(round(cfg.tfr.frequency_bands(:,1))), '-',num2str(round(cfg.tfr.frequency_bands(:,2))), ' Hz'));
 
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
        if nargin > 4
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
        results_file = fullfile(fldr, [plot_names{figr} '_' data.site_ID '_' con_info(cn).label ' ' PlotMethod]);   
        mtit([plottitle ' ' PlotMethod],'interpreter','none')
        export_fig(h(figr), results_file, '-pdf');
    end
    close all
end
end

