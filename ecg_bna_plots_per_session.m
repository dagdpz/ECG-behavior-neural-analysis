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
if strcmp(cfg.shuffle_normalization_method, 'zscore') % can change to baseline_method
    % 'zscore' - P_norm(t,f) = ( mean(real) - mean(shuffled) ) / std(shuffled)
    %     cbtitle = 'Z-score'; % for baseline normalization
    cbtitle = '(P - \mu) / std'; % for shuffle predictor normalization
    %     imscale = [-1, 1];
elseif strcmp(cfg.shuffle_normalization_method, 'division')
    cbtitle = 'P / \mu';
elseif strcmp(cfg.shuffle_normalization_method, 'subtraction')
    cbtitle = 'P - \mu';
elseif strcmp(cfg.shuffle_normalization_method, 'relchange')
    cbtitle = '(P - \mu) / \mu';
elseif strcmp(cfg.shuffle_normalization_method, 'none')
    cbtitle = 'not-Normalized';
end

% number of subplots required
nhandlabels = length(cfg.compare.reach_hands);
nspacelabels = length(cfg.compare.reach_spaces);
plot_names={'POW','ITPC','Power_BP','ITPC_BP','LFP_Evoked'};
nsubplots=numel(plot_names);

nrows=ceil(sqrt((nsubplots*nhandlabels*nspacelabels)));
ncolumns=2; %nrows;
results_folder=[cfg.sites_lfp_fldr filesep];

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
    
    % create figure
    h = figure('units','normalized','position',[0 0 1 1]);
    
    con_data=data.condition(cn).state_hs;
    states_valid=[];
    collim{1}=[];
    collim{2}=[];
    % loop through handspace
    for hs = 1:size(con_data, 2)
        if isempty(cat(3, con_data(:, hs).(PlotMethod).pow.mean)) % this is a strange break condition to be honest
            continue;
        end
        % concatenate tfs for different state windows for plotting
        concat.pow = [];
        concat.itpc = [];
        concat.itpcbp = [];
        concat.powbp = [];
        concat.lfp = [];
        concat.pow_sgnf = [];
        concat.itpc_sgnf = [];
        concat.itpcbp_sgnf = [];
        concat.powbp_sgnf = [];
        concat.lfp_sgnf = [];
        concat.tfr_time = [];
        concat.lfp_time = [];
        
        concat.freq  = con_data(1, hs).(PlotMethod).freq; %% this is actually in the settings...
        concat.label = con_data(1, hs).hs_label; %% this is in the settings too
        
        state_info = struct();
        for st = 1:size(con_data, 1)
            % state timing information
            % state onset sample number
            con=con_data(st, hs).(PlotMethod);
            con.pow_sgnf=con_data(st, hs).significance.pow;
            con.itpc_sgnf=con_data(st, hs).significance.itpc;
            con.itpcbp_sgnf=con_data(st, hs).significance.itpcbp;
            con.powbp_sgnf =con_data(st, hs).significance.powbp;
            con.lfp_sgnf=con_data(st, hs).significance.lfp;
            con.shuffled=con_data(st, hs).shuffled;
%             if ~isempty(con_data(st, hs).shuffled)
%                 con.shuffled=con_data(st, hs).shuffled;
%             else
%                 con.shuffled = [];
%             end
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
            con.powbp.mean(isnan(con.itpcbp.mean))=0;
            
            % concatenate across states with a NaN separation in between
            concat.pow      = cat(3, concat.pow,   con.pow.mean,   nan(size(con.pow.mean, 1),   size(con.pow.mean, 2),   100/25));
            concat.itpc     = cat(3, concat.itpc,  con.itpc.mean, nan(size(con.itpc.mean, 1), size(con.itpc.mean, 2), 100/25));
            concat.itpcbp   = cat(3, concat.itpcbp,con.itpcbp.mean,  nan(size(con.itpcbp.mean, 1),    size(con.itpcbp.mean, 2),    100));
            concat.powbp    = cat(3, concat.powbp,con.powbp.mean,  nan(size(con.powbp.mean, 1),    size(con.powbp.mean, 2),    100));
            concat.lfp      = cat(2, concat.lfp,   con.lfp.mean,     nan(size(con.lfp.mean, 1),        100));
            concat.tfr_time = [concat.tfr_time, con.tfr_time, nan(1, 100/25)];
            concat.lfp_time = [concat.lfp_time, con.time,     nan(1, 100)];
            concat.pow_sgnf = cat(3, concat.pow_sgnf,   con.pow_sgnf,   nan(size(con.pow_sgnf, 1),   size(con.pow_sgnf, 2),   100/25));
            concat.itpc_sgnf = cat(3, concat.itpc_sgnf,   con.itpc_sgnf,   nan(size(con.itpc_sgnf, 1),   size(con.itpc_sgnf, 2),   100/25));
            concat.itpcbp_sgnf = cat(3, concat.itpcbp_sgnf,   con.itpcbp_sgnf,   nan(size(con.itpcbp_sgnf, 1),   size(con.itpcbp_sgnf, 2),   100));
            concat.powbp_sgnf = cat(3, concat.powbp_sgnf,   con.powbp_sgnf,   nan(size(con.powbp_sgnf, 1),   size(con.powbp_sgnf, 2),   100));
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
        state_nRpeaks = (con_data(states_valid, hs).real.nRpeaks);
        
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
        for sp=1:2 % frequency spectra
            sph(sp,hs)=subplot(nrows, ncolumns, (nsubplots)*(hs-1)+sp);
            image(1:size(toplot{sp},3), 1:numel(concat.freq), squeeze(toplot{sp}),'CDataMapping','scaled');
            set(gca,'YDir','normal');
            hold on;
            % calculate the significance here:
            significance = double(squeeze(sigplot{sp})); %+repmat(randi([0 1], 60,19),[1,11]); checking the plot results
            %significance(significance==0)=NaN;
            % plotting a contour around the significant parts:
            contour(1:size(toplot{sp},3),1:numel(con.freq),significance,1,'linecolor','k')
            nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
            collim{sp}=[min([collim{sp}(:); nonnan(:)]) max([collim{sp}(:); nonnan(:)])];
            
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
            set(gca,'xticklabelrotation', 45)
            % add 0.5 since the time value is the center of the bin
            % add 0 at the beginning to make the y-axis visible
            set(gca, 'xlim', [0 state_samples(end) + 0.5]);
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            title([plot_names{sp},' - ',subplottitle, ' - nShuffles = ',num2str(cfg.n_permutations),' - nRpeaks = ',num2str(state_nRpeaks)]);
            
            clear significance;
            
        end
        
        %% Bandpassed POWER
        % Smoothing of the itpcbp here:
        jnk = [];
        win = 1:cfg.smoothWin; win=win-(numel(win)+1)/2;
        half_win = ceil(size(win,2)/2)-1;
        concat_powbp = cat(3,(concat.powbp(:,:,half_win:-1:1)),(concat.powbp(:,:,:)));
        concat_powbp = cat(3,concat_powbp, (concat.powbp(:,:,end:-1:end-half_win+1)));
        for m = 1: size(concat.powbp,1)
            for k=1:size(concat.powbp,2)
                jnk(m,k,:)= conv(squeeze(concat_powbp(m,k,:)), normpdf(win,0,numel(win)/6),'same');%./max(conv(ones(100,1), gausswin(win)));
            end
        end
        clear concat_powbp
        con_powbp_smooth = jnk(:,:,half_win+1:end-half_win);
        
        
        sp=3;
        sph(sp,hs)=subplot(nrows, ncolumns, (nsubplots)*(hs-1)+sp);             %% change color order to something nicer
        hold on;
        set(gca,'ColorOrder',jet(size(concat.powbp,2)));
        plot(repmat(concat.lfp_time,size(concat.powbp,2),1)', squeeze(con_powbp_smooth)')
        line([0 0], ylim, 'color', 'k');
        xlabel('Time(s)'); ylabel('Power (W)');

        %         if strcmp(PlotMethod,'real')
        % adding the signifiance horizontal lines:
        clear significance
        ylm = get(gca,'Ylim');
        stp = (ylm(2)-ylm(1))/100;
        significance = double(squeeze(concat.powbp_sgnf));         % i needed to create concat.itpcbp_sgnf, it basically appends Nans for a (potential) separator with a second alignment
        significance(significance==0)=NaN;                          % replacing zeros with Nans means once we plot, lines will be discontinoous there
        multiplicator= ylm(1)+(1:size(significance,1))*-1*stp;              % multiplicator basically defines position of significance line
        significance=significance.*repmat(multiplicator',1,size(significance,2));
        plot(repmat(concat.lfp_time,size(concat.powbp,2),1)', significance','linewidth',3);
        %         end
        legend({strcat(num2str(round(cfg.tfr.frequency_bands(:,1))), '-',num2str(round(cfg.tfr.frequency_bands(:,2))), ' Hz')},'fontsize',3);
        title([plot_names{sp},' - ',subplottitle, ' - nShuffles = ',num2str(cfg.n_permutations),...
            ' - nRpeaks = ',num2str(state_nRpeaks)],'fontsize',9,'interpreter','none');
        
        %% Bandpassed ITPC
        % Smoothing of the itpcbp here:
        jnk = [];
        win = 1:cfg.smoothWin; win=win-(numel(win)+1)/2;
        half_win = ceil(size(win,2)/2)-1;
        concat_itpcbp = cat(3,(concat.itpcbp(:,:,half_win:-1:1)),(concat.itpcbp(:,:,:)));
        concat_itpcbp = cat(3,concat_itpcbp, (concat.itpcbp(:,:,end:-1:end-half_win+1)));
        for m = 1: size(concat.itpcbp,1)
            for k=1:size(concat.itpcbp,2)
                jnk(m,k,:)= conv(squeeze(concat_itpcbp(m,k,:)), normpdf(win,0,numel(win)/6),'same');%./max(conv(ones(100,1), gausswin(win)));
            end
        end
        clear concat_itpcbp
        con_itpcbp_smooth = jnk(:,:,half_win+1:end-half_win);
        
        sp=4;
        sph(sp,hs)=subplot(nrows, ncolumns, (nsubplots)*(hs-1)+sp);             %% change color order to something nicer
        hold on;
        set(gca,'ColorOrder',jet(size(concat.itpcbp,2)));
        plot(repmat(concat.lfp_time,size(concat.itpcbp,2),1)', squeeze(con_itpcbp_smooth)')
        line([0 0], ylim, 'color', 'k');
        xlabel('Time(s)'); ylabel('ITPC value');
        
        %         if strcmp(PlotMethod,'real')
        % adding the signifiance horizontal lines:
        clear significance
        ylm = get(gca,'Ylim');
        stp = (ylm(2)-ylm(1))/100;
        significance = double(squeeze(concat.itpcbp_sgnf));         % i needed to create concat.itpcbp_sgnf, it basically appends Nans for a (potential) separator with a second alignment
        significance(significance==0)=NaN;                          % replacing zeros with Nans means once we plot, lines will be discontinoous there
        multiplicator= ylm(1)+(1:size(significance,1))*-1*stp;              % multiplicator basically defines position of significance line
        significance=significance.*repmat(multiplicator',1,size(significance,2));
        plot(repmat(concat.lfp_time,size(concat.itpcbp,2),1)', significance','linewidth',3);
        %         end
        legend({strcat(num2str(round(cfg.tfr.frequency_bands(:,1))), '-',num2str(round(cfg.tfr.frequency_bands(:,2))), ' Hz')},'fontsize',3);
        title([plot_names{sp},' - ',subplottitle, ' - nShuffles = ',num2str(cfg.n_permutations),...
            ' - nRpeaks = ',num2str(state_nRpeaks)],'fontsize',9,'interpreter','none');
        
        %% Evoked LFP
        % Smoothing of the evoked LFP here:
        jnk = [];
        win = 1:cfg.smoothWin; win=win-(numel(win)+1)/2;
        half_win = ceil(size(win,2)/2)-1;
        for k=1:size(concat.lfp,1)
            jnk(k,:)=conv([concat.lfp(k,half_win:-1:1),concat.lfp(k,:),concat.lfp(k,end:-1:end-half_win+1)], normpdf(win,0,numel(win)/6),'same');
        end
        con_lfp_mean_smooth = jnk(:,half_win+1:end-half_win);
        
        sp=5;
        sph(sp,hs)=subplot(nrows, ncolumns, (nsubplots)*(hs-1)+sp);
        hold on;
        plot(concat.lfp_time', squeeze(con_lfp_mean_smooth)','linewidth',1.5)
        line([0 0], ylim, 'color', 'k');
        xlabel('Time(s)'); ylabel('Voltage (V)');
        title([plot_names{sp},' - ',subplottitle, ' - nShuffles = ',num2str(cfg.n_permutations),...
            ' - nRpeaks = ',num2str(state_nRpeaks)],'fontsize',9,'interpreter','none');
        
        if strcmp(PlotMethod,'real')
            lineprops={};
            %             shadedErrorBar(con.time, con.lfp.mean,con.lfp.std./(size(con.lfp.std,2).^0.5),lineprops,1);
            shadedErrorBar(con.time, con.shuffled.lfp.mean,con.shuffled.lfp.std,lineprops,1);
        end
        clear significance
        ylm = get(gca,'Ylim');
        significance = double(squeeze(concat.lfp_sgnf));
        significance(significance==0)=NaN;
        significance=significance.*ylm(1);
        % adding the signifiance horizontal lines:
        plot(repmat(concat.lfp_time,size(concat.itpcbp,2),1)', significance','linewidth',3);
        %         end
        
    end
    
    %% format spectra colors
    for sp=1:2
        for hs = 1:size(con_data, 2)
            subplot(sph(sp,hs));
            set(gca,'CLim',collim{sp})
        end
        cm = colormap('jet');
        if nargin > 4
            cm = colormap(varargin{1});
        end
        cb = colorbar;
        if strcmp(PlotMethod,'normalized')
            set(get(cb,'title'),'string', cbtitle, 'fontsize',8);
        else
            if sp==1
                set(get(cb,'title'),'string', 'power(W)', 'fontsize',8);
            elseif sp==2
                set(get(cb,'title'),'string', 'ITPC', 'fontsize',8);
            end
        end
        colormap(cm);
    end
    
    results_file = fullfile(results_folder, [data.site_ID '_' con_info(cn).label ' ' PlotMethod]);
    if strcmp(PlotMethod,'normalized')
        mtit([plottitle ' ' PlotMethod ' (' cfg.shuffle_normalization_method ')'],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
    else
        mtit([plottitle ' ' PlotMethod],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
    end
    export_fig(h, results_file, '-pdf');
    close all
end
end

