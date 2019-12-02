function ecg_bna_plot_evoked_lfp( evoked_lfp, ecg_bna_cfg, plottitle, results_file, varargin )
%ecg_bna_plot_evoked_lfp  - Plots the LFP evoked response
%averages for different hand-space conditions to be compared
%
% USAGE:
%   ecg_bna_plot_evoked_lfp( evoked_lfp, ecg_bna_cfg, plottitle, results_file )
%
% INPUTS:
%       evoked_lfp       - average LFP power spectrum for different
%       hand-space conditions to be compared
%		ecg_bna_cfg      - struct containing the required settings
%           Required Fields: see ecg_bna_settings
%               1. 
%       plottitle        - title for the plot
%       results_file     - path to filename to store the resulting image
%
% REQUIRES:	
%
% See also ecg_bna_settings, ecg_bna_plot_site_evoked_LFP, 
% ecg_bna_avg_evoked_LFP_across_sites, 
% ecg_bna_avg_evoked_LFP_across_sessions
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

    h = figure;
    set(h, 'position', [100, 100,900, 675]);
    %set(h, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    
    yaxislabel = 'LFP amplitude';
    if nargin > 4
        settings = struct(varargin{:});
        yaxislabel = settings.ylabel;
    end
    
    ploterr = 1;

    % number of offset samples to divide between time windows
    % noffset = 100;
    
    % number of subplots required
    nhandlabels = length(ecg_bna_cfg.compare.reach_hands);
    nspacelabels = length(ecg_bna_cfg.compare.reach_spaces);
    
    % loop through handspace
    for hs = 1:size(evoked_lfp, 2)
        if ~isempty([evoked_lfp(:,hs).mean]) &&  ~isempty([evoked_lfp(:,hs).std])
            
            % now plot
            subplot(nhandlabels, nspacelabels, hs)
            hold on;
            colors = ['b'; 'r'; 'g'; 'y'; 'm'; 'c'; 'k'];
            if isfield(evoked_lfp(1, hs), 'colors') && ...
                    ~isempty(evoked_lfp(1, hs).colors)
                colors = evoked_lfp(1, hs).colors;
            end
            for i = 1:size(evoked_lfp(1, hs).mean, 1)
                plot(evoked_lfp(1, hs).time, evoked_lfp(1, hs).mean(i, :), ...
                    'Color', colors(i,:));
            end
            if isfield(evoked_lfp(1, hs), 'legend')
                legend(evoked_lfp(1, hs).legend);
            end
            if ploterr
                for i = 1:size(evoked_lfp(1, hs).mean, 1)
                    plot(evoked_lfp(1, hs).time, evoked_lfp(1, hs).mean(i, :) + evoked_lfp(1, hs).std(i, :), [colors(mod(i, length(colors))) ':']);
                    plot(evoked_lfp(1, hs).time, evoked_lfp(1, hs).mean(i, :) - evoked_lfp(1, hs).std(i, :), [colors(mod(i, length(colors))) ':']);
                end
            end
            
            if isfield(evoked_lfp(1, hs), 'shuffled_mean') && ...
                    isfield(evoked_lfp(1, hs), 'shuffled_std')
                for i = 1:size(evoked_lfp(1, hs).shuffled_mean, 1)
                    plot(evoked_lfp(1, hs).time, ...
                        evoked_lfp(1, hs).shuffled_mean(i, :), ...
                        'Color', colors(i));
                    plot(evoked_lfp(1, hs).time, ...
                        evoked_lfp(1, hs).shuffled_mean(i, :) + evoked_lfp(1, hs).shuffled_std(i, :), ...
                        [colors(mod(i, length(colors))) ':']);
                    plot(evoked_lfp(1, hs).time, ...
                        evoked_lfp(1, hs).shuffled_mean(i, :) - evoked_lfp(1, hs).shuffled_std(i, :), ...
                        [colors(mod(i, length(colors))) ':']);
                end
            end
            
            % mark state onsets
            %if isfield(evoked_lfp, 'state_name')
            %set(gca,'xtick', unique(state_samples))
            line([0 0], ylim, 'color', 'k'); 
%             for so = state_onsets
%                 line([so so], ylim, 'color', 'k'); 
%                 if isfield(evoked_lfp(state_onsets == so, hs), 'state_name') && ...
%                         ~isempty(evoked_lfp(state_onsets == so, hs).state_name)
%                     state_name = evoked_lfp(state_onsets == so, hs).state_name;
%                     ypos = ylim;
%                     ypos = ypos(1) + (ypos(2) - ypos(1))*0.2;
%                     text(so+1, ypos, state_name, 'fontsize', 10);
%                 end
%             end
%             %end
%             set(gca,'xticklabels', round(concat_states_lfp.time(unique(state_samples)), 2), 'fontsize', 10)
%             set(gca, 'xticklabelrotation', 45);
            xlabel('Time(s)');
            ylabel(yaxislabel);
            
            subplottitle = [evoked_lfp(1, hs).hs_label{1}];
            if isfield(evoked_lfp(1, hs), 'nsessions')
                subplottitle = [subplottitle ' (nsessions = ' ...
                    num2str(evoked_lfp(1, hs).nsessions) ')'];
            end
            if isfield(evoked_lfp(1, hs), 'nsites')
                subplottitle = [subplottitle ' (nsites = ' ...
                    num2str(evoked_lfp(1, hs).nsites) ')'];
            end
            if isfield(evoked_lfp(1, hs), 'ntrials') && ...
                    ~isempty(evoked_lfp(1, hs).ntrials)
                subplottitle = [subplottitle ' (ntrials = ' ...
                    num2str(evoked_lfp(1, hs).ntrials) ')'];
            end
            if isfield(evoked_lfp(1, hs), 'npeaks') && ...
                    ~isempty(evoked_lfp(1, hs).npeaks)
                subplottitle = [subplottitle ' (npeaks = ' ...
                    num2str(evoked_lfp(1, hs).npeaks) ')'];            
            end
            if isfield(evoked_lfp(1, hs), 'nshuffles') && ...
                    ~isempty(evoked_lfp(1, hs).nshuffles)
                subplottitle = [subplottitle ' (nshuffles = ' ...
                    num2str(evoked_lfp(1, hs).nshuffles) ')'];            
            end
            title(subplottitle);
        end
    end
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    export_fig(h, results_file);

end

