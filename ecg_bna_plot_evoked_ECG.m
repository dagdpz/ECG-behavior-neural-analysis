function ecg_bna_plot_evoked_ECG( evoked_lfp, ecg_bna_cfg, plottitle, results_file, varargin )
%ecg_bna_plot_evoked_lfp  - Plots the LFP/ECG evoked response
%averages for different conditions to be compared
%
% USAGE:
%   ecg_bna_plot_evoked_lfp( evoked_lfp, ecg_bna_cfg, plottitle, results_file )
%   ecg_bna_plot_evoked_lfp( evoked_lfp, ecg_bna_cfg, plottitle,
%   results_file, 'ylabel', 'ECG amplitude' )
%
% INPUTS:
%       evoked_lfp       - average LFP power spectrum for different
%       hand-space conditions to be compared
%		ecg_bna_cfg      - struct containing the required settings
%           Required Fields: see ecg_bna_settings
%               compare.reach_hands     - reach hands
%               compare.reach_spaces    - reach spaces
%       plottitle        - title for the plot
%       results_file     - path to filename to store the resulting image
%       varargin         - other settings for plotting given as name-value
%       pairs. Currently only supports name 'ylabel' (if a user specified
%       label has to be assigned for y axis, by default y axis label is
%       'LFP amplitude')
%
% See also ecg_bna_compute_session_Rpeak_evoked_LFP,
% ecg_bna_avg_sessions_Rpeak_evoked_LFP, ecg_bna_compute_session_evoked_ECG,
% ecg_bna_avg_sessions_Rpeak_evoked_ECG
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
nhandlabels = 1; nspacelabels = 1;
if isfield(ecg_bna_cfg.compare, 'reach_hands')
    nhandlabels = length(ecg_bna_cfg.compare.reach_hands);
end
if isfield(ecg_bna_cfg.compare, 'reach_spaces')
    nspacelabels = length(ecg_bna_cfg.compare.reach_spaces);
end

n_conditions=nhandlabels*nspacelabels;
n_states=size(evoked_lfp, 1);

% loop through handspace
for st=1:size(evoked_lfp, 1)
    for hs = 1:n_conditions
        if ~isempty([evoked_lfp(st,hs).mean]) &&  ~all(all(isnan([evoked_lfp(st,hs).mean]))) && ~isempty([evoked_lfp(st,hs).std])
            
            % now plot
            subplot(n_conditions, n_states, (hs-1)*n_states+st)
            hold on;
            colors = ['b'; 'r'; 'g'; 'y'; 'm'; 'c'; 'k'];
            if isfield(evoked_lfp(st, hs), 'colors') && ~isempty(evoked_lfp(st, hs).colors)
                colors = evoked_lfp(st, hs).colors;
            end
            for i = 1:size(evoked_lfp(st, hs).mean, 1)
                plot(evoked_lfp(st, hs).time, evoked_lfp(st, hs).mean(i, :), 'Color', colors(i,:));
            end
            if isfield(evoked_lfp(st, hs), 'legend')
                legend(evoked_lfp(st, hs).legend);
            end
            if ploterr
                for i = 1:size(evoked_lfp(st, hs).mean, 1)
                    plot(evoked_lfp(st, hs).time, evoked_lfp(st, hs).mean(i, :) + evoked_lfp(st, hs).std(i, :), [colors(mod(i, length(colors))) ':']);
                    plot(evoked_lfp(st, hs).time, evoked_lfp(st, hs).mean(i, :) - evoked_lfp(st, hs).std(i, :), [colors(mod(i, length(colors))) ':']);
                end
            end
            
            if isfield(evoked_lfp(st, hs), 'shuffled_mean') && isfield(evoked_lfp(st, hs), 'shuffled_std')
                for i = 1:size(evoked_lfp(st, hs).shuffled_mean, 1)
                    plot(evoked_lfp(st, hs).time, evoked_lfp(st, hs).shuffled_mean(i, :), 'Color', colors(i),'linestyle','-.');
                    plot(evoked_lfp(st, hs).time, evoked_lfp(st, hs).shuffled_mean(i, :) + evoked_lfp(st, hs).shuffled_std(i, :), [colors(mod(i, length(colors))) ':']);
                    plot(evoked_lfp(st, hs).time, evoked_lfp(st, hs).shuffled_mean(i, :) - evoked_lfp(st, hs).shuffled_std(i, :), [colors(mod(i, length(colors))) ':']);
                end
            end
            
            % mark state onsets
            %if isfield(evoked_lfp, 'state_name')
            %set(gca,'xtick', unique(state_samples))
            line([0 0], ylim, 'color', 'k');
            
            
            subplottitle = [evoked_lfp(st, hs).hs_label{1}];
            if isfield(evoked_lfp(1, hs), 'nsessions')
                subplottitle = [subplottitle ' (nsessions = ' num2str(evoked_lfp(1, hs).nsessions) ')'];
            end
            if isfield(evoked_lfp(1, hs), 'nsites')
                subplottitle = [subplottitle ' (nsites = ' num2str(evoked_lfp(1, hs).nsites) ')'];
            end
            if isfield(evoked_lfp(1, hs), 'ntrials') && ~isempty(evoked_lfp(1, hs).ntrials)
                subplottitle = [subplottitle ' (ntrials = ' num2str(evoked_lfp(1, hs).ntrials) ')'];
            end
            if isfield(evoked_lfp(1, hs), 'npeaks') && ~isempty(evoked_lfp(1, hs).npeaks)
                subplottitle = [subplottitle ' (npeaks = ' num2str(evoked_lfp(1, hs).npeaks) ')'];
            end
            if isfield(evoked_lfp(1, hs), 'nshuffles') && ~isempty(evoked_lfp(1, hs).nshuffles)
                subplottitle = [subplottitle ' (nshuffles = ' num2str(evoked_lfp(1, hs).nshuffles) ')'];
            end
            
            if st==1
                ylabel({subplottitle;yaxislabel});
            end
            if hs==1
                title(ecg_bna_cfg.analyse_Rpeak_states{st,2});
            end
            if hs==n_conditions
                xlabel('Time(s)');
            end
            
        end
    end
end
ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_'), 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
export_fig(h, results_file,'-pdf');

end

