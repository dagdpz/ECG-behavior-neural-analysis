function ecg_bna_plot_evoked_R2Rt( evoked_R2Rt, ecg_bna_cfg, plottitle, results_file, varargin )
%ecg_bna_plot_evoked_lfp  - Plots the ECG R2R interval averages for given
%time windows and conditions to be compared 
%
% USAGE:
%   ecg_bna_plot_evoked_R2Rt( evoked_R2Rt, ecg_bna_cfg, plottitle, results_file )
%   ecg_bna_plot_evoked_R2Rt( evoked_R2Rt, ecg_bna_cfg, plottitle,
%   results_file, 'ylabel', 'Norm. R2Rt', 'err', 'stdev' )  
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
%       pairs. Names can be 
%           'ylabel'  - if a user specified label has to be assigned for y
%           axis, value must be a string) 
%           'err'     - error measure to be used. Can be 'stdev' for
%           standard deviation(default) or 'stderr' for standard error. 
%
% See also ecg_bna_compute_session_evoked_ECG_R2Rt,
% ecg_bna_avg_sessions_ECGb2bt_evoked

    % defaults
    if ecg_bna_cfg.normalize_R2Rt
        yaxislabel = 'Relative R2R time';
    else
        yaxislabel = 'R2R time (s)';
    end
    err = 'stdev'; % standard deviation
    
    % get settings
    if nargin > 4
        settings = struct(varargin{:});
        if isfield(settings, 'ylabel')
            yaxislabel = settings.ylabel;
        end
        if isfield(settings, 'err')
            err = settings.err;
        end
        
    end
    
    h1 = figure;
    
    %set(h1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    set(h1, 'position', [100, 100,900, 675]);
    
    % number of offset samples to divide between time windows
    noffset = 150;
    
    % number of subplots required
    nhandlabels = 1; nspacelabels = 1;
    if isfield(ecg_bna_cfg.compare, 'reach_hands')
        nhandlabels = length(ecg_bna_cfg.compare.reach_hands);
    end
    if isfield(ecg_bna_cfg.compare, 'reach_spaces')
        nspacelabels = length(ecg_bna_cfg.compare.reach_spaces);
    end
    
    % loop through handspace
    for hs = 1:size(evoked_R2Rt, 2)
        if ~isempty([evoked_R2Rt(:,hs).mean]) &&  ~isempty([evoked_R2Rt(:,hs).std])
            % concatenate states
            concat_states_R2Rt = struct();
            concat_states_R2Rt.trial = {};
            concat_states_R2Rt.mean = [];
            concat_states_R2Rt.err = [];
            concat_states_R2Rt.time = [];

            state_info = struct();             
            
            nstates = size(evoked_R2Rt, 1);
            ax = cell(1, numel(evoked_R2Rt));
            for st = 1:nstates
                
                if isempty(evoked_R2Rt(st,hs).mean) && isempty(evoked_R2Rt(st,hs).std)
                    continue;
                end
                
                state_info(st).onset_s = find(...
                    evoked_R2Rt(st, hs).time <= 0, 1, 'last');
                state_info(st).onset_t = 0;
                state_info(st).start_s = 1;
                state_info(st).start_t = evoked_R2Rt(st, hs).time(1);
                state_info(st).finish_s = length(evoked_R2Rt(st, hs).time);
                state_info(st).finish_t = evoked_R2Rt(st, hs).time(end);                    

                if st > 1
                    state_info(st).start_s = length(concat_states_R2Rt.time) + ...
                        state_info(st).start_s;
                    state_info(st).finish_s = length(concat_states_R2Rt.time) + ...
                        state_info(st).finish_s;
                    state_info(st).onset_s = length(concat_states_R2Rt.time) + ...
                        state_info(st).onset_s;
                end

                % concatenate mean, std and time of evoked R2Rt for
                % different states
                if isfield(evoked_R2Rt(1, hs), 'ntrials') && ...
                     ~isempty(evoked_R2Rt(1, hs).ntrials)
                    state_evoked_R2Rt = evoked_R2Rt(st, hs).ecg_b2bt;
                    if ecg_bna_cfg.normalize_R2Rt
                           state_evoked_R2Rt = state_evoked_R2Rt ./ ...
                               repmat(evoked_R2Rt(st, hs).trials_mean', ...
                               [1 size(state_evoked_R2Rt, 2)]);
                    end
                    concat_states_R2Rt.trial = [concat_states_R2Rt.trial, ...
                    horzcat(state_evoked_R2Rt, nan(size(evoked_R2Rt(st, hs).ecg_b2bt, 1), noffset))];
                end
                
                state_evoked_mean = evoked_R2Rt(st, hs).mean;
                state_evoked_err = evoked_R2Rt(st, hs).std;
                if strcmp(err, 'stderr')
                    state_evoked_err = evoked_R2Rt(st, hs).std ...
                        / sqrt(size(evoked_R2Rt(st, hs).ecg_b2bt, 1));
                elseif strcmp(err, 'stdev')
                    state_evoked_err = evoked_R2Rt(st, hs).std;
                end                
                
                concat_states_R2Rt.mean = [concat_states_R2Rt.mean, ...
                    state_evoked_mean, nan(size(evoked_R2Rt(st, hs).mean, 1), noffset)];
                concat_states_R2Rt.err = [concat_states_R2Rt.err, ...
                    state_evoked_err, nan(size(evoked_R2Rt(st, hs).std, 1), noffset)];
                concat_states_R2Rt.time = [concat_states_R2Rt.time, ...
                    evoked_R2Rt(st, hs).time, nan(1, noffset)];
                concat_states_R2Rt.label = evoked_R2Rt(st, hs).hs_label;
                if isfield(evoked_R2Rt(st, hs), 'legend')
                    concat_states_R2Rt.legend = evoked_R2Rt(st, hs).legend;
                end
                
                % now plot
                ax{st} = subplot(nhandlabels*nspacelabels*2, nstates, (((hs-1)*nstates) + st));
                hold on;
                colors = ['b'; 'r'; 'g'; 'y'; 'm'; 'c'; 'k'];

                % plot individual trials
                if isfield(evoked_R2Rt(st, hs), 'ntrials') && ~isempty(evoked_R2Rt(st, hs).ntrials)
                    xx = [0];
                    %for s = 1:length(evoked_R2Rt(st, hs).ecg_b2bt)
                    xx = linspace(1, size(evoked_R2Rt(st, hs).ecg_b2bt, 2), ...
                        size(evoked_R2Rt(st, hs).ecg_b2bt, 2)) + xx(end);
                    state_evoked_R2Rt = evoked_R2Rt(st, hs).ecg_b2bt;
                    if ecg_bna_cfg.normalize_R2Rt
                           state_evoked_R2Rt = state_evoked_R2Rt ./ ...
                               repmat(evoked_R2Rt(st, hs).trials_mean', ...
                               [1 size(state_evoked_R2Rt, 2)]);
                    end
                    plot(ax{st}, evoked_R2Rt(st, hs).time, state_evoked_R2Rt, 'Color', [0.6, 0.6, 0.6])
                    %end
                end
                for i = 1:size(evoked_R2Rt(st, hs).mean, 1)
                    shadedErrorBar(evoked_R2Rt(st, hs).time, evoked_R2Rt(st, hs).mean(i, :), evoked_R2Rt(st, hs).std(i, :), 'lineprops', colors(i, :));
                end
                if isfield(evoked_R2Rt(st, hs), 'legend')
                    legend(ax{st}, evoked_R2Rt(st, hs).legend);
                end


                line([0 0], ylim, 'color', 'k'); 
                if isfield(evoked_R2Rt(st, hs), 'state_name') && ...
                        ~isempty(evoked_R2Rt(st, hs).state_name)
                    state_name = evoked_R2Rt(st, hs).state_name;
                    plottxt = state_name;
                    if isfield(evoked_R2Rt(st, hs), 'ntrials') && ...
                        ~isempty(evoked_R2Rt(st, hs).ntrials)
                        ntrials = (evoked_R2Rt(st, hs).ntrials);
                        plottxt = sprintf('%s (%g)', plottxt, ntrials);
                    end
                    if isfield(evoked_R2Rt(st, hs), 'nsessions')
                        nsessions = (evoked_R2Rt(st, hs).nsessions);
                        plottxt = sprintf('%s (nsessions = %g)', plottxt, nsessions);
                    end
                    ypos = ylim;
                    ypos = ypos(1) + (ypos(2) - ypos(1))*0.2;
                    title(plottxt);
                end

                
                xlabel(ax{st}, 'Time(s)');
                ylabel(ax{st}, yaxislabel);

            end
            
            linkaxes([ax{:}],'y');
            
        end
    end
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    export_fig(h1, [results_file '.png']);
    print(h1, '-depsc', [results_file '.ai']);

end

