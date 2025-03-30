function [ out_triggered_site_data] = ecg_bna_compute_cluster_sig_withinSites_between_cond( triggered_site_data, data_label, events, cfg )
% Assuming the data is loaded in the following variables:
% data_task: data for the task condition (freq x time x trials)
% data_rest: data for the rest condition (freq x time x trials)
% tfr_time: A vector representing the time points (1 x time)
% tfr_freq: A vector representing the frequency points (1 x freq)
path_to_triggered_site_data = fullfile(cfg.sites_lfp_fldr);
if ~exist(path_to_triggered_site_data, 'dir')
    mkdir(path_to_triggered_site_data);
end

cond = {[triggered_site_data.condition(1).label,'-',triggered_site_data.condition(2).label]};
% events = {triggered_site_data.condition(1).event.event_name};
for fi = 1:length(data_label)
    Fin = data_label{fi};
    for e = 1: length(events)
        if strcmp(triggered_site_data.condition(1).event(e).event_name , 'Cue')
            continue
        end
        if strcmp(Fin,'pow') ||strcmp(Fin,'itpc')
            data_rest = triggered_site_data.condition(1).event(e).real.(Fin).complete;
            data_task =  triggered_site_data.condition(2).event(e).real.(Fin).complete;
            % ================================================================
% % %             % removing the field here??
% % %             triggered_site_data.condition(1).event(e).real.(Fin) = rmfield(triggered_site_data.condition(1).event(e).real.(Fin), "complete");
% % %             triggered_site_data.condition(2).event(e).real.(Fin) = rmfield(triggered_site_data.condition(2).event(e).real.(Fin), "complete");
% % %             
            
            tfr_freq = cfg.lfp.foi; %logspace(log10(2),log10(120),60);
            tfr_time = triggered_site_data.condition(1).event(e).tfr_time;
            
            freq_task = [];
            freq_task.label     = {Fin};  % Dummy label (no real channels here)
            freq_task.dimord    = 'rpt_freq_time';  % Trials, frequencies, and time
            freq_task.freq      = tfr_freq;  % Frequency vector
            freq_task.time      = tfr_time;  % Time vector
            freq_task.powspctrm = data_task;  % [trials x frequencies x time]
            
            freq_rest = [];
            freq_rest.label     = {Fin};  % Dummy label (no real channels here)
            freq_rest.dimord    = 'rpt_freq_time';  % Trials, frequencies, and time
            freq_rest.freq      = tfr_freq;  % Frequency vector
            freq_rest.time      = tfr_time;  % Time vector
            freq_rest.powspctrm = data_rest;  % [trials x frequencies x time]
            
            
            cfg_perm = [];
            cfg_perm.method           = 'montecarlo';             % Monte Carlo permutation method
            cfg_perm.statistic        = 'ft_statfun_indepsamplesT'; % Independent samples t-test
            %     cfg_perm.statistic        = 'ft_statfun_depsamplesT'; % dependent samples t-test
            cfg_perm.correctm         = 'cluster';                % Cluster-based correction
            cfg_perm.clusteralpha     = 0.05;                     % Cluster-forming threshold
            cfg_perm.clusterstatistic = 'maxsum';                 % Max-sum clustering or 'maxsize', 'wcm' (weighted cluster mass), 'wcs' (weighted cluster size)
            cfg_perm.minnbchan        = 0;                        % No neighboring channels (no spatial information)
            cfg_perm.tail             = 0;                        % Two-tailed test
            cfg_perm.clustertail      = 0;                        % Two-tailed for clustering
            cfg_perm.alpha            = 0.05;                     % Alpha level for the permutation test
            cfg_perm.correcttail      = 'alpha';                  % FDR correction for two-tailed T-test
            cfg_perm.numrandomization = 1000;                     % Number of randomizations
            cfg_perm.neighbours            = [];
            
            % a dummy connectivity matrix for the frequencies
            num_freqs = length(tfr_freq);  % Number of frequency bins
            cfg_perm.connectivity = ones(num_freqs, num_freqs);  % Fully connected between frequencies
            
            nTaskTrials = size(freq_task.powspctrm,1);
            nRestTrials = size(freq_rest.powspctrm,1);
            
            % Design matrix (2 rows: condition & trial index)
            design = zeros(2, nTaskTrials + nRestTrials);
            
            % First row: Condition (1 = Task, 2 = Rest)
            design(1, 1:nTaskTrials) = 1;  % Task condition
            design(1, nTaskTrials+1:end) = 2;  % Rest condition
            
            % Second row: Trial index
            design(2, 1:nTaskTrials) = 1:nTaskTrials;
            design(2, nTaskTrials+1:end) = 1:nRestTrials;  % Separate indexing
            
            cfg_perm.design = design;
            cfg_perm.ivar   = 1;  % Independent variable (Task vs. Rest)
            cfg_perm.uvar   = []; % No subject pairing (independent samples)
            
            [stat_tfs] = ft_freqstatistics(cfg_perm, freq_task, freq_rest);
            
            triggered_site_data.condition(1).event(e).stats.(Fin).val = stat_tfs;
            triggered_site_data.condition(2).event(e).stats.(Fin).val = stat_tfs;
            
            triggered_site_data.condition(1).event(e).stats.(Fin).cfg_perm = cfg_perm;
            triggered_site_data.condition(2).event(e).stats.(Fin).cfg_perm = cfg_perm;
            
            
            
            %=================================================================
            % plotting the significance
            h(e) = figure;
            imagesc(stat_tfs.time, 1:numel(stat_tfs.freq),stat_tfs.stat);
            set(gca,'YDir','normal');
            %     set(gca, 'YScale', 'log');  % If you want to show the y-axis in logarithmic scale
            fbandstart = unique(cfg.lfp.frequency_bands(:))';
            set(gca,'Xlim',[-.25 .25]);
            line([0 0], ylim, 'color', 'k');
            % horizontal lines to separate frequency bands
            fbandstart_idx = zeros(size(fbandstart));
            for f = fbandstart
                f_idx = find(abs(stat_tfs.freq - f) == min(abs(stat_tfs.freq - f)), 1, 'first');
                yline(f_idx, 'color', 'k', 'linestyle', '--');
                fbandstart_idx(fbandstart == f) = f_idx;
            end
            set(gca,'TickDir','out')
            set(gca, 'ytick', fbandstart_idx);
            set(gca, 'yticklabel', fbandstart);
            set(gca, 'ylim', [0.5,numel(stat_tfs.freq) + 0.5]);
            hold on
            contour(stat_tfs.time, 1:numel(stat_tfs.freq), stat_tfs.mask, [0.5, 0.5], 'LineColor', 'k', 'LineWidth', 2); % Overlay significant clusters
            if isfield(stat_tfs,"posclusterslabelmat") || ~isempty(stat_tfs.posclusterslabelmat)
                contour(stat_tfs.time, 1:numel(stat_tfs.freq), stat_tfs.posclusterslabelmat, [0.5, 0.5], 'LineColor', 'k', 'LineWidth', 2); % Overlay significant clusters
            end
            % clabel(C, h);
            colormap(jet);  % Change the colormap if desired
            colorbar;  % Show colorbar
            title([triggered_site_data.site_ID,'-',triggered_site_data.target,'-',cond{1},'-',Fin,'- Significant difference - numPerm =',...
                num2str(cfg_perm.numrandomization)],'fontsize', 8,'Interpreter', 'none');
            
            
            results_file{e} = fullfile([path_to_triggered_site_data,filesep,triggered_site_data.site_ID,'-',triggered_site_data.target,'-',...
                triggered_site_data.condition(1).event(e).event_name,' triggered_',Fin,'_',cond{1},'_Significant difference']);
            export_fig(h(e),[results_file{e},'.pdf']);
            
            close all,
            
        elseif strcmp(Fin,'lfp')
            data_rest = triggered_site_data.condition(1).event(e).real.(Fin).complete;
            data_task =  triggered_site_data.condition(2).event(e).real.(Fin).complete;
            % ================================================================================
% %             % removing the field here??
% %             triggered_site_data.condition(1).event(e).real.(Fin) = rmfield(triggered_site_data.condition(2).event(e).real.(Fin), "complete");
% %             triggered_site_data.condition(2).event(e).real.(Fin) = rmfield(triggered_site_data.condition(2).event(e).real.(Fin), "complete");
% %             
            time = triggered_site_data.condition(1).event(e).time;
            
            data_task = squeeze(data_task); % Convert to (Trials × Time)
            
            [nTrials, ~] = size(data_task);
            
            LFP_task = [];
            LFP_task.label     = {'LFP'};  % Meaningful channel label
            LFP_task.fsample   = 1000; % Adjust as needed
            LFP_task.time      = repmat({time}, 1, nTrials); % Time should be a cell array (1 per trial)
            LFP_task.trial     = arrayfun(@(i) data_task(i,:), 1:nTrials, 'UniformOutput', false); % Convert each trial to cell
            LFP_task.trialinfo = ones(nTrials, 1); % Assign condition labels (1 = Task)
            
            data_rest = squeeze(data_rest); % Convert to (Trials × Time)
            
            [nTrials, ~] = size(data_rest);
            
            LFP_rest = [];
            LFP_rest.label     = {'LFP'};  % Meaningful channel label
            LFP_rest.fsample   = 1000; % Adjust as needed
            LFP_rest.time      = repmat({time}, 1, nTrials); % Time should be a cell array (1 per trial)
            LFP_rest.trial     = arrayfun(@(i) data_rest(i,:), 1:nTrials, 'UniformOutput', false); % Convert each trial to cell
            LFP_rest.trialinfo = ones(nTrials, 1); % Assign condition labels (1 = Task)
            
            
            cfg_tl = [];
            cfg_tl.keeptrials = 'yes'; % Keep individual trials
            LFP_task_erp = ft_timelockanalysis(cfg_tl, LFP_task);
            LFP_rest_erp = ft_timelockanalysis(cfg_tl, LFP_rest);
            
            
            cfg_perm = [];
            cfg_perm.method           = 'montecarlo';               % Monte Carlo permutation method
            cfg_perm.statistic        = 'indepsamplesT';            %'ft_statfun_indepsamplesT'; % Independent samples t-test
            cfg_perm.correctm         = 'cluster';                  % Cluster-based correction
            cfg_perm.clusteralpha     = 0.05;                       % Cluster-forming threshold
            cfg_perm.clusterstatistic = 'maxsum';                   % Max-sum clustering or 'maxsize', 'wcm' (weighted cluster mass), 'wcs' (weighted cluster size)
            cfg_perm.minnbchan        = 0;                          % No neighboring channels (no spatial information)
            cfg_perm.tail             = 0;                          % Two-tailed test
            cfg_perm.clustertail      = 0;                          % Two-tailed for clustering
            cfg_perm.alpha            = 0.05;                       % Alpha level for the permutation test
            cfg_perm.correcttail      = 'alpha';                    % FDR correction for two-tailed T-test
            cfg_perm.numrandomization = 1000;                       % Number of randomizations
            
            
            nTaskTrials = size(data_task, 1);
            nRestTrials = size(data_rest, 1);
            
            cfg_perm.design = [ones(1,nTaskTrials), ones(1,nRestTrials)*2]; % design matrix;
            cfg_perm.ivar = 1; % Independent variable (Task vs. Rest)
            
            [stat_lfp] = ft_timelockstatistics(cfg_perm, LFP_task_erp, LFP_rest_erp);
            
            triggered_site_data.condition(1).event(e).stats.(Fin).val = stat_lfp;
            triggered_site_data.condition(2).event(e).stats.(Fin).val = stat_lfp;
            
            triggered_site_data.condition(1).event(e).stats.(Fin).cfg_perm = cfg_perm;
            triggered_site_data.condition(2).event(e).stats.(Fin).cfg_perm = cfg_perm;
            
%             pos_cluster_pvals = [stat_lfp.posclusters(:).prob];
%             neg_cluster_pvals = [stat_lfp.negclusters(:).prob];
            
        end
    end
end
out_triggered_site_data = triggered_site_data;
end