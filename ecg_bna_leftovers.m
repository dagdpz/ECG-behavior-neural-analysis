
% Baseline power calculation - this needs to move to somewhere entirely different (?)
% site_lfp = lfp_tfa_compute_site_baseline( site_lfp, session_info, lfp_tfa_cfg );
%     allsites_lfp = [allsites_lfp, site_lfp];
% end

% %% calculate cross power spectrum between sites and sync measure spectrogram - PROBABLY NOT WORKING at the moment
% if any(strcmp(lfp_tfa_cfg.analyses, 'sync')) || any(strcmp(lfp_tfa_cfg.analyses, 'syncsp'))
%     % prepare results folder
%     results_fldr = fullfile(session_info.proc_results_fldr, 'crossspectrum');
%     if ~exist(results_fldr, 'dir')
%         mkdir(results_fldr);
%     end
%     % loop through each site
%     for i = 1:length(allsites_lfp)-1
%         site1_lfp = allsites_lfp(i);
%         % pair a site
%         for j = i+1:length(allsites_lfp)
%             site2_lfp = allsites_lfp(j);
%             fprintf('Computing cross power spectrum for site pair %s - %s\n', ...
%                 site1_lfp.site_ID, site2_lfp.site_ID);
%             sitepair_crosspow = lfp_tfa_compute_sitepair_csd(site1_lfp, site2_lfp, lfp_tfa_cfg);
%             % save data
%             results_mat = fullfile(results_fldr, ['sites_crosspow_', sitepair_crosspow.sites{1} '-' sitepair_crosspow.sites{2} '.mat']);
%             save(results_mat, 'sitepair_crosspow', '-v7.3');
%
%         end
%     end
% end

%% calculate sync measure spectrum
%     if any(strcmp(lfp_tfa_cfg.analyses, 'syncspctrm'))
%         % loop through each site
%         for i = 1:length(allsites_lfp)-1
%             site1_lfp = allsites_lfp(i);
%             % pair a site
%             for j = i+1:length(allsites_lfp)
%                 site2_lfp = allsites_lfp(j);
%                 fprintf('Computing sync spectrum for site pair %s - %s\n', ...
%                     site1_lfp.site_ID, site2_lfp.site_ID);
%                 % compute ppc spectrum between sitepair
%                 % get the trial conditions for this session
%                 conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg, ...
%                     {session_info.Preinj_blocks, session_info.Postinj_blocks});
%                 sitepair_syncspctrm = lfp_tfa_sitepair_avg_syncspctrum(site1_lfp, site2_lfp, conditions, lfp_tfa_cfg);
%             end
%         end
%     end






%%
% 
% % Average across sites for a session
% session_avg = struct();
% % targets for this session
% targets = unique({session_proc_lfp.target});
% % average each target separately
% for t = 1:length(targets)
%     session_avg(t).target = targets{t};
%     
%     % conditions
%     for cn = 1:length(site_conditions)
%         session_avg(t).condition(cn).state_hs = struct();
%         session_avg(t).condition(cn).cfg_condition = site_conditions(cn);
%         session_avg(t).condition(cn).label = site_conditions(cn).label;
%         session_avg(t).condition(cn).condition = site_conditions(cn);
%         session_avg(t).condition(cn).label = site_conditions(cn).label;
%         session_avg(t).condition(cn).session = session_proc_lfp(1).session;
%         session_avg(t).condition(cn).target = session_proc_lfp(1).target;
%         
%         for st = 1:size(analyse_states, 1)
%             for hs = 1:length(hs_labels)
%                 session_avg(t).condition(cn).state_hs(st, hs).nsites = 0;
%                 session_avg(t).condition(cn).state_hs(st, hs).lfp = [];
%                 session_avg(t).condition(cn).state_hs(st, hs).mean = [];
%                 session_avg(t).condition(cn).state_hs(st, hs).std = [];
%                 session_avg(t).condition(cn).state_hs(st, hs).time = [];
%                 session_avg(t).condition(cn).state_hs(st, hs).shuffled_lfp = [];
%             end
%         end
%         for i = 1:length(sites_data)
%             % if the site's target is same as target being considered
%             if ~strcmp(session_proc_lfp(i).target, targets{t})
%                 continue;
%             end
%             if sites_data(i).use_for_avg
%                 % calculate the average evoked LFP across sites for this condition
%                 if isempty(sites_data(i).condition(cn).state_hs) || ~isfield(sites_data(i).condition(cn).state_hs, 'mean')
%                     continue;
%                 end
%                 for hs = 1:size(sites_data(i).condition(cn).state_hs, 2)
%                     for st = 1:size(sites_data(i).condition(cn).state_hs, 1)
%                         if isempty(sites_data(i).condition(cn).state_hs(st, hs).mean)
%                             continue;
%                         end
%                         session_avg(t).condition(cn).state_hs(st, hs).nsites = session_avg(t).condition(cn).state_hs(st, hs).nsites + 1;
%                         if session_avg(t).condition(cn).state_hs(st, hs).nsites == 1
%                             % subtracting shuffled mean
%                             session_avg(t).condition(cn).state_hs(st, hs).shuffled_lfp = [session_avg(t).condition(cn).state_hs(st, hs).shuffled_lfp;...
%                                 sites_data(i).condition(cn).state_hs(st, hs).shuffled] ;
%                             session_avg(t).condition(cn).state_hs(st, hs).lfp = [session_avg(t).condition(cn).state_hs(st, hs).lfp;...
%                                 sites_data(i).condition(cn).state_hs(st, hs).mean - sites_data(i).condition(cn).state_hs(st, hs).shuffled] ;
%                             session_avg(t).condition(cn).state_hs(st, hs).time = sites_data(i).condition(cn).state_hs(st, hs).time;
%                             
%                         else % nsamples differs (hopefully only slightly (?) across sites
%                             nsamples = length(session_avg(t).condition(cn).state_hs(st, hs).time);
%                             if nsamples > length(sites_data(i).condition(cn).state_hs(st, hs).time)
%                                 nsamples = length(sites_data(i).condition(cn).state_hs(st, hs).time);
%                             end
%                             % subtracting shuffled mean
%                             session_avg(t).condition(cn).state_hs(st, hs).shuffled_lfp  = [session_avg(t).condition(cn).state_hs(st, hs).shuffled_lfp(:,1:nsamples); ...
%                                 sites_data(i).condition(cn).state_hs(st, hs).shuffled(1:nsamples)] ;
%                             session_avg(t).condition(cn).state_hs(st, hs).lfp  = [session_avg(t).condition(cn).state_hs(st, hs).lfp(:,1:nsamples); ...
%                                 sites_data(i).condition(cn).state_hs(st, hs).mean(1:nsamples)- sites_data(i).condition(cn).state_hs(st, hs).shuffled(1:nsamples)] ;
%                             session_avg(t).condition(cn).state_hs(st, hs).time = session_avg(t).condition(cn).state_hs(st, hs).time(1:nsamples) ;
%                             
%                         end
%                         % struct to store average evoked LFP across sites
%                         session_avg(t).condition(cn).state_hs(st, hs).hs_label = sites_data(i).condition(cn).state_hs(st, hs).hs_label;
%                         if isfield(sites_data(i).condition(cn).state_hs(st, hs), 'state') && isfield(sites_data(i).condition(cn).state_hs(st, hs), 'state_name')
%                             session_avg(t).condition(cn).state_hs(st, hs).state = sites_data(i).condition(cn).state_hs(st, hs).state;
%                             session_avg(t).condition(cn).state_hs(st, hs).state_name = sites_data(i).condition(cn).state_hs(st, hs).state_name;
%                         end
%                     end
%                     
%                 end
%             end
%         end
%         % average LFP across sites for a session
%         for hs = 1:size(session_avg(t).condition(cn).state_hs, 2)
%             for st = 1:size(session_avg(t).condition(cn).state_hs, 1)
%                 session_avg(t).condition(cn).state_hs(st, hs).mean = nanmean(session_avg(t).condition(cn).state_hs(st, hs).lfp, 1);
%                 session_avg(t).condition(cn).state_hs(st, hs).std = nanstd(session_avg(t).condition(cn).state_hs(st, hs).lfp, 0, 1);
%                 session_avg(t).condition(cn).state_hs(st, hs).shuffled = nanmean(session_avg(t).condition(cn).state_hs(st, hs).shuffled_lfp, 1);
%                 session_avg(t).condition(cn).state_hs(st, hs).shuffled = nanstd(session_avg(t).condition(cn).state_hs(st, hs).shuffled_lfp, 0, 1);
%             end
%         end
%         % plot average evoked LFP across sites for this session
%         plottitle = ['Session: ', session_avg(t).condition(cn).session ', Target = ' session_avg(t).condition(cn).target ' (ref_' ecg_bna_cfg.ref_hemisphere '), ' site_conditions(cn).label ', '];
%         if site_conditions(cn).choice == 0
%             plottitle = [plottitle 'Instructed trials'];
%         elseif site_conditions(cn).choice == 1
%             plottitle = [plottitle 'Choice trials'];
%         end
%         result_file = fullfile(results_folder_evoked, ['LFP_Evoked_' session_avg(t).condition(cn).session '_' session_avg(t).condition(cn).target '_' site_conditions(cn).label]);
%         ecg_bna_plot_evoked_lfp (session_avg(t).condition(cn).state_hs, ecg_bna_cfg, plottitle, result_file);
%     end
%     
%     
%     % difference between conditions
%     session_avg(t).difference = [];
%     for diff = 1:size(ecg_bna_cfg.diff_condition, 2)
%         diff_condition = ecg_bna_cfg.diff_condition{diff};
%         diff_color = []; diff_legend = [];
%         if isfield(ecg_bna_cfg, 'diff_color')
%             diff_color = ecg_bna_cfg.diff_color{diff};
%         end
%         if isfield(ecg_bna_cfg, 'diff_legend')
%             diff_legend = ecg_bna_cfg.diff_legend{diff};
%         end
%         session_avg(t).difference = [session_avg(t).difference, ...
%             ecg_bna_compute_diff_condition_average('Rpeak_evoked_LFP', session_avg(t).condition, diff_condition, diff_color, diff_legend)];
%     end
%     % plot Difference TFR
%     for dcn = 1:length(session_avg(t).difference)
%         if ~isempty(session_avg(t).difference(dcn).state_hs)
%             if isfield(session_avg(t).difference(dcn).state_hs,'mean')
%                 plottitle = ['Target ', ecg_bna_cfg.compare.targets{t}, ' (ref_', ecg_bna_cfg.ref_hemisphere, ') ', session_avg(t).difference(dcn).label];
%                 result_file = fullfile(results_folder_evoked, ['LFP_DiffEvoked_' ecg_bna_cfg.compare.targets{t} '_' 'diff_condition' num2str(dcn)]);
%                 %session_avg(t).difference(dcn).label '.png']);
%                 ecg_bna_plot_evoked_lfp(session_avg(t).difference(dcn).state_hs, ecg_bna_cfg, plottitle, result_file);
%             end
%         end
%     end
%     
% end
% 
% close all;
% 
% % store session average data
% session_data.session_avg = session_avg;
% 
% save mat files
% save(fullfile(session_result_folder, ['Rpeak_triggered_session_' session_data.session '.mat']), 'session_data');

