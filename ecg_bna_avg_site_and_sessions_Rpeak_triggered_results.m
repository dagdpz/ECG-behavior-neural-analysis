function out = ecg_bna_avg_site_and_sessions_Rpeak_triggered_results(session_raw, ecg_bna_cfg,across)
% Condition-based Rpeak Triggered response average across many sites from
% a single session or multiple sessions
%
% USAGE:
%	session_Rpeak_triggered_raw = ecg_bna_avg_site_and_sessions_Rpeak_triggered_results(session_Rpeak_triggered_raw, ecg_bna_cfg,across)
%
% INPUTS:
%		session_Rpeak_triggered_raw	- struct containing the condition-based
%		Rpeak Triggered response for indiviual sites
%
%       Required Fields:
%               session.sites - 1xM struct containing condition-based
%               Rpeak Triggered responses for M sites
%
%       ecg_bna_cfg         - struct containing the required settings
%       Required Fields:
%               conditions          - trial conditions to compare, see
%               lfp_tfa_settings.m and lfp_tfa_compare_conditions.m
%               root_results_fldr   - root folder where results are saved
%               compare.targets     - targets to compare, see lfp_tfa_settings.m
%               ref_hemisphere      - reference hemisphere for ipsi and
%               contra labeling
%               diff_condition      - conditions to compare, the plot
%               for compared conditions would be shown one on top of the
%               other
%       Optional Fields:
%               diff_color          - color to be used for plotting the
%               compared conditions
%               diff_legend         - legend to be used while plotting the
%               compared conditions
% OUTPUTS:
%		sites_avg           - structure containing condition-based
%		Rpeak Triggered response averaged across multiple sites
%
%       futhur addition:
%       session_avg           - structure containing condition-based
%		Rpeak Triggered response averaged across sessions
%
% REQUIRES:	lfp_tfa_plot_evoked_lfp
%
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%


% results folder
results_fldr = fullfile(ecg_bna_cfg.analyse_lfp_folder, 'ECG_triggered_avg');
if ~exist(results_fldr, 'dir')
    mkdir(results_fldr);
end

avg = struct();
targets = {session_raw.sites(:).target};
targets = unique(targets);
for t = 1:length(targets)
    
    avg.(across)(t).target = targets{t};
    avg.(across)(t).session = session_raw.session;
    avg.(across)(t).site_ID = [session_raw.session,'_Population_of_',num2str(length(session_raw.sites)),'_Sites'];
    
    if ~strcmp(targets{t}, ecg_bna_cfg.compare.targets)
        continue;
    end
    
    for cn = 1:length(ecg_bna_cfg.conditions)
        
        fprintf('averaging Condition %s across %s\n', ecg_bna_cfg.conditions(cn).label, across);
        % initializing
        concat.real.pow = [];
        concat.real.itpc = [];
        concat.real.itpcbp = [];
        concat.real.powbp = [];
        concat.real.lfp = [];
        concat.real.nRpeaks = [];
        concat.ntrials = [];
        concat.normalized.pow = [];
        concat.normalized.itpc = [];
        concat.normalized.itpcbp = [];
        concat.normalized.powbp = [];
        concat.normalized.lfp = [];
        
        
        avg.(across)(t).condition(cn).label = ecg_bna_cfg.conditions(cn).label;
        avg.(across)(t).condition(cn).cfg_condition = ecg_bna_cfg.conditions(cn);       
        avg.(across)(t).condition(cn).state_hs = struct();
        
        if (session_raw.sites(1).condition(cn).ntrials ==0)
            avg.(across)(t).condition(cn).ntrials = 0;
            continue;
        end
        avg.(across)(t).condition(cn).state_hs.state = session_raw.sites(1).condition(cn).state_hs.state;
        avg.(across)(t).condition(cn).state_hs.state_name = session_raw.sites(1).condition(cn).state_hs.state_name;
        avg.(across)(t).condition(cn).state_hs.hs_label = session_raw.sites(1).condition(cn).state_hs.hs_label;
        
        avg.(across)(t).condition(cn).state_hs.real = struct();
        avg.(across)(t).condition(cn).state_hs.shuffled = struct();
        avg.(across)(t).condition(cn).state_hs.normalized = struct();
        avg.(across)(t).condition(cn).state_hs.significance = struct();
        % real 
        avg.(across)(t).condition(cn).state_hs.real.time = session_raw.sites(1).condition(cn).state_hs.real.time;
        avg.(across)(t).condition(cn).state_hs.real.tfr_time = session_raw.sites(1).condition(cn).state_hs.real.tfr_time;
        avg.(across)(t).condition(cn).state_hs.real.freq = session_raw.sites(1).condition(cn).state_hs.real.freq;
        avg.(across)(t).condition(cn).state_hs.real.state = session_raw.sites(1).condition(cn).state_hs.state;
        avg.(across)(t).condition(cn).state_hs.real.state_name = session_raw.sites(1).condition(cn).state_hs.state_name;
        
        % normalized
        avg.(across)(t).condition(cn).state_hs.normalized.time = session_raw.sites(1).condition(cn).state_hs.normalized.time;
        avg.(across)(t).condition(cn).state_hs.normalized.tfr_time = session_raw.sites(1).condition(cn).state_hs.normalized.tfr_time;
        avg.(across)(t).condition(cn).state_hs.normalized.freq = session_raw.sites(1).condition(cn).state_hs.normalized.freq;
        avg.(across)(t).condition(cn).state_hs.normalized.state = session_raw.sites(1).condition(cn).state_hs.state;
        avg.(across)(t).condition(cn).state_hs.normalized.state_name = session_raw.sites(1).condition(cn).state_hs.state_name;

        % having significance and shuffle in place for plotting function to work 
        % ==> it's now only rand 0,1
        avg.(across)(t).condition(cn).state_hs.significance.pow = randi([0 1],size(session_raw.sites(1).condition(cn).state_hs.significance.pow));
        avg.(across)(t).condition(cn).state_hs.significance.itpc = randi([0 1],size(session_raw.sites(1).condition(cn).state_hs.significance.itpc));
        avg.(across)(t).condition(cn).state_hs.significance.itpcbp = randi([0 1],size(session_raw.sites(1).condition(cn).state_hs.significance.itpcbp));
        avg.(across)(t).condition(cn).state_hs.significance.powbp = randi([0 1],size(session_raw.sites(1).condition(cn).state_hs.significance.powbp));
        avg.(across)(t).condition(cn).state_hs.significance.lfp = randi([0 1],size(session_raw.sites(1).condition(cn).state_hs.significance.lfp));
        % shuffled
        avg.(across)(t).condition(cn).state_hs.shuffled.pow = randi([0 1],size(session_raw.sites(1).condition(cn).state_hs.shuffled.pow));
        avg.(across)(t).condition(cn).state_hs.shuffled.itpc = randi([0 1],size(session_raw.sites(1).condition(cn).state_hs.shuffled.itpc));
        avg.(across)(t).condition(cn).state_hs.shuffled.itpcbp = randi([0 1],size(session_raw.sites(1).condition(cn).state_hs.shuffled.itpcbp));
        avg.(across)(t).condition(cn).state_hs.shuffled.powbp = randi([0 1],size(session_raw.sites(1).condition(cn).state_hs.shuffled.powbp));
        avg.(across)(t).condition(cn).state_hs.shuffled.lfp = randi([0 1],size(session_raw.sites(1).condition(cn).state_hs.shuffled.lfp));
        avg.(across)(t).condition(cn).state_hs.shuffled.time = session_raw.sites(1).condition(cn).state_hs.real.time;
        avg.(across)(t).condition(cn).state_hs.shuffled.tfr_time = session_raw.sites(1).condition(cn).state_hs.real.tfr_time;
        avg.(across)(t).condition(cn).state_hs.shuffled.freq = session_raw.sites(1).condition(cn).state_hs.real.freq;
        
        

        for i = 1:length(session_raw.sites)
            if ~strcmp(targets{t},session_raw.sites(i).target)
                continue;
            end
            
            % concatnating the normalized variables of each sites:
            concat.real.nRpeaks = [session_raw.sites(i).condition(cn).state_hs.real.nRpeaks];
            concat.ntrials = [session_raw.sites(i).condition(cn).state_hs.ntrials];
            % real
            concat.real.pow = cat(4,concat.real.pow,session_raw.sites(i).condition(cn).state_hs.real.pow.mean);
            concat.real.itpc = cat(4,concat.real.itpc,session_raw.sites(i).condition(cn).state_hs.real.itpc.mean);
            concat.real.powbp = cat(4,concat.real.powbp,session_raw.sites(i).condition(cn).state_hs.real.powbp.mean);
            concat.real.itpcbp = cat(4,concat.real.itpcbp, session_raw.sites(i).condition(cn).state_hs.real.itpcbp.mean);
            concat.real.lfp = cat(3,concat.real.lfp, session_raw.sites(i).condition(cn).state_hs.real.lfp.mean);
            % normalized
            concat.normalized.pow = cat(4,concat.normalized.pow,session_raw.sites(i).condition(cn).state_hs.normalized.pow.mean);
            concat.normalized.itpc = cat(4,concat.normalized.itpc,session_raw.sites(i).condition(cn).state_hs.normalized.itpc.mean);
            concat.normalized.powbp = cat(4,concat.normalized.powbp,session_raw.sites(i).condition(cn).state_hs.normalized.powbp.mean);
            concat.normalized.itpcbp = cat(4,concat.normalized.itpcbp, session_raw.sites(i).condition(cn).state_hs.normalized.itpcbp.mean);
            concat.normalized.lfp = cat(3,concat.normalized.lfp, session_raw.sites(i).condition(cn).state_hs.normalized.lfp.mean);


        end
        
        % number of total trials and Rpeaks
        avg.(across)(t).condition(cn).ntrials = sum(concat.ntrials);
        avg.(across)(t).condition(cn).state_hs.real.nRpeaks = sum(concat.real.nRpeaks);

        
        % mean
        avg.(across)(t).condition(cn).state_hs.real.pow.mean = nanmean(concat.real.pow,4);
        avg.(across)(t).condition(cn).state_hs.real.itpc.mean = nanmean(concat.real.itpc,4);
        avg.(across)(t).condition(cn).state_hs.real.powbp.mean = nanmean(concat.real.powbp,4);
        avg.(across)(t).condition(cn).state_hs.real.itpcbp.mean = nanmean(concat.real.itpcbp,4);
        avg.(across)(t).condition(cn).state_hs.real.lfp.mean = nanmean(concat.real.lfp,3);
        % %
        avg.(across)(t).condition(cn).state_hs.normalized.pow.mean = nanmean(concat.normalized.pow,4);
        avg.(across)(t).condition(cn).state_hs.normalized.itpc.mean = nanmean(concat.normalized.itpc,4);
        avg.(across)(t).condition(cn).state_hs.normalized.powbp.mean = nanmean(concat.normalized.powbp,4);
        avg.(across)(t).condition(cn).state_hs.normalized.itpcbp.mean = nanmean(concat.normalized.itpcbp,4);
        avg.(across)(t).condition(cn).state_hs.normalized.lfp.mean = nanmean(concat.normalized.lfp,3);
        % std
        avg.(across)(t).condition(cn).state_hs.real.pow.std = nanstd(concat.real.pow,0,4);
        avg.(across)(t).condition(cn).state_hs.real.itpc.std = nanstd(concat.real.itpc,0,4);
        avg.(across)(t).condition(cn).state_hs.real.powbp.std = nanstd(concat.real.powbp,0,4);
        avg.(across)(t).condition(cn).state_hs.real.itpcbp.std = nanstd(concat.real.itpcbp,0,4);
        avg.(across)(t).condition(cn).state_hs.real.lfp.std = nanstd(concat.real.lfp,0,3);
        % %
        avg.(across)(t).condition(cn).state_hs.normalized.pow.std = nanstd(concat.normalized.pow,0,4);
        avg.(across)(t).condition(cn).state_hs.normalized.itpc.std = nanstd(concat.normalized.itpc,0,4);
        avg.(across)(t).condition(cn).state_hs.normalized.powbp.std = nanstd(concat.normalized.powbp,0,4);
        avg.(across)(t).condition(cn).state_hs.normalized.itpcbp.std = nanstd(concat.normalized.itpcbp,0,4);
        avg.(across)(t).condition(cn).state_hs.normalized.lfp.std = nanstd(concat.normalized.lfp,0,3);
               
    end
    % condition to plot:
    cond_cfg = [avg.(across)(t).condition(:).cfg_condition];
    % plotting
    ecg_bna_cfg.sites_lfp_fldr   = results_fldr;
    data2plot = avg.(across)(t);
    ecg_bna_plots_per_session( data2plot, cond_cfg, ecg_bna_cfg, 'normalized')
    
end

out = avg.(across);

end