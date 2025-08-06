function [ triggered_site_data ] = ecg_bna_compute_triggered_variables_ECG( site , triggers, trials,cfg)

% suppress warning for xticklabel
warning ('off', 'MATLAB:hg:willberemoved');

switch cfg.to_trigger
    case 'LFP'
        FN ={'pow','powbp','lfp','pha','itpc','itpcbp'};
        samples_FN='LFP_samples';
    case 'MUA'
        FN ={'mua'};
        samples_FN='MUA_samples';
    case 'ECG'
        FN ={'ecg'};
        samples_FN='ECG_samples';
end

% folder to save sitewise results
site_results_folder = fullfile(cfg.sites_fldr);
if ~exist(site_results_folder, 'dir')
    mkdir(site_results_folder);
end
resampling_factor=site.tfs.resampling_factor;

u_blocks=unique(site.block);
trial_end_samples=0;
for b=u_blocks
    blocktrials=site.block==b;
    prev_block_end=trial_end_samples(end);
    trial_end_samples=[trial_end_samples prev_block_end+round(cumsum(site.(samples_FN)(blocktrials))*resampling_factor)];
    trial_end_samples(end)=prev_block_end+floor(sum(site.(samples_FN)(blocktrials))*resampling_factor);
end
trial_end_samples(1)=[];
trial_start_samples=[0 trial_end_samples(1:end-1)]+1;

%% convert to ipsi/contra
sitetrials=ph_get_unit_trials(site,trials);
positions=[sitetrials.tar_pos]-[sitetrials.fix_pos];
hemifields=num2cell(sign(real(positions)));
fixations=num2cell([sitetrials.fix_pos]);
positions=num2cell(positions);
[sitetrials.position]=deal(positions{:});
[sitetrials.hemifield]=deal(hemifields{:});
[sitetrials.fixation]=deal(fixations{:});
sitetrials=ph_LR_to_CI(cfg,site,sitetrials);  %% convert...

%% take over IDs
trig.condition = struct();
trig.site_ID = site.site_ID;
trig.session = site.session;
trig.target  = site.target;

% loop through conditions

for cn = 1:length(cfg.condition)
    % store details of analysed condition
    trig.condition(cn).label    = cfg.condition(cn).name;
    
    % get trial indices for the given condition
    cond_trials = ecg_bna_get_condition_trials(sitetrials, cfg.condition(cn));
    
    if sum(cond_trials) == 0
        continue;
    end
    % fprintf('Condition Trials = %i , Num Trial_idx = %i, CIX = %i\n',numel(cond_trials),sum(trial_idx), sum(cix))
    
    % loop through time windows around the states to analyse
    for e = 1:size(cfg.analyse_states, 1)
        state=cfg.analyse_states(e,:);
        event_name=state{1};
        width_in_samples=[floor([state{4}]*site.tfs.sr) ceil([state{5}]*site.tfs.sr)]; %% time to samples!
        
        con_end_samples=trial_end_samples(cond_trials);
        con_start_samples=trial_start_samples(cond_trials);
        time=[width_in_samples(1):width_in_samples(2)]/site.tfs.sr;
        
        trig.condition(cn).event(e).trials = find(cond_trials); %find(cond_trials);
        trig.condition(cn).event(e).ntrials = sum(cond_trials);
        trig.condition(cn).event(e).event_name = event_name;
        trig.condition(cn).event(e).time=time;
        trig.condition(cn).event(e).tfr_time=time;
        
        trig_con_s=triggers.([event_name '_real']);
        trig_con_s(trig_con_s<con_start_samples(1)-width_in_samples(1))=0;
        trig_con_s(trig_con_s>con_end_samples(end)-width_in_samples(2)-2)=0;
        
        inbetween=false(size(trig_con_s));
        for t=1:numel(con_end_samples)-1
            inbetween=inbetween | trig_con_s<con_start_samples(t+1) & trig_con_s>con_end_samples(t);
        end
        trig_con_s(inbetween)=0;
        
        realD = ecg_bna_get_triggered_parameters(site,trig_con_s, width_in_samples,cfg);
        
        %% compute shuffled power spectra, ITPC spectra, lfp, and bandpassed ITPC:
        trig_con_s=triggers.([event_name '_shuffled']);
        trig_con_s(trig_con_s<con_start_samples(1)-width_in_samples(1))=0;
        trig_con_s(trig_con_s>con_end_samples(end)-width_in_samples(2))=0;
        
        
        inbetween=false(size(trig_con_s));
        for t=1:numel(con_end_samples)-1
            inbetween=inbetween | trig_con_s<con_start_samples(t+1) & trig_con_s>con_end_samples(t);
        end
        trig_con_s(inbetween)=0;
        [shuffledD, significance]= ecg_bna_get_triggered_parameters(site,trig_con_s,width_in_samples,cfg,realD);
        
        normalized = ecg_bna_compute_shufflePredictor_normalization_general(realD,shuffledD,cfg);
        %significance = ecg_bna_compute_significance(realD,shuffledD,cfg);
        
        trig.condition(cn).event(e).real=realD;
        trig.condition(cn).event(e).shuffled=shuffledD;
        trig.condition(cn).event(e).normalized=normalized;
        trig.condition(cn).event(e).significance=significance;
        
    end
end
triggered_site_data = trig;
if isfield(cfg.lfp, 'TaskRest_SigClust') && cfg.lfp.TaskRest_SigClust == 1
    [ out_triggered_site_data] = ecg_bna_compute_cluster_sig_withinSites_between_cond( triggered_site_data, {'pha','pow','lfp'},cfg.analyse_states(:,1), cfg );
    
    triggered_site_data = out_triggered_site_data;
    clear out_triggered_site_data
end
% plots - if we don't shuffle, there will be no shuffled!
methods= {'real','shuffled','normalized'};
for mt = 1: numel(methods)
    switch cfg.to_trigger
        case 'LFP'
            ecg_bna_plots_per_site( trig, cfg, methods{mt}) % per site!
        case 'MUA'
            ecg_bna_plots_per_mua_site( trig, cfg, methods{mt}) % per site!
        case 'ECG'
            ecg_bna_plots_per_ecg_site( trig, cfg, methods{mt}) % per site!
    end
    
    % Note: ===> last input could be 'real', 'shuffled', or 'normalized'
end

if isfield(cfg.lfp, 'removeComplete') &&cfg.lfp.removeComplete==1
    for cn = 1:length(cfg.condition)
        for e = 1:size(cfg.analyse_states, 1)
            if ~isempty(triggered_site_data.condition(cn).event)
                for fin = 1:length(FN)
                    Fin = FN{fin};
                    triggered_site_data.condition(cn).event(e).real.(Fin)     = rmfield(triggered_site_data.condition(cn).event(e).real.(Fin), 'complete');
                    triggered_site_data.condition(cn).event(e).shuffled.(Fin) = rmfield(triggered_site_data.condition(cn).event(e).shuffled.(Fin), 'complete');
                end
            else
                continue;
            end
        end
    end
end


save(fullfile(site_results_folder, [trig.site_ID '.mat']), 'triggered_site_data');
% save(fullfile(site_results_folder, [trig.site_ID '.mat']), 'triggered_site_data','-v7.3');
close all;
end