function [ triggered_site_data ] = ecg_bna_compute_triggered_variables( site , triggers, trials,cfg)

% suppress warning for xticklabel
warning ('off', 'MATLAB:hg:willberemoved');

switch cfg.to_trigger
    case 'LFP';
        FN ={'pow','powbp','evoked','pha','itpc','itpcbp'};
    case 'MUA';
        FN ={'evoked'};
end
datatype=lower(cfg.to_trigger); % mua or lfp
site_results_folder=cfg.fldr.([cfg.to_trigger '_sites']);
samples_FN=[cfg.to_trigger '_samples'];

resampling_factor=site.(datatype).resampling_factor;

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

conditions_valid=true(size(cfg.condition));

% loop through events
for e = 1:size(cfg.analyse_states, 1)
    state=cfg.analyse_states(e,:);
    event_name=state{1};
    width_in_samples=[floor([state{4}]*site.(datatype).sr) ceil([state{5}]*site.(datatype).sr)]; %% time to samples!
    time=[width_in_samples(1):width_in_samples(2)]/site.(datatype).sr;
    
    % loop through conditions
    for cn = 1:length(cfg.condition)
        % store details of analysed condition
        trig.condition(cn).label    = cfg.condition(cn).name;
        
        % get trial indices for the given condition
        cond_trials = ecg_bna_get_condition_trials(sitetrials, cfg.condition(cn));
        if sum(cond_trials) == 0
            conditions_valid(cn)=false;
            continue;
        end
        % fprintf('Condition Trials = %i , Num Trial_idx = %i, CIX = %i\n',numel(cond_trials),sum(trial_idx), sum(cix))
        
        
        con_end_samples=trial_end_samples(cond_trials);
        con_start_samples=trial_start_samples(cond_trials);
        
        trig.condition(cn).event(e).trials = find(cond_trials); %find(cond_trials);
        trig.condition(cn).event(e).ntrials = sum(cond_trials);
        trig.condition(cn).event(e).event_name = event_name;
        trig.condition(cn).event(e).time=time;
        trig.condition(cn).event(e).tfr_time=time;
        
        trig_con_s=triggers.([event_name '_observed']);
        trig_con_s(trig_con_s<con_start_samples(1)-width_in_samples(1))=0;
        trig_con_s(trig_con_s>con_end_samples(end)-width_in_samples(2)-2)=0;
        
        inbetween=false(size(trig_con_s));
        for t=1:numel(con_end_samples)-1
            inbetween=inbetween | trig_con_s<con_start_samples(t+1) & trig_con_s>con_end_samples(t);
        end
        trig_con_s(inbetween)=0;
        
        observedD = ecg_bna_get_triggered_parameters(site,trig_con_s, width_in_samples,cfg);
        
        %% compute surrogate power spectra, ITPC spectra, lfp, and bandpassed ITPC:
        trig_con_s=triggers.([event_name '_surrogate']);
        trig_con_s(trig_con_s<con_start_samples(1)-width_in_samples(1))=0;
        trig_con_s(trig_con_s>con_end_samples(end)-width_in_samples(2))=0;
        
        
        inbetween=false(size(trig_con_s));
        for t=1:numel(con_end_samples)-1
            inbetween=inbetween | trig_con_s<con_start_samples(t+1) & trig_con_s>con_end_samples(t);
        end
        trig_con_s(inbetween)=0;
        tic
        [surrogateD, significance, cluster_z_values]= ecg_bna_get_triggered_parameters(site,trig_con_s,width_in_samples,cfg,observedD);
        toc
        normalized = ecg_bna_normalize(observedD,surrogateD,cfg);
        %significance = ecg_bna_compute_significance(observedD,surrogateD,cfg);
        
        trig.condition(cn).event(e).observed=observedD;
        trig.condition(cn).event(e).surrogate=surrogateD;
        trig.condition(cn).event(e).normalized=normalized;
        trig.condition(cn).event(e).significance=significance;
        trig.condition(cn).event(e).cluster_z_values=cluster_z_values;
        
%         if isfield(observedD,'evoked')
%             pre_samples=1:abs(width_in_samples(1))+1;
%             post_samples=max(pre_samples)+1:width_in_samples(2)+max(pre_samples);
%             
%             pre=mean(squeeze(observedD.evoked.complete(:,1,pre_samples)),2);
%             post=mean(squeeze(observedD.evoked.complete(:,1,post_samples)),2);
%             [h,p]=ttest(pre,post);
%             trig.condition(cn).event(e).prevspost.p=p;
%             trig.condition(cn).event(e).prevspost.h=h;
%         end
    end
    if ~isempty(cfg.lfp.compare_conditions)
        for ncomp=1:numel(cfg.lfp.compare_conditions)
            c1=cfg.lfp.compare_conditions{ncomp}(1);
            c2=cfg.lfp.compare_conditions{ncomp}(2);
            cond = [trig.condition(c1).label,'_',trig.condition(c2).label];
            if all (conditions_valid([c1,c2]))
                ev_temp = ecg_bna_compare_per_site(trig.condition(c1).event(e),trig.condition(c2).event(e),cfg);
                fntemp=fieldnames(ev_temp);
                for f=1:numel(fntemp)
                    fn=fntemp{f};
                    trig.(cond).event(e).(fn)=ev_temp.(fn);
                end
                
                normalized = ecg_bna_normalize(trig.(cond).event(e).observed,trig.(cond).event(e).surrogate,cfg);
                trig.(cond).event(e).normalized=normalized;
            else
                trig.(cond).event(e).ntriggers=0;
            end
            trig.(cond).event(e).event_name=cfg.analyse_states{e,1};
            trig.(cond).label=cond;
        end
    end
end

% plots - if we don't shuffle, there will be no surrogate!
if cfg.plot_per_site
    if ~isempty(cfg.lfp.compare_conditions) %% && both conditions exist...
        for ncomp=1:numel(cfg.lfp.compare_conditions)
            c1=cfg.lfp.compare_conditions{ncomp}(1);
            c2=cfg.lfp.compare_conditions{ncomp}(2);
            if all (conditions_valid([c1,c2]))
                switch cfg.to_trigger
                    case 'LFP';
                        ecg_bna_plots_per_site( trig, cond, cfg, 'normalized') % per site!
                    case 'MUA';
                        ecg_bna_plots_per_mua_site( trig, cond, cfg, 'normalized') % per site!
                end
            end
        end
    end
    methods= {'observed','surrogate','normalized'};
    for mt = 1: numel(methods)
        switch cfg.to_trigger
            case 'LFP';
                ecg_bna_plots_per_site( trig, 'condition',cfg, methods{mt}) % per site!
            case 'MUA';
                ecg_bna_plots_per_mua_site( trig,'condition', cfg, methods{mt}) % per site!
        end
        % Note: ===> last input could be 'observed', 'surrogate', or 'normalized'
    end
end
if isfield(cfg.lfp, 'removeComplete') &&cfg.lfp.removeComplete==1
    for cn = 1:length(cfg.condition)
        for e = 1:size(cfg.analyse_states, 1)
            if ~isempty(trig.condition(cn).event)
                for fin = 1:length(FN)
                    Fin = FN{fin};
                    trig.condition(cn).event(e).observed.(Fin)     = rmfield(trig.condition(cn).event(e).observed.(Fin), 'complete');
                    trig.condition(cn).event(e).surrogate.(Fin) = rmfield(trig.condition(cn).event(e).surrogate.(Fin), 'complete');
                end
            else
                continue;
            end
        end
    end
end

triggered_site_data = trig;
save(fullfile(site_results_folder, [trig.site_ID '.mat']), 'triggered_site_data');
close all;
end