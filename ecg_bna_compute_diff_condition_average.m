function [ diff_cond_avg ] = ecg_bna_compute_diff_condition_average( analysis, cond_avg, diff_condition, diff_color, diff_legend )
%ecg_bna_compute_diff_condition_average - function to compare the ECG 
%related neural and behavior analysis averages for given conditions
%
% USAGE:
%	[ diff_evoked ] = ecg_bna_compute_diff_condition_average( analysis,
%	cond_avg, diff_condition, diff_color, diff_legend ) 
%
% INPUTS:
%       analysis        - the kind of analysis whose results are to be
%       averaged
%       cond_avg        - 1xnC struct containing the average analysis results
%       for nC conditions
%       diff_condition  - conditions to be compared
%       diff_color      - colors to be used for plotting the different
%       conditions being compared
%       diff_legend     - legend to be used for conditions being compared
% 
% OUTPUTS:
%		diff_cond_avg    - structure containing results of comparison of
%		given conditions 
%
%


    conditions = [cond_avg.cfg_condition];
    
    cfg_condition1 = struct();
    cfg_condition2 = struct();
    
    for i = 1:length(diff_condition)/2
        diff_cond_avg = struct();
        compare.field = diff_condition{1, 2*(i-1) + 1};
        compare.values = diff_condition{1, 2*(i-1) + 2};
        % check if conditions to compare exist
        if strcmp(compare.field, 'perturbation')
            if sum([compare.values{:}] == unique([conditions.perturbation])) < 2
                continue;
            end
        elseif strcmp(compare.field, 'choice')
            if sum([compare.values{:}] == unique([conditions.choice])) < 2
                continue;
            end
        elseif strcmp(compare.field, 'success')
            if sum([compare.values{:}] == unique([conditions.success])) < 2
                continue;
            end
        elseif strcmp(compare.field, 'type_eff')
            if sum(ismember(vertcat(compare.values{:}), ...
                unique([conditions.type; conditions.effector]', 'rows'), 'rows')) < 2
                continue;
            end
        end
    
        dcn = 0;
        traversed_idx = [];
        for cn = 1:length(cond_avg)
            condition_found = false;
            if strcmp(compare.field, 'choice')
                condition_found = cond_avg(cn).cfg_condition.choice == compare.values{1};

            elseif strcmp(compare.field, 'perturbation')
                condition_found = cond_avg(cn).cfg_condition.perturbation == compare.values{1};
            
            elseif strcmp(compare.field, 'success')
                condition_found = cond_avg(cn).cfg_condition.success == compare.values{1};
                
            elseif strcmp(compare.field, 'type_eff')
                condition_found = cond_avg(cn).cfg_condition.type == compare.values{1}(1) ...
                    & cond_avg(cn).cfg_condition.effector == compare.values{1}(2);

            end
            % initially load the pre-injection data structure
            if condition_found
                traversed_idx = [traversed_idx cn];            
            else
                continue;
            end
            for d = 1:length(cond_avg)
                if any(traversed_idx == d), continue; end
                comparison_pair_found = false;

                if strcmp(compare.field, 'choice')
                    comparison_pair_found = cond_avg(d).cfg_condition.type == cond_avg(cn).cfg_condition.type ...
                        & cond_avg(d).cfg_condition.effector == cond_avg(cn).cfg_condition.effector ...
                        & cond_avg(d).cfg_condition.choice == compare.values{2} ...
                        & cond_avg(d).cfg_condition.perturbation == cond_avg(cn).cfg_condition.perturbation ...
                        & cond_avg(d).cfg_condition.success == cond_avg(cn).cfg_condition.success;

                elseif strcmp(compare.field, 'perturbation')
                    comparison_pair_found = cond_avg(d).cfg_condition.type == cond_avg(cn).cfg_condition.type ...
                        & cond_avg(d).cfg_condition.effector == cond_avg(cn).cfg_condition.effector ...
                        & cond_avg(d).cfg_condition.choice == cond_avg(cn).cfg_condition.choice ...
                        & cond_avg(d).cfg_condition.perturbation == compare.values{2} ...
                        & cond_avg(d).cfg_condition.success == cond_avg(cn).cfg_condition.success;
                    
                elseif strcmp(compare.field, 'success')
                    comparison_pair_found = cond_avg(d).cfg_condition.type == cond_avg(cn).cfg_condition.type ...
                        & cond_avg(d).cfg_condition.effector == cond_avg(cn).cfg_condition.effector ...
                        & cond_avg(d).cfg_condition.choice == cond_avg(cn).cfg_condition.choice ...
                        & cond_avg(d).cfg_condition.success == compare.values{2} ...
                        & cond_avg(d).cfg_condition.perturbation == cond_avg(cn).cfg_condition.perturbation;
                    
                elseif strcmp(compare.field, 'type_eff')
                    comparison_pair_found = cond_avg(d).cfg_condition.type == compare.values{2}(1) ...
                        & cond_avg(d).cfg_condition.effector == compare.values{2}(2) ...
                        & cond_avg(d).cfg_condition.choice == cond_avg(cn).cfg_condition.choice ...
                        & cond_avg(d).cfg_condition.success == cond_avg(cn).cfg_condition.success ...
                        & cond_avg(d).cfg_condition.perturbation == cond_avg(cn).cfg_condition.perturbation;
                end
                if comparison_pair_found

                    dcn = dcn + 1;
                    % pre injection
                    cond1_avg = cond_avg(cn);
                    % post injection
                    cond2_avg = cond_avg(d);
                                       
                    diff_cond_avg.difference(dcn) = cond2_avg;
                    
                                      
                    diff_cond_avg.difference(dcn).cfg_condition = cond2_avg.cfg_condition;
                    
                    plot_legend = cell(1,2);
%                     
                    if nargin > 3 && ~isempty(diff_legend)
                        plot_legend = diff_legend;
                    elseif strcmp(compare.field, 'choice')                        
                        cfg_condition1.choice = cond2_avg.cfg_condition.choice;
                        cfg_condition2.choice = cond1_avg.cfg_condition.choice;
                        plot_legend{1} = lfp_tfa_get_condition_label(cfg_condition1, 'long');
                        plot_legend{2} = lfp_tfa_get_condition_label(cfg_condition2, 'long');
                    elseif strcmp(compare.field, 'perturbation')
                        cfg_condition1.perturbation = cond2_avg.cfg_condition.perturbation;
                        cfg_condition2.perturbation = cond1_avg.cfg_condition.perturbation;
                        plot_legend{1} = lfp_tfa_get_condition_label(cfg_condition1, 'long');
                        plot_legend{2} = lfp_tfa_get_condition_label(cfg_condition2, 'long');
                    elseif strcmp(compare.field, 'success')
                        cfg_condition1.success = cond2_avg.cfg_condition.success;
                        cfg_condition2.success = cond1_avg.cfg_condition.success;
                        plot_legend{1} = lfp_tfa_get_condition_label(cfg_condition1, 'long');
                        plot_legend{2} = lfp_tfa_get_condition_label(cfg_condition2, 'long');
                    elseif strcmp(compare.field, 'type_eff')
                        diff_cond_avg.difference(dcn).cfg_condition.type_eff = nan;
                        cfg_condition1.type = cond2_avg.cfg_condition.type;
                        cfg_condition1.effector = cond2_avg.cfg_condition.effector;
                        cfg_condition2.type = cond1_avg.cfg_condition.type;
                        cfg_condition2.effector = cond1_avg.cfg_condition.effector;
                        plot_legend{1} = lfp_tfa_get_condition_label(cfg_condition1, 'long');
                        plot_legend{2} = lfp_tfa_get_condition_label(cfg_condition2, 'long');
                    end
                    
                    % change the condition label
                    switch (compare.field)
                        case ('choice')
                            diff_cond_avg.difference(dcn).cfg_condition.choice = -1;
                        case ('perturbation')
                            diff_cond_avg.difference(dcn).cfg_condition.perturbation = -1;
                        case ('success')
                            diff_cond_avg.difference(dcn).cfg_condition.success = -1;
                        case ('type_eff')
                            diff_cond_avg.difference(dcn).cfg_condition.type_eff = -1;
                    end                    
                    diff_cond_avg.difference(dcn).label = ecg_bna_get_condition_label(...
                        diff_cond_avg.difference(dcn).cfg_condition, 'long'); 
                    
                    if strcmp(analysis, 'Rpeak_evoked_LFP') || ...
                            strcmp(analysis, 'Rpeak_evoked_ECG') || ...
                            strcmp(analysis, 'Event_evoked_ECG_R2Rt')
                        
                        if ~isfield(cond1_avg.hs_tuned_evoked, 'mean')
                            continue;
                        end
                        if ~isfield(cond2_avg.hs_tuned_evoked, 'mean')
                            continue;
                        end                    
                        % loop through handspace tunings
                        diff_cond_avg.difference(dcn).hs_tuned_evoked = cond2_avg.hs_tuned_evoked;
                        for hs = 1:size(cond2_avg.hs_tuned_evoked, 2)
                            for st = 1:size(cond2_avg.hs_tuned_evoked, 1)

                                if isfield(cond1_avg.hs_tuned_evoked(st, hs), 'mean') && ...
                                        isfield(cond2_avg.hs_tuned_evoked(st, hs), 'mean') && ...
                                        ~isempty(cond1_avg.hs_tuned_evoked(st, hs).mean) && ...
                                    ~isempty(cond2_avg.hs_tuned_evoked(st, hs).mean)
                                    ntimebins = min([length(cond2_avg.hs_tuned_evoked(st, hs).time), ...
                                        length(cond1_avg.hs_tuned_evoked(st, hs).time)]);

                                    diff_cond_avg.difference(dcn).hs_tuned_evoked(st, hs).mean = ...
                                        [cond2_avg.hs_tuned_evoked(st, hs).mean(:,1:ntimebins); ...
                                        cond1_avg.hs_tuned_evoked(st, hs).mean(:,1:ntimebins)];
                                    diff_cond_avg.difference(dcn).hs_tuned_evoked(st, hs).std = ...
                                        [cond2_avg.hs_tuned_evoked(st, hs).std(:,1:ntimebins); ...
                                        cond1_avg.hs_tuned_evoked(st, hs).std(:,1:ntimebins)];
                                    diff_cond_avg.difference(dcn).hs_tuned_evoked(st, hs).time = ...
                                        cond2_avg.hs_tuned_evoked(st, hs).time(1:ntimebins);  
                                    if isfield(cond2_avg.hs_tuned_evoked(st, hs), 'ntrials')
                                        diff_cond_avg.difference(dcn).hs_tuned_evoked(st, hs).ntrials = ...
                                            [];
                                    elseif isfield(cond2_avg.hs_tuned_evoked(st, hs), 'npeaks')
                                        diff_cond_avg.difference(dcn).hs_tuned_evoked(st, hs).npeaks = ...
                                            [];
                                    end
                                    diff_cond_avg.difference(dcn).hs_tuned_evoked(st, hs).legend = ...
                                        plot_legend;

                                    if nargin > 2 && ~isempty(diff_color)
                                        diff_cond_avg.difference(dcn).hs_tuned_evoked(st, hs).color = ...
                                            diff_color;
                                    end

                                else
                                    %diff_evoked.difference(dcn).hs_tuned_evoked(st, hs).lfp = [];
                                    diff_cond_avg.difference(dcn).hs_tuned_evoked(st, hs).time = [];
                                    diff_cond_avg.difference(dcn).hs_tuned_evoked(st, hs).mean = [];
                                    diff_cond_avg.difference(dcn).hs_tuned_evoked(st, hs).std = [];
                                    diff_cond_avg.difference(dcn).hs_tuned_evoked(st, hs).legend = [];
                                    diff_cond_avg.difference(dcn).hs_tuned_evoked(st, hs).color = [];
                                end
                            end
                        end
                    elseif strcmp(analysis, 'Rpeak_evoked_states_onset')
                        if ~isfield(cond1_avg.Rpeak_evoked, 'abs_timeprob') ...
                                || ~isfield(cond1_avg.Rpeak_evoked, 'rel_timeprob')
                            continue;
                        end
                        if ~isfield(cond2_avg.Rpeak_evoked, 'abs_timeprob') ...
                                || ~isfield(cond2_avg.Rpeak_evoked, 'rel_timeprob')
                            continue;
                        end
                        % loop through handspace tunings
                        diff_cond_avg.difference(dcn).Rpeak_evoked = cond2_avg.Rpeak_evoked;
                        for hs = 1:size(cond2_avg.Rpeak_evoked, 2)
                            for st = 1:size(cond2_avg.Rpeak_evoked, 1)

                                if isfield(cond1_avg.Rpeak_evoked(st, hs), 'abs_timeprob') && ...
                                        isfield(cond2_avg.Rpeak_evoked(st, hs), 'abs_timeprob') && ...
                                        ~isempty(cond1_avg.Rpeak_evoked(st, hs).abs_timeprob) && ...
                                        ~isempty(cond2_avg.Rpeak_evoked(st, hs).abs_timeprob) && ...
                                        isfield(cond1_avg.Rpeak_evoked(st, hs), 'rel_timeprob') && ...
                                        isfield(cond2_avg.Rpeak_evoked(st, hs), 'rel_timeprob') && ...
                                        ~isempty(cond1_avg.Rpeak_evoked(st, hs).rel_timeprob) && ...
                                        ~isempty(cond2_avg.Rpeak_evoked(st, hs).rel_timeprob)
                                    diff_cond_avg.difference(dcn).Rpeak_evoked(st, hs).abs_timefromRpeak = ...
                                        {cond2_avg.Rpeak_evoked(st, hs).abs_timefromRpeak, ...
                                        cond1_avg.Rpeak_evoked(st, hs).abs_timefromRpeak};
                                    diff_cond_avg.difference(dcn).Rpeak_evoked(st, hs).abs_timeprob = ...
                                        [cond2_avg.Rpeak_evoked(st, hs).abs_timeprob, ...
                                        cond1_avg.Rpeak_evoked(st, hs).abs_timeprob];
                                    diff_cond_avg.difference(dcn).Rpeak_evoked(st, hs).rel_timefromRpeak = ...
                                        {cond2_avg.Rpeak_evoked(st, hs).rel_timefromRpeak, ...
                                        cond1_avg.Rpeak_evoked(st, hs).rel_timefromRpeak};
                                    diff_cond_avg.difference(dcn).Rpeak_evoked(st, hs).rel_timeprob = ...
                                        [cond2_avg.Rpeak_evoked(st, hs).rel_timeprob, ...
                                        cond1_avg.Rpeak_evoked(st, hs).rel_timeprob];
                                    if isfield(cond2_avg.Rpeak_evoked(st, hs), 'ntrials')
                                        diff_cond_avg.difference(dcn).Rpeak_evoked(st, hs).ntrials = ...
                                            [];
                                    end
                                    diff_cond_avg.difference(dcn).Rpeak_evoked(st, hs).legend = ...
                                        plot_legend;
                                    if nargin > 2 && ~isempty(diff_color)
                                        diff_cond_avg.difference(dcn).Rpeak_evoked(st, hs).colors = ...
                                            diff_color;
                                    end

                                else
                                    diff_cond_avg.difference(dcn).Rpeak_evoked(st, hs).abs_timefromRpeak = [];
                                    diff_cond_avg.difference(dcn).Rpeak_evoked(st, hs).abs_timeprob = [];
                                    diff_cond_avg.difference(dcn).Rpeak_evoked(st, hs).rel_timefromRpeak = [];
                                    diff_cond_avg.difference(dcn).Rpeak_evoked(st, hs).rel_timeprob = [];
                                    diff_cond_avg.difference(dcn).Rpeak_evoked(st, hs).ntrials = [];
                                end
                            end
                        end 
                    elseif strcmp(analysis, 'Rpeak_evoked_LFP_TFS')
                        if isempty(cond1_avg.hs_tuned_tfs) 
                            continue;
                        end
                        if isempty(cond2_avg.hs_tuned_tfs)
                            continue;
                        end

                        if ~isfield(cond2_avg.hs_tuned_tfs, 'powspctrm') || ~isfield(cond1_avg.hs_tuned_tfs, 'powspctrm')
                            continue;
                        end
                        % loop through handspace tunings
                        diff_cond_avg.difference(dcn).hs_tuned_tfs = cond2_avg.hs_tuned_tfs;
                        for hs = 1:size(cond2_avg.hs_tuned_tfs, 2)
                            for st = 1:size(cond2_avg.hs_tuned_tfs, 1)

                                if isfield(cond1_avg.hs_tuned_tfs(st, hs), 'powspctrm') && ...
                                        isfield(cond2_avg.hs_tuned_tfs(st, hs), 'powspctrm') && ...
                                        ~isempty(cond1_avg.hs_tuned_tfs(st, hs).powspctrm) && ...
                                    ~isempty(cond2_avg.hs_tuned_tfs(st, hs).powspctrm)
                                    ntimebins = min([size(cond2_avg.hs_tuned_tfs(st, hs).powspctrm, 3), ...
                                        size(cond1_avg.hs_tuned_tfs(st, hs).powspctrm, 3)]);
                                    % calculate the difference
                                    diff_cond_avg.difference(dcn).hs_tuned_tfs(st, hs).powspctrm = ...
                                        cond2_avg.hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins) - ...
                                        cond1_avg.hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins);
                                    diff_cond_avg.difference(dcn).hs_tuned_tfs(st, hs).time = ...
                                        cond2_avg.hs_tuned_tfs(st, hs).time(1:ntimebins);  
                                    if isfield(cond2_avg.hs_tuned_tfs(st, hs), 'ntrials')
                                        diff_cond_avg.difference(dcn).hs_tuned_tfs(st, hs).ntrials = ...
                                            [];
                                    end

                                else
                                    diff_cond_avg.difference(dcn).hs_tuned_tfs(st, hs).powspctrm = [];
                                    diff_cond_avg.difference(dcn).hs_tuned_tfs(st, hs).time = [];
                                    diff_cond_avg.difference(dcn).hs_tuned_tfs(st, hs).freq = [];
                                end
                            end
                        end
                    end
                else
                    continue;
                end
            end
        end
        if isfield(diff_cond_avg, 'difference')
            cond_avg = diff_cond_avg.difference;
        end
    end
    
    % generate return variable
    if isfield(diff_cond_avg, 'difference')
        diff_cond_avg = diff_cond_avg.difference;
    else
        diff_cond_avg = [];
    end
end

