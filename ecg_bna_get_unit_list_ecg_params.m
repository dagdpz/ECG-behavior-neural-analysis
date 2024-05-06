function ecg_bna_get_unit_list_ecg_params(cfg)

list_of_lists = dir([cfg.unit_lists filesep '*.mat']);

list_of_lists_600          = dir([cfg.unit_lists filesep '*_600.mat']); % find files after exclusion
list_of_lists_noCB         = dir([cfg.unit_lists filesep '*noCB*.mat']);
list_of_lists_withCB       = dir([cfg.unit_lists filesep '*withCB.mat']);
list_of_lists_CB_corr      = dir([cfg.unit_lists filesep '*noCB_corr.mat']);
list_of_lists_CB_corr_ccs  = dir([cfg.unit_lists filesep '*noCB_corr_ccs.mat']);
list_of_lists_CB_excl      = dir([cfg.unit_lists filesep '*_excluded*.mat']);
list_of_lists_before       = dir([cfg.unit_lists filesep '*_before_exclusion.mat']);

lists2drop = ismember({list_of_lists.name}, {list_of_lists_600.name}) | ...
    ismember({list_of_lists.name}, {list_of_lists_noCB.name}) | ...
    ismember({list_of_lists.name}, {list_of_lists_withCB.name}) | ...
    ismember({list_of_lists.name}, {list_of_lists_CB_corr.name}) | ...
    ismember({list_of_lists.name}, {list_of_lists_CB_corr_ccs.name}) | ...
    ismember({list_of_lists.name}, {list_of_lists_CB_excl.name}) | ...
    ismember({list_of_lists.name}, {list_of_lists_before.name});

list_of_lists = list_of_lists(~lists2drop);

for listNum = 1:length(list_of_lists)
    
    filename = [list_of_lists(listNum).folder filesep list_of_lists(listNum).name];
    
    % stable for at least one task block OR one rest block
    load(filename, 'unit_ids', 'targets', 'sites')
    
    % figure out number of R-peaks
    dt = ecg_bna_load_variables(cfg,unit_ids, 'per_unit', 'Output', {'NrEvents'});
    
    enoughRpeaks = zeros(length(cfg.condition), length(unit_ids));
    for conNum = 1:length(cfg.condition)
        L = cfg.condition(conNum).name;
        
        enoughRpeaks(conNum,:) = dt.(L).NrEvents > cfg.spk.unit_exclusion.nCardiacCycles;
        
    end
    if strfind(list_of_lists(listNum).name, 'stableTaskAndRest')
        ids_enough_Rpeaks = all(enoughRpeaks);
    else
        ids_enough_Rpeaks = any(enoughRpeaks);
    end
    
    clear dt enoughRpeaks
                                             
    dt = ecg_bna_load_variables(cfg, unit_ids, 'cardioballistic', 'data', {'distance2thr', 'AMP_MI', 'cc_PSTH_feature', 'pp_PSTH_feature'});
    enough_bins = zeros(length(cfg.condition), length(unit_ids));
    for conNum = 1:length(cfg.condition)
        L = cfg.condition(conNum).name;
        
        enough_bins(conNum,:) = dt.(L).distance2thr > 1;
        
    end
    ids_enough_bins = all(enough_bins);
    sum(ids_enough_bins)
    
    % Considerations of exclusion of units by fit significance and the
    % goodness-of-fit parameters:
    % 1. There are fits with p < 0.05 that I consider significant: some of
    % them are true positive, some of them are not.
    % 2. I want to take (mostly) true positives and drop most false
    % positives using a threshold for R-squared
    % 3. I will consider significant those fits that have p-value <
    % 0.05 and explain considerable variance of the data (R-sq. > 0.3)
    
%     [~, h_AMP_MI] = bonf_holm([out_rest.AMP_MI(:,2), out_task.AMP_MI(:,2)]); % p < 0.05 corresponds to significant cosine fit implying a significant cardioballistic effect
    
    h_AMP_MI = zeros(length(cfg.condition), length(unit_ids));
    high_Rsq = zeros(length(cfg.condition), length(unit_ids));
    AMP_cc   = zeros(length(cfg.condition), length(unit_ids));
    for conNum = 1:length(cfg.condition)
        
        L = cfg.condition(conNum).name;
        
        h_AMP_MI(conNum,:) = dt.(L).AMP_MI(:,2) < 0.01;
        high_Rsq(conNum,:) = dt.(L).AMP_MI(:,4) > 0.3; % to drop false positives
        AMP_cc(conNum,:)   = dt.(L).pp_PSTH_feature;
        
    end
    
    % p < 0.05 - significant correlation between phase PSTH and phase
    % dynamic of AMP feature
    [~, h_AMP_cc] = bonf_holm(AMP_cc);
    
%     with_CBE = any(h_AMP_MI' & high_Rsq');
    with_CBE = any(h_AMP_MI);
    
%     no_cardioballistic_effect = ~any(h_AMP_MI' & high_Rsq'); % if a unit is affected in at least one condition we'll get rid of it
    no_cardioballistic_effect = ~any(h_AMP_MI);
    
    no_cardioballistic_effect_corrected = no_cardioballistic_effect | ids_enough_bins;
    
    no_cardioballistic_effect_corrected_ccs = (no_cardioballistic_effect | ids_enough_bins) & ~any(h_AMP_cc);
    % find unit ids and targets with enough R-peak counts
    unit_ids_600 = unit_ids(ids_enough_Rpeaks);
    targets_600  = targets(ids_enough_Rpeaks);
    sites_600    = sites(ids_enough_Rpeaks);
    % without cardioballistic effect
    unit_ids_noCB    = unit_ids(ids_enough_Rpeaks & no_cardioballistic_effect);
    targets_noCB     = targets(ids_enough_Rpeaks & no_cardioballistic_effect);
    sites_noCB       = sites(ids_enough_Rpeaks & no_cardioballistic_effect);
    % with cardioballistic effect but still having enough R-peaks for the
    % analysis
    unit_ids_withCB    = unit_ids(ids_enough_Rpeaks & ~no_cardioballistic_effect);
    targets_withCB     = targets(ids_enough_Rpeaks & ~no_cardioballistic_effect);
    sites_withCB       = sites(ids_enough_Rpeaks & ~no_cardioballistic_effect);
    % without cardioballistic effect and high spike amplitude
    unit_ids_noCB_corr = unit_ids(ids_enough_Rpeaks & no_cardioballistic_effect_corrected);
    targets_noCB_corr  = targets(ids_enough_Rpeaks & no_cardioballistic_effect_corrected);
    sites_noCB_corr    = sites(ids_enough_Rpeaks & no_cardioballistic_effect_corrected);
    % without cardioballistic effect, high spike amplitude, and with
    % non-significant correlation coefficient between phase PSTH and AMP 
    % phase dynamic
    unit_ids_noCB_corr_ccs = unit_ids(ids_enough_Rpeaks & no_cardioballistic_effect_corrected_ccs);
    targets_noCB_corr_ccs  = targets(ids_enough_Rpeaks & no_cardioballistic_effect_corrected_ccs);
    sites_noCB_corr_ccs    = sites(ids_enough_Rpeaks & no_cardioballistic_effect_corrected_ccs);
    
    % save unit list for 600
    unit_ids = unit_ids_600;
    targets  = targets_600;
    sites    = sites_600;
    save([filename(1:end-4) '_600.mat'], 'unit_ids', 'targets', 'sites')
    % noCB
    unit_ids = unit_ids_noCB;
    targets  = targets_noCB;
    sites    = sites_noCB;
    save([filename(1:end-4) '_noCB.mat'], 'unit_ids', 'targets', 'sites')
    % withCB
    unit_ids = unit_ids_withCB;
    targets  = targets_withCB;
    sites    = sites_withCB;
    save([filename(1:end-4) '_withCB.mat'], 'unit_ids', 'targets', 'sites')
    % noCB + amp criterion
    unit_ids = unit_ids_noCB_corr;
    targets  = targets_noCB_corr;
    sites    = sites_noCB_corr;
    save([filename(1:end-4) '_noCB_corr.mat'], 'unit_ids', 'targets', 'sites')
    % noCB + amp corterion + phase PSTH vs. feature dynamic correlation
    unit_ids = unit_ids_noCB_corr_ccs;
    targets  = targets_noCB_corr_ccs;
    sites    = sites_noCB_corr_ccs;
    save([filename(1:end-4) '_noCB_corr_ccs.mat'], 'unit_ids', 'targets', 'sites')
    
    % 4. save the table with unit counts
    % 600
    [unique_areas, ~, ic_600] = unique(targets_600);
    unit_counts = accumarray(ic_600,1);
    % noCB
    [~, ~, ic_noCB] = unique(targets_noCB);
    unit_counts_noCB = accumarray(ic_noCB,1);
    % withCB
    [unqTargets_withCB, ~, ic_withCB] = unique(targets_withCB);
    unit_counts_withCB = accumarray(ic_withCB,1);
    [Lia, Locb] = ismember(unique_areas, unqTargets_withCB);
    unit_counts_withCB_tmp = double(Lia);
    unit_counts_withCB_tmp(Lia) = unit_counts_withCB;
    unit_counts_withCB = unit_counts_withCB_tmp;
    % noCB + amp
    [~, ~, ic_noCB_corr] = unique(targets_noCB_corr);
    unit_counts_noCB_corr = accumarray(ic_noCB_corr,1);
    % noCB + amp + high/low cc PSTH vs. feature dynamic
    [~, ~, ic_noCB_corr_ccs] = unique(targets_noCB_corr_ccs);
    unit_counts_noCB_corr_ccs = accumarray(ic_noCB_corr_ccs,1);
    
    T = table(unique_areas, unit_counts, unit_counts_noCB, unit_counts_withCB, unit_counts_noCB_corr, unit_counts_noCB_corr_ccs);
    writetable(T, [filename(1:end-4) '_600.xls'])
    %     clear T unique_areas ic unit_counts filename unit_ids_after_exclusion targets_after_exclusion
    
end

end
