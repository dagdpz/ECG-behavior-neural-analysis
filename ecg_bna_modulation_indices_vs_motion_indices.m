function ecg_bna_modulation_indices_vs_motion_indices(cfg)

load('Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_Magnus_TaskRest_20240105\unit_lists\unitInfo_after_exclusion_stableTaskAndRest_600.mat', 'unit_ids')

perUnitDataFolder = ...
    'Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_Magnus_TaskRest_20240105\per_unit_stable_600\';
cardioballisticDataFolder = ...
    'Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_Magnus_TaskRest_20240105\cardioballistic_stable_600\';

for unitNum = 1:length(unit_ids)
    disp(['Loading unit ' num2str(unitNum) ' out of ' num2str(length(unit_ids))])
    % figure out files
    currFileName_perUnit = dir([perUnitDataFolder unit_ids{unitNum} '_*.mat']);
    currFileName_cardioballistic = dir([cardioballisticDataFolder unit_ids{unitNum} '_*.mat']);
    
    % load files, drop unnecessary fields
    load([currFileName_perUnit.folder filesep currFileName_perUnit.name], 'Output')
    if ~isempty(currFileName_cardioballistic)
        load([currFileName_cardioballistic.folder filesep currFileName_cardioballistic.name], 'data')
    end
    
    dt.unitID{unitNum} = Output.unit_ID;
    dt.target{unitNum} = Output.target;
    for c = 1:length(cfg.condition)
        L=cfg.condition(c).name;
        
        dt.(L).sig_n_bins(unitNum)           = Output.(L).sig_n_bins;
        dt.(L).FR_ModIndex_SubtrSDP(unitNum) = Output.(L).FR_ModIndex_SubtrSDP;
        dt.(L).FR_ModIndex_PcS(unitNum)      = Output.(L).FR_ModIndex_PcS;
        if exist('data', 'var')
            dt.(L).histogram_MI(unitNum)         = data.(L).histogram_MI;
            dt.(L).histogram_p(unitNum)          = data.(L).histogram_p;
            dt.(L).AMP_MI(:, unitNum)            = data.(L).AMP_MI';
            dt.(L).HW_MI(:, unitNum)             = data.(L).HW_MI';
            dt.(L).TPW_MI(:, unitNum)            = data.(L).TPW_MI';
            dt.(L).REP_MI(:, unitNum)            = data.(L).REP_MI';
        else
            dt.(L).histogram_MI(unitNum)         = NaN;
            dt.(L).histogram_p(unitNum)          = NaN;
            dt.(L).AMP_MI(:, unitNum)            = NaN(5,1);
            dt.(L).HW_MI(:, unitNum)             = NaN(5,1);
            dt.(L).TPW_MI(:, unitNum)            = NaN(5,1);
            dt.(L).REP_MI(:, unitNum)            = NaN(5,1);
        end
        
    end
    
    clear Output data
    
end

TargetBrainArea = dt.target;
Ana_TargetBrainArea = cfg.targets;

if cfg.combine_hemispheres
    TargetBrainArea     = cellfun(@(x) x(1:end-2),TargetBrainArea,'UniformOutput',false);
    Ana_TargetBrainArea = cellfun(@(x) x(1:end-2),Ana_TargetBrainArea,'UniformOutput',false);
end

Ana_TargetBrainArea = {'VPL', 'dPul', 'MD'};
N_Areas = length(Ana_TargetBrainArea);
N_conditions = 2;

%% restructure data by area and condition
for a = 1: N_Areas
    T=Ana_TargetBrainArea{a};
    for c=1:N_conditions
        L=cfg.condition(c).name;
        dat_fieldnames = fieldnames(dt.(L));
        
        for fn=1:numel(dat_fieldnames)
            N=dat_fieldnames{fn};
            CD.(T).(L).(N) = dt.(L).(N)(:, ismember(TargetBrainArea,T));
        end
        
        dat_fieldnames={'unitID','target'};
        for fn=1:numel(dat_fieldnames)
            N=dat_fieldnames{fn};
            CD.(T).(L).(N) = dt.(N)(ismember(TargetBrainArea,T));
        end
    end
end

%% modulation indices (cosine-fitted) vs. motion indices
for a = 1: N_Areas
    T=Ana_TargetBrainArea{a};
    
    figure,
    set(gcf, 'Position', [702 102 985 883])
    for c=1:N_conditions
        L=cfg.condition(c).name;
        
        % take data for all units
        curr_mod_idx    = 100 * CD.(T).(L).histogram_MI;
        tmp = CD.(T).(L).AMP_MI(:, 1);
        if ~isempty(tmp)
            curr_motion_idx = 100 * CD.(T).(L).AMP_MI(1, :);
        end
        
        [r_tmp, p_tmp] = corrcoef(curr_mod_idx, curr_motion_idx);
        r = r_tmp(2,1);
        p = p_tmp(2,1);
        clear r_tmp p_tmp
        
        subplot(2,2,c)
        scatter(curr_mod_idx, curr_motion_idx, [], cfg.condition(c).color)
        lsline
        box on
        title({[T ' (N = ' num2str(length(curr_mod_idx)) '): ' L ' - all units'], ...
            ['cc = ' num2str(r) '; p = ' num2str(p)]})
        
        if c == 1
            xlabel('Modulation Index (Cosine Fitted), %')
            ylabel('Motion Index, %')
        end
        
        % take data for only significant units without cardioballistic
        % effect
        curr_mod_p = CD.(T).(L).histogram_p;
        [~, sig_idx] = bonf_holm(curr_mod_p, 0.05);
        curr_motion_p = CD.(T).(L).AMP_MI(2, :);
        [~, sig_CB_idx] = bonf_holm(curr_motion_p, 0.05); % no cardioballistic effect
        
        curr_mod_idx    = curr_mod_idx(sig_idx & ~sig_CB_idx);
        curr_motion_idx = curr_motion_idx(sig_idx & ~sig_CB_idx);
        
        [r_tmp, p_tmp] = corrcoef(curr_mod_idx, curr_motion_idx);
        r = r_tmp(2,1);
        p = p_tmp(2,1);
        clear r_tmp p_tmp
        
        subplot(2,2,c+2)
        scatter(curr_mod_idx, curr_motion_idx, [], cfg.condition(c).color)
        lsline
        box on
        title({'Significantly modulated units without cardioballistic effect', ...
            [T ' (N = ' num2str(length(curr_mod_idx)) '): ' L], ...
            ['cc = ' num2str(r) '; p = ' num2str(p)]})
        
        if c == 1
            xlabel('Modulation Index (Cosine Fitted), %')
            ylabel('Motion Index, %')
        end
        
    end
end

%% modulation indices vs. motion indices
for a = 1: N_Areas
    T=Ana_TargetBrainArea{a};
    
    figure,
    for c=1:N_conditions
        L=cfg.condition(c).name;
        
        % take data for all units
        curr_mod_idx = CD.(T).(L).FR_ModIndex_PcS;
        tmp = CD.(T).(L).AMP_MI(:, 1);
        if ~isempty(tmp)
            curr_motion_idx = 100 * CD.(T).(L).AMP_MI(1, :);
        end
        
        [r_tmp, p_tmp] = corrcoef(curr_mod_idx, curr_motion_idx);
        r = r_tmp(2,1);
        p = p_tmp(2,1);
        clear r_tmp p_tmp
        
        subplot(2,2,c)
        scatter(curr_mod_idx, curr_motion_idx, [], cfg.condition(c).color)
        lsline
        box on
        title({[T ' (N = ' num2str(length(curr_mod_idx)) '): ' L ' - all units'], ...
            ['cc = ' num2str(r) '; p = ' num2str(p)]})
        
        if c == 1
            xlabel('Modulation Index, %')
            ylabel('Motion Index, %')
        end
        
        % take data for only significant units without cardioballistic
        % effect
        curr_mod_sig = CD.(T).(L).sig_n_bins;
        sig_idx = curr_mod_sig > 4;
        curr_motion_p = CD.(T).(L).AMP_MI(2, :);
        [~, sig_CB_idx] = bonf_holm(curr_motion_p, 0.05); % no cardioballistic effect
        
        curr_mod_idx    = curr_mod_idx(sig_idx & ~sig_CB_idx);
        curr_motion_idx = curr_motion_idx(sig_idx & ~sig_CB_idx);
        
        if length(curr_mod_idx) > 2
            [r_tmp, p_tmp] = corrcoef(curr_mod_idx, curr_motion_idx);
            r = r_tmp(2,1);
            p = p_tmp(2,1);
            clear r_tmp p_tmp
        else
            r = NaN;
            p = NaN;
        
        subplot(2,2,c)
        scatter(curr_mod_idx, curr_motion_idx, [], cfg.condition(c).color)
        lsline
        box on
        title({'Significantly modulated units without cardioballistic effect', ...
            [T ' (N = ' num2str(length(curr_mod_idx)) '): ' L], ...
            ['cc = ' num2str(r) '; p = ' num2str(p)]})
        
        if c == 1
            xlabel('Modulation Index, %')
            ylabel('Motion Index, %')
        end
        
    end
end


end