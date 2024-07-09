function dt = ecg_bna_load_variables(cfg,unit_list, folder, data_struct, var_names)

cond_list = {cfg.condition.name};

for varNum = 1:length(var_names)
    
    if strfind(var_names{varNum}, '.')
        dotId = strfind(var_names{varNum}, '.');
        beforeDot = var_names{varNum}(1:dotId-1);
        afterDot  = var_names{varNum}(dotId+1:end);
        
        for condNum = 1:length(cfg.condition)
            L = cfg.condition(condNum).name;
            
            dt.(L).(beforeDot).(afterDot) = nan(length(unit_list),1);
        end
        
    elseif strcmp(var_names{varNum}, 'pearson_r') || strcmp(var_names{varNum}, 'pearson_p') || strcmp(var_names{varNum}, 'n_cycles') || ...
            strcmp(var_names{varNum}, 'FR_pearson_r') || strcmp(var_names{varNum}, 'FR_pearson_p') || strcmp(var_names{varNum}, 'FR_n_cycles') || ...
            strcmp(var_names{varNum}, 'RR_pearson_r') || strcmp(var_names{varNum}, 'RR_pearson_p') || strcmp(var_names{varNum}, 'RR_n_cycles')
        
        for condNum = 1:length(cfg.condition)
            L = cfg.condition(condNum).name;
            
            dt.(L).(var_names{varNum}) = nan(25,length(unit_list));
        end
        
    else
        
        for condNum = 1:length(cfg.condition)
            L = cfg.condition(condNum).name;
            
            dt.(L).(var_names{varNum}) = nan(length(unit_list),1);
        end
        
    end
    
end

folder=[cfg.SPK_root_results_fldr filesep folder filesep];

for fileNum = 1:length(unit_list)
    
    currFileInfo = dir([folder unit_list{fileNum} '_*.mat']);
    
    if isempty(currFileInfo)
        continue
    elseif ~isempty(currFileInfo)
        
        disp(['Reading file ' num2str(fileNum) ' out of ' num2str(length(unit_list))])
        A = load([folder filesep currFileInfo.name], data_struct);
        
        % drop all condition-related fields
        B = rmfield(A.(data_struct), cond_list);
        
        % move non-condition-related fields into the main data structure
        for fn = fieldnames(B)'
            if fileNum == 1
                if strcmp(fn{1}, 'criteria')
                    dt.(fn{1}).SNR_F = [];
                    dt.(fn{1}).SNR_V = [];
                    
                    dt.(fn{1}).stability_F = [];
                    dt.(fn{1}).stability_V = [];
                    
                    for subfield = fieldnames(B.criteria)'
                        dt.(fn{1}).(subfield{1})(fileNum) = B.criteria.(subfield{1});
                    end
                    
                elseif ischar(B.(fn{1}))
                    dt.(fn{1}) = {B.(fn{1})};
                else
                    dt.(fn{1}) = B.(fn{1});
                end
            else
                if strcmp(fn{1}, 'criteria')
                    for subfield = fieldnames(B.criteria)'
                        dt.(fn{1}).(subfield{1})(fileNum) = B.criteria.(subfield{1});
                    end
                elseif strcmp(fn{1}, 'thresholds_microV')
                    dt.(fn{1}) = [dt.(fn{1}) B.(fn{1})];
                elseif ischar(B.(fn{1})) || isstruct(B.(fn{1}))
                    dt.(fn{1}) = [dt.(fn{1}); B.(fn{1})];
                end
            end
            
        end
        
        for varNum = 1:length(var_names)
            for condNum = 1:length(cfg.condition)
                L = cfg.condition(condNum).name;
            
                if strfind(var_names{varNum}, '.')
                    dotId = strfind(var_names{varNum}, '.');
                    beforeDot = var_names{varNum}(1:dotId-1);
                    afterDot  = var_names{varNum}(dotId+1:end);
                    
                    var_len = length(A.(data_struct).(L).(beforeDot).(afterDot));
                    if strcmp(var_names{varNum}, 'vonMises') && var_len > 4
                        var_len = 4;
                    elseif strcmp(var_names{varNum}, 'cosine') && var_len > 3
                        var_len = 3;
                    end
                    
                    if ~isempty(A.(data_struct).(L).(beforeDot).(afterDot))
                        dt.(L).(beforeDot).(afterDot)(fileNum,1:var_len) = A.(data_struct).(L).(beforeDot).(afterDot)(1:var_len);
                    end
                    
                elseif strcmp(var_names{varNum}, 'AMP_MI')
                    dt.(L).(var_names{varNum})(fileNum, 1:5) = A.(data_struct).(L).(var_names{varNum});
                elseif strcmp(var_names{varNum}, 'cc_PSTH_feature') || strcmp(var_names{varNum}, 'pp_PSTH_feature')
                    dt.(L).(var_names{varNum})(fileNum) = A.(data_struct).(L).(var_names{varNum})(1);
                elseif strcmp(var_names{varNum}, 'NrEvents')
                    dt.(L).(var_names{varNum})(fileNum) = A.(data_struct).(L).(var_names{varNum})';
                elseif strcmp(var_names{varNum}, 'pearson_r') || strcmp(var_names{varNum}, 'pearson_p') || strcmp(var_names{varNum}, 'n_cycles') || ...
                        strcmp(var_names{varNum}, 'FR_pearson_r') || strcmp(var_names{varNum}, 'FR_pearson_p') || strcmp(var_names{varNum}, 'FR_n_cycles') || ...
                        strcmp(var_names{varNum}, 'RR_pearson_r') || strcmp(var_names{varNum}, 'RR_pearson_p') || strcmp(var_names{varNum}, 'RR_n_cycles')
                    dt.(L).(var_names{varNum})(:, fileNum) = A.(data_struct).(L).(var_names{varNum});
                else
                    if ~isempty(A.(data_struct).(L).(var_names{varNum}))
                        dt.(L).(var_names{varNum})(fileNum,:) = A.(data_struct).(L).(var_names{varNum})';
                    end
                end
            end
            
        end
        
        clear Output
    end
    
end

end