
% sort the data according to brain areas
for f_brain = 1: length(Ana_TargetBrainArea)
    T=Ana_TargetBrainArea{f_brain};
    
    idx_brain = cellfun(@(x) contains(lower(x.target), lower(T)), SPK_PSTH, 'Uniformoutput', false);
    
    for c=1:numel(cfg.condition)
        L=cfg.condition(c).name;
        unit_ID_tmp = cellfun(@(x,y) x.unit_ID(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).unit_ID = [unit_ID_tmp{:}]';
        
        target_tmp = cellfun(@(x,y) x.target(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).target = [target_tmp{:}]';
        
        quantSNR_tmp = cellfun(@(x,y) x.quantSNR(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).quantSNR = [quantSNR_tmp{:}]';
        
        single_rating_tmp = cellfun(@(x,y) x.Single_rating(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).Single_rating = [single_rating_tmp{:}]';
        
        stability_rating_tmp = cellfun(@(x,y) x.stability_rating(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).stability_rating = [stability_rating_tmp{:}]';
        
        SD_tmp = cellfun(@(x,y) x.(L).SD(y,:)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).SD = [SD_tmp{:}]';
        
        SDP_tmp = cellfun(@(x,y) x.(L).SDP(y,:)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).SDP = [SDP_tmp{:}]';
        
        SD_SEM_tmp = cellfun(@(x,y) x.(L).SD_SEM(y,:)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).SD_SEM = [SD_SEM_tmp{:}]';
        
        sig_tmp = cellfun(@(x,y) x.(L).sig(y,:)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).sig_period = [sig_tmp{:}]';
        
        sig_FR_diff_tmp = cellfun(@(x,y) x.(L).sig_FR_diff(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).sig_FR_diff = [sig_FR_diff_tmp{:}]';
        
        sig_n_bins_tmp = cellfun(@(x,y) x.(L).sig_n_bins(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).sig_n_bins = [sig_n_bins_tmp{:}]';
        
        sig_sign_tmp = cellfun(@(x,y) x.(L).sig_sign(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).sig_sign = [sig_sign_tmp{:}]';
        
        sig_time_tmp = cellfun(@(x,y) x.(L).sig_time(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).sig_time = [sig_time_tmp{:}]';
        
        FR_tmp = cellfun(@(x,y) x.(L).FR(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).FR = [FR_tmp{:}]';
        
        NrEvents_tmp = cellfun(@(x,y) x.(L).NrEvents(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).NrEvents = [NrEvents_tmp{:}]';
        
        NrTrials_tmp = cellfun(@(x,y) x.(L).NrTrials(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).NrTrials = [NrTrials_tmp{:}]';
        
        SDsubstractedSDP_tmp = cellfun(@(x,y) x.(L).SDsubstractedSDP(y,:)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).SDsubstractedSDP = [SDsubstractedSDP_tmp{:}]';
        
        SDsubstractedSDP_normalized_tmp = cellfun(@(x,y) x.(L).SDsubstractedSDP_normalized(y,:)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).SDsubstractedSDP_normalized = [SDsubstractedSDP_normalized_tmp{:}]';
        
        FR_ModIndex_SubtrSDP_tmp = cellfun(@(x,y) x.(L).FR_ModIndex_SubtrSDP(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).FR_ModIndex_SubtrSDP = [FR_ModIndex_SubtrSDP_tmp{:}]';
        
        FR_ModIndex_PcS_tmp = cellfun(@(x,y) x.(L).FR_ModIndex_PcS(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(T).(L).FR_ModIndex_PcS = [FR_ModIndex_PcS_tmp{:}]';
        
    end
end


%% Crieria for the dataset
% if OnlyUnits_withRestANDTask
%     for i_BrArea = 1: length(fieldnames(Out))
%         Idx_Units_RestTask = []; Idx_NoSpikes_T = [];  Idx_NoSpikes_R = [];
%         
%         %  Task = NAN and Rest = NAN - should not exist?
%         Idx_NoSpikes_R = find(any(isnan([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{1}).SD]),2));
%         Idx_NoSpikes_T = find(any(isnan([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{2}).SD]),2));
%         Idx_NoSpikes_RestTask = intersect( Idx_NoSpikes_T,Idx_NoSpikes_R);
%         Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{1}).unit_ID(Idx_NoSpikes_RestTask,:)       %What happenend here?
%         
%         Idx_Spikes_R = ~any(isnan([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{1}).SD]),2);
%         Idx_Spikes_T = ~any(isnan([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{2}).SD]),2);
%         Idx_Units_RestTask = ~(Idx_Spikes_T & Idx_Spikes_R);
%         
%         disp([(Ana_TargetBrainArea{i_BrArea}),': ', num2str(length(Idx_Units_RestTask))])
%         for c=1:numel(cfg.condition)
%             L=cfg.condition(c).name;
%             if i_tsk == 1;  Idx = Idx_NoSpikes_R ; else Idx = Idx_NoSpikes_T ; end
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).SD(Idx, :)           = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).SD_SEM(Idx, :)       = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).SDP(Idx, :)          = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).sig_period(Idx, :)   = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).FR(Idx)              = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).quantSNR(Idx)        = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).Single_rating(Idx)   = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).stability_rating(Idx)= [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).sig_n_bins(Idx)      = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).sig_time(Idx)        = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).sig_sign(Idx)        = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).sig_FR_diff(Idx)     = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).unit_ID(Idx)         = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).target(Idx)          = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).SDsubstractedSDP(Idx, :)            = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).SDsubstractedSDP_normalized(Idx, :) = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).FR_ModIndex_SubtrSDP(Idx)           = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).FR_ModIndex_PcS(Idx)                = [] ;
%         end
%         
%         Idx_Spikes_R = find(~any(isnan([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{1}).SD]),2));
%         Idx_Spikes_T = find(~any(isnan([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{2}).SD]),2));
%         Idx_Units_RestTask = intersect(Idx_Spikes_R, Idx_Spikes_T);
%         disp([(Ana_TargetBrainArea{i_BrArea}),': ', num2str(length(Idx_Units_RestTask))])
%         
%     end
% end
%% Exclusion criteria - here we exclude only by the number of R-peaks
% for i_BrArea = 1: length(fieldnames(Out))
%     
%     fieldList = fieldnames(Out.(Ana_TargetBrainArea{i_BrArea}).Rest);
%     
%     ids_included = ...
%         [Out.(Ana_TargetBrainArea{i_BrArea}).Rest.NrEvents] > ecg_bna_cfg.unit_exclusion.nCardiacCycles & ... % find INVALID units with below # of heart cycles
%         [Out.(Ana_TargetBrainArea{i_BrArea}).Task.NrEvents] > ecg_bna_cfg.unit_exclusion.nCardiacCycles;
%     
%     for c=1:numel(cfg.condition)
%         L=cfg.condition(c).name;
%         
%         for fieldNum = 1:length(fieldList)
%             
%             if strcmp(fieldList{fieldNum}, 'unit_ID') || strcmp(fieldList{fieldNum}, 'target')
%                 Out_excluded.(Ana_TargetBrainArea{i_BrArea}).(L).(fieldList{fieldNum}) = ...
%                     Out.(Ana_TargetBrainArea{i_BrArea}).(L).(fieldList{fieldNum})(~ids_included);
%                 Out.(Ana_TargetBrainArea{i_BrArea}).(L).(fieldList{fieldNum}) = ...
%                     Out.(Ana_TargetBrainArea{i_BrArea}).(L).(fieldList{fieldNum})(ids_included);
%             else
%                 Out_excluded.(Ana_TargetBrainArea{i_BrArea}).(L).(fieldList{fieldNum}) = ...
%                     Out.(Ana_TargetBrainArea{i_BrArea}).(L).(fieldList{fieldNum})(~ids_included,:);
%                 Out.(Ana_TargetBrainArea{i_BrArea}).(L).(fieldList{fieldNum}) = ...
%                     Out.(Ana_TargetBrainArea{i_BrArea}).(L).(fieldList{fieldNum})(ids_included,:);
%             end
%         end
%         
%     end
%     
% end


%% Selection criteria - their values are in the ecg_bna_cfg now
% for i_BrArea = 1: length(fieldnames(Out))
%     NanUnits_BeforeExclusion_idx = ...
%         isnan([Out.(Ana_TargetBrainArea{i_BrArea}).Rest.SD(:,1) ...
%         Out.(Ana_TargetBrainArea{i_BrArea}).Task.SD(:,1)]);
%     
%     NanUnits_BeforeExclusion_idx = all(NanUnits_BeforeExclusion_idx,2);
%     
%     for c=1:numel(cfg.condition)
%         L=cfg.condition(c).name;
%         
%         fieldList = fieldnames(Out.(Ana_TargetBrainArea{i_BrArea}).(L));
%         for fieldNum = 1:length(fieldList)
%             if strcmp(fieldList{fieldNum}, 'unit_ID') || strcmp(fieldList{fieldNum}, 'target')
%                 Out.(Ana_TargetBrainArea{i_BrArea}).(L).(fieldList{fieldNum}) = ...
%                     Out.(Ana_TargetBrainArea{i_BrArea}).(L).(fieldList{fieldNum})(~NanUnits_BeforeExclusion_idx);
%             else
%                 Out.(Ana_TargetBrainArea{i_BrArea}).(L).(fieldList{fieldNum}) = ...
%                     Out.(Ana_TargetBrainArea{i_BrArea}).(L).(fieldList{fieldNum})(~NanUnits_BeforeExclusion_idx,:);
%             end
%         end
%     end
% end
% 
% for i_BrArea = 1: length(fieldnames(Out))
%     disp([Ana_TargetBrainArea{i_BrArea} ' Units before exclusion criteria '])
%     unitList_beforeExclusion.(Ana_TargetBrainArea{i_BrArea}) = ...
%         union(Out.(Ana_TargetBrainArea{i_BrArea}).Rest.unit_ID, Out.(Ana_TargetBrainArea{i_BrArea}).Task.unit_ID);
%     length(unitList_beforeExclusion.(Ana_TargetBrainArea{i_BrArea}))
% end
% 
% disp('Saved list of units before exclusion')
% save([basepath_to_save filesep 'units_before_exclusion.mat'], 'unitList_beforeExclusion')
% 
% % Tab_ExcludedUnits = [];
% for i_BrArea = 1: length(fieldnames(Out))
%     InVal_unit_ID_Rest = [];
%     for c=1:numel(cfg.condition)
%         L=cfg.condition(c).name;
%         NanUnits_BeforeExclusion_idx = isnan(Out.(Ana_TargetBrainArea{i_BrArea}).(L).SD(:,1) );
%         NanUnits_BeforeExclusion =  sum(NanUnits_BeforeExclusion_idx);
%         
%         %         fieldList = fieldnames(Out.(Ana_TargetBrainArea{i_BrArea}).(L));
%         %         for fieldNum = 1:length(fieldList)
%         %             Out.(Ana_TargetBrainArea{i_BrArea}).(L).(fieldList{fieldNum}) = ...
%         %                 Out.(Ana_TargetBrainArea{i_BrArea}).(L).(fieldList{fieldNum})(~NanUnits_BeforeExclusion_idx,:);
%         %         end
%         
%         ValidUnits_BeforeExclusion =  sum(~NanUnits_BeforeExclusion_idx);
%         
%         InVal_idx1 = [Out.(Ana_TargetBrainArea{i_BrArea}).(L).FR] > ecg_bna_cfg.unit_exclusion.FR_thresholds(2) | ...
%             [Out.(Ana_TargetBrainArea{i_BrArea}).(L).FR] < ecg_bna_cfg.unit_exclusion.FR_thresholds(1); % find INVALID units with FR below FR threshold
%         
%         InVal_idx2 = [Out.(Ana_TargetBrainArea{i_BrArea}).(L).NrEvents] < ecg_bna_cfg.unit_exclusion.nCardiacCycles; % find INVALID units with below # of heart cycles
%         
%         InVal_idx3 = ...
%             [Out.(Ana_TargetBrainArea{i_BrArea}).(L).quantSNR] > ecg_bna_cfg.unit_exclusion.SNR_thresholds(2) | ...
%             [Out.(Ana_TargetBrainArea{i_BrArea}).(L).quantSNR] < ecg_bna_cfg.unit_exclusion.SNR_thresholds(1); % find INVALID units with SNR below threshold
%         
%         InVal_idx4 = ...
%             [Out.(Ana_TargetBrainArea{i_BrArea}).(L).stability_rating] > ecg_bna_cfg.unit_exclusion.FR_stability_thresholds(2) | ...
%             [Out.(Ana_TargetBrainArea{i_BrArea}).(L).stability_rating] < ecg_bna_cfg.unit_exclusion.FR_stability_thresholds(1); % find INVALID units with FR stability below threshold
%         
% %         InVal_idx = InVal_idx1 | InVal_idx2 | InVal_idx3 | InVal_idx4;
%         InVal_idx = InVal_idx1 | InVal_idx2 | InVal_idx3;
%         
%         if sum(InVal_idx)
%             InVal_unit_ID = Out.(Ana_TargetBrainArea{i_BrArea}).(L).unit_ID(InVal_idx);
%             if i_tsk == 1  % task
%                 InVal_unit_ID_Rest = InVal_unit_ID;
%                 idx_IdentRest_Task = 0;
%                 InVal_Nr = repmat(length(InVal_unit_ID), length(InVal_idx), 1);
%                 FR = Out.(Ana_TargetBrainArea{i_BrArea}).(L).FR(InVal_idx);
%             else % rest
%                 idx_IdentRest_Task =  ismember(InVal_unit_ID, InVal_unit_ID_Rest);
%                 InVal_Nr = repmat(length(InVal_unit_ID(~idx_IdentRest_Task)), length(InVal_idx), 1);
%                 FR = Out.(Ana_TargetBrainArea{i_BrArea}).(L).FR(InVal_idx);
%                 
%             end
%             
%             TaskType = repmat((TaskTyp(i_tsk)), length(InVal_idx), 1);
%             %             if strcmp((Ana_TargetBrainArea{i_BrArea}), 'VPL_R')||   strcmp((Ana_TargetBrainArea{i_BrArea}), 'mdT_L') ||   strcmp((Ana_TargetBrainArea{i_BrArea}), 'mdT_R')
%             %                 BrainArea = repmat({[ '_', (Ana_TargetBrainArea{i_BrArea})]}, length(InVal_idx), 1);
%             %             else
%             BrainArea = repmat((Ana_TargetBrainArea(i_BrArea)), length(InVal_idx), 1);
%             %             end
%             %             Criterium_SpkPerSec = repmat(Criterium_SpkPerSec, length(InVal_idx), 1);
%             %             Criterium_NrCardiacCycles = repmat(Criterium_NrCardiacCycles, length(InVal_idx), 1);
%             Nr_InVal_idx2 = repmat(length(InVal_idx2), length(InVal_idx), 1);
%             Nr_InVal_idx1 = repmat(length(InVal_idx1), length(InVal_idx), 1);
%             NanUnits_BeforeExclusion = repmat(NanUnits_BeforeExclusion, length(InVal_idx), 1);
%             ValidUnits_BeforeExclusion = repmat(ValidUnits_BeforeExclusion, length(InVal_idx), 1);
%             
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).quantSNR(InVal_idx)        = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).Single_rating(InVal_idx)   = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).stability_rating(InVal_idx)= [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).unit_ID(InVal_idx)         = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).target(InVal_idx)          = [] ;
%             
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).SD(InVal_idx, :)           = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).SD_SEM(InVal_idx, :)       = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).SDP(InVal_idx, :)          = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).sig_period(InVal_idx, :)   = [] ;
%             
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).FR(InVal_idx)              = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).sig_n_bins(InVal_idx)      = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).sig_time(InVal_idx)        = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).sig_sign(InVal_idx)        = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).sig_FR_diff(InVal_idx)     = [] ;
%             
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).NrEvents(InVal_idx)        = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).NrTrials(InVal_idx)        = [] ;
%             
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).SDsubstractedSDP(InVal_idx, :)            = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).SDsubstractedSDP_normalized(InVal_idx, :) = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).FR_ModIndex_SubtrSDP(InVal_idx)           = [] ;
%             Out.(Ana_TargetBrainArea{i_BrArea}).(L).FR_ModIndex_PcS(InVal_idx)                = [] ;
%             
%             NanUnits_AfterExclusion    =  sum(isnan(Out.(Ana_TargetBrainArea{i_BrArea}).(L).SD(:,1) ));
%             ValidUnits_AfterExclusion  =  sum(~isnan(Out.(Ana_TargetBrainArea{i_BrArea}).(L).SD(:,1) ));
%             NanUnits_AfterExclusion    = repmat(NanUnits_AfterExclusion, length(InVal_idx), 1);
%             ValidUnits_AfterExclusion  = repmat(ValidUnits_AfterExclusion, length(InVal_idx), 1);
%             
%             %             ExcludedUnits = table(NanUnits_BeforeExclusion,NanUnits_AfterExclusion, ValidUnits_BeforeExclusion, ValidUnits_AfterExclusion, Criterium_SpkPerSec,Nr_InVal_idx1, FR, Criterium_NrCardiacCycles, Nr_InVal_idx2, TaskType,BrainArea, InVal_unit_ID,InVal_idx,InVal_Nr  );
%             %             Tab_ExcludedUnits = [Tab_ExcludedUnits;ExcludedUnits ];
%         end
%     end
% end
% 
% for i_BrArea = 1: length(fieldnames(Out))
%     disp([Ana_TargetBrainArea{i_BrArea} ' Units after exclusion criteria '])
%     unitList_afterExclusion.(Ana_TargetBrainArea{i_BrArea}) = ...
%         union(Out.(Ana_TargetBrainArea{i_BrArea}).Rest.unit_ID, Out.(Ana_TargetBrainArea{i_BrArea}).Task.unit_ID);
%     length(unitList_afterExclusion.(Ana_TargetBrainArea{i_BrArea}))
% end
% 
% disp('Saved list of units after exclusion')
% save([basepath_to_save filesep 'units_after_exclusion.mat'], 'unitList_afterExclusion')

% if saveTable == 1
% save([basepath_to_save,filesep , 'Table_excludedUnits' ],'Tab_ExcludedUnits');
% filename = [basepath_to_save,filesep , 'Table_excludedUnits2.xlsx' ];
% % writetable(Tab_ExcludedUnits,filename,'Sheet',1,  'Range' ,'A1' )
% disp(['SAVED   ', basepath_to_save,filesep , 'Table_excludedUnits' ])
% end


%% overview about all units - change to include the
ThreeTiming = {'T<-50', '-50>T<50', 'T>50'} ;

hf = figure('Name',sprintf('BarPlot'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
for a = 1: N_Areas
    for c=1:N_conditions
        L=cfg.condition(c).name;
        for i_Time = 1: length(ThreeTiming)
            
            
            out = [Out.(Ana_TargetBrainArea{a}).(L)];
            % from all nicht NAN units - how many were significant?
            Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
            Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
            
            idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
            idx_SigDec      = (out.sig_sign == -1);
            idx_SigInc      = (out.sig_sign == 1);
            idx_SigTime_BeforeMinus50 = (out.sig_time < -50 );
            idx_SigTime_Around0       = (out.sig_time > -50 ) & (out.sig_time < 50) ;
            idx_SigTime_After50       = (out.sig_time > 50);
            idx_Time = [];
            switch i_Time
                case 1 %{ 'BeforeMinus50'}
                    idx_Time_Before = (idx_SigTime_BeforeMinus50 & idx_SigInc)+  (idx_SigTime_After50 & idx_SigDec) ;
                    
                case 2 %'Around0'
                    idx_Time_During = idx_SigTime_Around0 ;
                    
                case 3 %'After50'
                    idx_Time_End = (idx_SigTime_BeforeMinus50 & idx_SigDec)+  (idx_SigTime_After50 & idx_SigInc) ;
                    
            end
            
            if i_tsk == 1
                Pc_SignFR_rest(a,:) = ([sum(out.sig_sign(idx_sig & idx_Time) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] / sum(Idx_Units_NonNaN)) *100;
                Nb_SignFR_rest(a,:) = ([sum(out.sig_sign(idx_sig) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] ) ;
                Pc_SignFR_rest2(a,:) = round(([sum(idx_sig), (sum(~idx_sig) -Idx_Units_NaN)] / sum(Idx_Units_NonNaN)) *100);
                % sum(Nb_SignFR_rest, 2)
            else
                Pc_SignFR_task(a,:) = ([sum(out.sig_sign(idx_sig) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] / sum(Idx_Units_NonNaN)) *100;
                Nb_SignFR_task(a,:) = ([sum(out.sig_sign(idx_sig) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] ) ;
                Pc_SignFR_task2(a,:) = round(([sum(idx_sig),(sum(~idx_sig) -Idx_Units_NaN),] / sum(Idx_Units_NonNaN)) *100);
                % sum(Nb_SignFR_task, 2)
            end
            
        end
    end
end
ha1 = subplot(1,2,1);% pie plot how many
barpairs =  [Pc_SignFR_rest];
b = bar(barpairs,'stacked', 'Facecolor','flat' );
title('Rest: non-sig.yellow,iFR-blue,dFR-green','interpreter','none');
set(gca,'XTickLabel',fieldnames(Out),'fontsize',10);

ha1 = subplot(1,2,2);% pie plot how many
barpairs =  [Pc_SignFR_task];
b = bar(barpairs,'stacked', 'Facecolor','flat' );
title('Task: non-sig.yellow,iFR-blue,dFR-green','interpreter','none');
set(gca,'XTickLabel',fieldnames(Out),'fontsize',10);




filename= 'Pc_CardiacRelatedUnits';
if savePlot;
    export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
    close all;
end
