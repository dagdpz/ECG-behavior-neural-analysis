function ecg_bna_avg_spike_histogram(SPK_PSTH,session_info, ecg_bna_cfg)
% Here comes some sort of across population plot i assume?

savePlot = 1;
saveTable= 0;
OnlyUnits_withRestANDTask = 0;
Graph_SelectionCriterion = 1;
colors = distinguishable_colors(25);

basepath_to_save=[session_info(1).SPK_fldr filesep 'Population'];
if ~exist(basepath_to_save,'dir')
    mkdir(basepath_to_save);
end

TaskTyp = {'Rest', 'Task'};
condition_colors = {[0 0 1], [1 0 0]};
lineProps = {{'color','b','linewidth',4}, ...
    {'color','r','linewidth',4}};

TargetBrainArea = cellfun(@(x) [x.target]', SPK_PSTH, 'Uniformoutput', false);
TargetBrainArea = [TargetBrainArea{:}];
Ana_TargetBrainArea = unique(TargetBrainArea);
mdt_idx = cellfun(@(x) strcmp(x, 'mdT_R') || strcmp(x, 'mdT_L'), Ana_TargetBrainArea);
Ana_TargetBrainArea(mdt_idx) = deal({'MD'});
dpul_idx = cellfun(@(x) strcmp(x, 'dPul_R') || strcmp(x, 'dPul_L'), Ana_TargetBrainArea);
Ana_TargetBrainArea(dpul_idx) = deal({'dPul'});
% Ana_TargetBrainArea = unique(Ana_TargetBrainArea);

Ana_TargetBrainArea = {'VPL', 'dPul', 'MD'};

%% Create function to concatenate the variables per TargetBrainArea
Out = [];

% sort the data according to brain areas
for f_brain = 1: length(Ana_TargetBrainArea)
    
    idx_brain = cellfun(@(x) contains(lower(x.target), lower(Ana_TargetBrainArea{f_brain})), SPK_PSTH, 'Uniformoutput', false);
    
    for i_tsk = 1: numel(TaskTyp)
        
        unit_ID_tmp = cellfun(@(x,y) x.unit_ID(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).unit_ID = [unit_ID_tmp{:}]';
        
        target_tmp = cellfun(@(x,y) x.target(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).target = [target_tmp{:}]';
        
        quantSNR_tmp = cellfun(@(x,y) x.quantSNR(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).quantSNR = [quantSNR_tmp{:}]';
        
        single_rating_tmp = cellfun(@(x,y) x.Single_rating(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).Single_rating = [single_rating_tmp{:}]';
        
        stability_rating_tmp = cellfun(@(x,y) x.stability_rating(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).stability_rating = [stability_rating_tmp{:}]';
        
        SD_tmp = cellfun(@(x,y) x.(TaskTyp{i_tsk}).SD(y,:)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SD = [SD_tmp{:}]';
        
        SDP_tmp = cellfun(@(x,y) x.(TaskTyp{i_tsk}).SDP(y,:)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SDP = [SDP_tmp{:}]';
        
        SD_SEM_tmp = cellfun(@(x,y) x.(TaskTyp{i_tsk}).SD_SEM(y,:)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SD_SEM = [SD_SEM_tmp{:}]';
        
        sig_tmp = cellfun(@(x,y) x.(TaskTyp{i_tsk}).sig(y,:)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_period = [sig_tmp{:}]';
        
        sig_FR_diff_tmp = cellfun(@(x,y) x.(TaskTyp{i_tsk}).sig_FR_diff(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_FR_diff = [sig_FR_diff_tmp{:}]';
        
        sig_n_bins_tmp = cellfun(@(x,y) x.(TaskTyp{i_tsk}).sig_n_bins(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_n_bins = [sig_n_bins_tmp{:}]';
        
        sig_sign_tmp = cellfun(@(x,y) x.(TaskTyp{i_tsk}).sig_sign(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_sign = [sig_sign_tmp{:}]';
        
        sig_time_tmp = cellfun(@(x,y) x.(TaskTyp{i_tsk}).sig_time(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_time = [sig_time_tmp{:}]';
        
        FR_tmp = cellfun(@(x,y) x.(TaskTyp{i_tsk}).FR(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).FR = [FR_tmp{:}]';
        
        NrEvents_tmp = cellfun(@(x,y) x.(TaskTyp{i_tsk}).NrEvents(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).NrEvents = [NrEvents_tmp{:}]';
        
        NrTrials_tmp = cellfun(@(x,y) x.(TaskTyp{i_tsk}).NrTrials(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).NrTrials = [NrTrials_tmp{:}]';
        
        SDsubstractedSDP_tmp = cellfun(@(x,y) x.(TaskTyp{i_tsk}).SDsubstractedSDP(y,:)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SDsubstractedSDP = [SDsubstractedSDP_tmp{:}]';
        
        SDsubstractedSDP_normalized_tmp = cellfun(@(x,y) x.(TaskTyp{i_tsk}).SDsubstractedSDP_normalized(y,:)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized = [SDsubstractedSDP_normalized_tmp{:}]';
        
        FR_ModIndex_SubtrSDP_tmp = cellfun(@(x,y) x.(TaskTyp{i_tsk}).FR_ModIndex_SubtrSDP(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).FR_ModIndex_SubtrSDP = [FR_ModIndex_SubtrSDP_tmp{:}]';
        
        FR_ModIndex_PcS_tmp = cellfun(@(x,y) x.(TaskTyp{i_tsk}).FR_ModIndex_PcS(y)', SPK_PSTH, idx_brain, 'Uniformoutput', false);
        Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).FR_ModIndex_PcS = [FR_ModIndex_PcS_tmp{:}]';
        
    end
end


%% Crieria for the dataset
if OnlyUnits_withRestANDTask
    for i_BrArea = 1: length(fieldnames(Out))
        Idx_Units_RestTask = []; Idx_NoSpikes_T = [];  Idx_NoSpikes_R = [];
        
        %  Task = NAN and Rest = NAN - should not exist?
        Idx_NoSpikes_R = find(any(isnan([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{1}).SD]),2));
        Idx_NoSpikes_T = find(any(isnan([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{2}).SD]),2));
        Idx_NoSpikes_RestTask = intersect( Idx_NoSpikes_T,Idx_NoSpikes_R);
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{1}).unit_ID(Idx_NoSpikes_RestTask,:)       %What happenend here?
        
        Idx_Spikes_R = ~any(isnan([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{1}).SD]),2);
        Idx_Spikes_T = ~any(isnan([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{2}).SD]),2);
        Idx_Units_RestTask = ~(Idx_Spikes_T & Idx_Spikes_R);
        
        disp([(Ana_TargetBrainArea{i_BrArea}),': ', num2str(length(Idx_Units_RestTask))])
        for i_tsk = 1: numel(TaskTyp)
            if i_tsk == 1;  Idx = Idx_NoSpikes_R ; else Idx = Idx_NoSpikes_T ; end
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD(Idx, :)           = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD_SEM(Idx, :)       = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDP(Idx, :)          = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_period(Idx, :)   = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR(Idx)              = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).quantSNR(Idx)        = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).Single_rating(Idx)   = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).stability_rating(Idx)= [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_n_bins(Idx)      = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_time(Idx)        = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_sign(Idx)        = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_FR_diff(Idx)     = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).unit_ID(Idx)         = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).target(Idx)          = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP(Idx, :)            = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized(Idx, :) = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_SubtrSDP(Idx)           = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_PcS(Idx)                = [] ;
        end
        
        Idx_Spikes_R = find(~any(isnan([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{1}).SD]),2));
        Idx_Spikes_T = find(~any(isnan([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{2}).SD]),2));
        Idx_Units_RestTask = intersect(Idx_Spikes_R, Idx_Spikes_T);
        disp([(Ana_TargetBrainArea{i_BrArea}),': ', num2str(length(Idx_Units_RestTask))])
        
    end
end
%% Selection criteria - their values are in the ecg_bna_cfg now
for i_BrArea = 1: length(fieldnames(Out))
    NanUnits_BeforeExclusion_idx = ...
        isnan([Out.(Ana_TargetBrainArea{i_BrArea}).Rest.SD(:,1) ...
        Out.(Ana_TargetBrainArea{i_BrArea}).Task.SD(:,1)]);
    
    NanUnits_BeforeExclusion_idx = all(NanUnits_BeforeExclusion_idx,2);
    
    for i_tsk = 1: numel(TaskTyp)
        
        fieldList = fieldnames(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}));
        for fieldNum = 1:length(fieldList)
            if strcmp(fieldList{fieldNum}, 'unit_ID') || strcmp(fieldList{fieldNum}, 'target')
                Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).(fieldList{fieldNum}) = ...
                    Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).(fieldList{fieldNum})(~NanUnits_BeforeExclusion_idx);
            else
                Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).(fieldList{fieldNum}) = ...
                    Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).(fieldList{fieldNum})(~NanUnits_BeforeExclusion_idx,:);
            end
        end
    end
end

for i_BrArea = 1: length(fieldnames(Out))
    disp([Ana_TargetBrainArea{i_BrArea} ' Units before exclusion criteria '])
    unitList_beforeExclusion.(Ana_TargetBrainArea{i_BrArea}) = ...
        union(Out.(Ana_TargetBrainArea{i_BrArea}).Rest.unit_ID, Out.(Ana_TargetBrainArea{i_BrArea}).Task.unit_ID);
    length(unitList_beforeExclusion.(Ana_TargetBrainArea{i_BrArea}))
end

disp('Saved list of units before exclusion')
save([basepath_to_save filesep 'units_before_exclusion.mat'], 'unitList_beforeExclusion')

% Tab_ExcludedUnits = [];
for i_BrArea = 1: length(fieldnames(Out))
    InVal_unit_ID_Rest = [];
    for i_tsk = 1: numel(TaskTyp)
        
        NanUnits_BeforeExclusion_idx = isnan(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD(:,1) );
        NanUnits_BeforeExclusion =  sum(NanUnits_BeforeExclusion_idx);
        
        %         fieldList = fieldnames(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}));
        %         for fieldNum = 1:length(fieldList)
        %             Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).(fieldList{fieldNum}) = ...
        %                 Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).(fieldList{fieldNum})(~NanUnits_BeforeExclusion_idx,:);
        %         end
        
        ValidUnits_BeforeExclusion =  sum(~NanUnits_BeforeExclusion_idx);
        
        InVal_idx1 = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR] > ecg_bna_cfg.unit_exclusion.FR_thresholds(2) | ...
            [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR] < ecg_bna_cfg.unit_exclusion.FR_thresholds(1); % find INVALID units with FR below FR threshold
        
        InVal_idx2 = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).NrEvents] < ecg_bna_cfg.unit_exclusion.nCardiacCycles; % find INVALID units with below # of heart cycles
        
        InVal_idx3 = ...
            [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).quantSNR] > ecg_bna_cfg.unit_exclusion.SNR_thresholds(2) | ...
            [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).quantSNR] < ecg_bna_cfg.unit_exclusion.SNR_thresholds(1); % find INVALID units with SNR below threshold
        
        InVal_idx4 = ...
            [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).stability_rating] > ecg_bna_cfg.unit_exclusion.FR_stability_thresholds(2) | ...
            [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).stability_rating] < ecg_bna_cfg.unit_exclusion.FR_stability_thresholds(1); % find INVALID units with FR stability below threshold
        
%         InVal_idx = InVal_idx1 | InVal_idx2 | InVal_idx3 | InVal_idx4;
        InVal_idx = InVal_idx1 | InVal_idx2 | InVal_idx3;
        
        if sum(InVal_idx)
            InVal_unit_ID = Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).unit_ID(InVal_idx);
            if i_tsk == 1  % task
                InVal_unit_ID_Rest = InVal_unit_ID;
                idx_IdentRest_Task = 0;
                InVal_Nr = repmat(length(InVal_unit_ID), length(InVal_idx), 1);
                FR = Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR(InVal_idx);
            else % rest
                idx_IdentRest_Task =  ismember(InVal_unit_ID, InVal_unit_ID_Rest);
                InVal_Nr = repmat(length(InVal_unit_ID(~idx_IdentRest_Task)), length(InVal_idx), 1);
                FR = Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR(InVal_idx);
                
            end
            
            TaskType = repmat((TaskTyp(i_tsk)), length(InVal_idx), 1);
            %             if strcmp((Ana_TargetBrainArea{i_BrArea}), 'VPL_R')||   strcmp((Ana_TargetBrainArea{i_BrArea}), 'mdT_L') ||   strcmp((Ana_TargetBrainArea{i_BrArea}), 'mdT_R')
            %                 BrainArea = repmat({[ '_', (Ana_TargetBrainArea{i_BrArea})]}, length(InVal_idx), 1);
            %             else
            BrainArea = repmat((Ana_TargetBrainArea(i_BrArea)), length(InVal_idx), 1);
            %             end
            %             Criterium_SpkPerSec = repmat(Criterium_SpkPerSec, length(InVal_idx), 1);
            %             Criterium_NrCardiacCycles = repmat(Criterium_NrCardiacCycles, length(InVal_idx), 1);
            Nr_InVal_idx2 = repmat(length(InVal_idx2), length(InVal_idx), 1);
            Nr_InVal_idx1 = repmat(length(InVal_idx1), length(InVal_idx), 1);
            NanUnits_BeforeExclusion = repmat(NanUnits_BeforeExclusion, length(InVal_idx), 1);
            ValidUnits_BeforeExclusion = repmat(ValidUnits_BeforeExclusion, length(InVal_idx), 1);
            
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).quantSNR(InVal_idx)        = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).Single_rating(InVal_idx)   = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).stability_rating(InVal_idx)= [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).unit_ID(InVal_idx)         = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).target(InVal_idx)          = [] ;
            
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD(InVal_idx, :)           = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD_SEM(InVal_idx, :)       = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDP(InVal_idx, :)          = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_period(InVal_idx, :)   = [] ;
            
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR(InVal_idx)              = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_n_bins(InVal_idx)      = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_time(InVal_idx)        = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_sign(InVal_idx)        = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_FR_diff(InVal_idx)     = [] ;
            
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).NrEvents(InVal_idx)        = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).NrTrials(InVal_idx)        = [] ;
            
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP(InVal_idx, :)            = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized(InVal_idx, :) = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_SubtrSDP(InVal_idx)           = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_PcS(InVal_idx)                = [] ;
            
            NanUnits_AfterExclusion    =  sum(isnan(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD(:,1) ));
            ValidUnits_AfterExclusion  =  sum(~isnan(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD(:,1) ));
            NanUnits_AfterExclusion    = repmat(NanUnits_AfterExclusion, length(InVal_idx), 1);
            ValidUnits_AfterExclusion  = repmat(ValidUnits_AfterExclusion, length(InVal_idx), 1);
            
            %             ExcludedUnits = table(NanUnits_BeforeExclusion,NanUnits_AfterExclusion, ValidUnits_BeforeExclusion, ValidUnits_AfterExclusion, Criterium_SpkPerSec,Nr_InVal_idx1, FR, Criterium_NrCardiacCycles, Nr_InVal_idx2, TaskType,BrainArea, InVal_unit_ID,InVal_idx,InVal_Nr  );
            %             Tab_ExcludedUnits = [Tab_ExcludedUnits;ExcludedUnits ];
        end
    end
end

for i_BrArea = 1: length(fieldnames(Out))
    disp([Ana_TargetBrainArea{i_BrArea} ' Units after exclusion criteria '])
    unitList_afterExclusion.(Ana_TargetBrainArea{i_BrArea}) = ...
        union(Out.(Ana_TargetBrainArea{i_BrArea}).Rest.unit_ID, Out.(Ana_TargetBrainArea{i_BrArea}).Task.unit_ID);
    length(unitList_afterExclusion.(Ana_TargetBrainArea{i_BrArea}))
end

disp('Saved list of units after exclusion')
save([basepath_to_save filesep 'units_after_exclusion.mat'], 'unitList_afterExclusion')

% if saveTable == 1
% save([basepath_to_save,filesep , 'Table_excludedUnits' ],'Tab_ExcludedUnits');
% filename = [basepath_to_save,filesep , 'Table_excludedUnits2.xlsx' ];
% % writetable(Tab_ExcludedUnits,filename,'Sheet',1,  'Range' ,'A1' )
% disp(['SAVED   ', basepath_to_save,filesep , 'Table_excludedUnits' ])
% end

%%  Calculations
for i_BrArea = 1: length(fieldnames(Out))
    for i_tsk = 1: numel(TaskTyp)
        out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
        
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_perECGTriggeredAverage = nanmean(out.SD,2);
        % mean(SD - SDP)
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean                    =   nanmean(out.SDsubstractedSDP);
        % standard error of the SDmean
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SEM                =  nanstd(out.SDsubstractedSDP)/ sqrt(length(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean )) ;
        
    end
end


%% overview about all units
hf = figure('Name',sprintf('BarPlot'),'Position',[1000 150 600 400],'PaperPositionMode', 'auto');

for i_BrArea = 1: length(fieldnames(Out))
    for i_tsk = 1: numel(TaskTyp)
        
        out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
        Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
        Idx_Units_NaN =  sum(~Idx_Units_NonNaN);
        idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
        
        if i_tsk == 1
            % increase, decrease, non-sign
            Pc_SignFR_rest(i_BrArea,:) = ([sum(out.sig_sign(idx_sig) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] / sum(Idx_Units_NonNaN)) *100;
            Nb_SignFR_rest(i_BrArea,:) = ([sum(out.sig_sign(idx_sig) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] ) ;
        else
            Pc_SignFR_task(i_BrArea,:) = ([sum(out.sig_sign(idx_sig) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] / sum(Idx_Units_NonNaN)) *100;
            Nb_SignFR_task(i_BrArea,:) = ([sum(out.sig_sign(idx_sig) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] ) ;
        end
        
    end
    
end

subplot(1,2,1);
barpairs =  [Pc_SignFR_rest];
b = bar(barpairs,'stacked', 'Facecolor','flat');
b(3).FaceColor = [1 1 1];
title('Rest','interpreter','none');
set(gca,'XTickLabel',fieldnames(Out),'fontsize',10);
ylim([0 100])
ylabel('Percentage of Units, %')

subplot(1,2,2);
barpairs =  [Pc_SignFR_task];
b = bar(barpairs,'stacked', 'Facecolor','flat');
b(3).FaceColor = [1 1 1];
title('Task','interpreter','none');
set(gca,'XTickLabel',fieldnames(Out),'fontsize',10);
ylim([0 100])
% legend(b, {'increase FR', 'decrease FR', 'non-significant'}, 'Location', 'Best')

filename= ['Pc_CardiacRelatedUnits'];

if savePlot;
    save([basepath_to_save,filesep ,filename '.mat'], 'Pc_SignFR_rest', 'Pc_SignFR_task')
    export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
    close all;
end







%% How much does the surrogate and mean firing divergate
if savePlot == 1
    hf = figure('Name',sprintf('Surrogate_MeanFr'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for i_BrArea = 1: length(fieldnames(Out))
        for i_tsk = 1: numel(TaskTyp)
            out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            
            ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);
            hold on;
            if numel(out.FR) > 1 && sum(~isnan(out.FR))
                scatter(out.FR, nanmean(out.SDP,2), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})  %,30, out.FR(idx_sig)/max(out.FR) , 'filled')
                set(gca,'ylim',[0,max(out.FR)]);
                set(gca,'xlim',[0,max(out.FR)]);
                x = linspace(0,max(out.FR)); y = linspace(0,max(out.FR)); hold on;
                plot(x,y, 'k');
                xlabel('average Firing rate','fontsize',14 );
                ylabel('mean FR of jittered data','fontsize',14 );
                axis square; box on;
                title([ (TaskTyp{i_tsk}), Ana_TargetBrainArea{i_BrArea}]);
            end
            
            ha1 = subplot(2,length(fieldnames(Out)), (i_BrArea + length(fieldnames(Out))));
            hold on;
            if size(out.SD, 1) > 1  && sum(~isnan(nanmean(out.SD,2)))
                scatter(nanmean(out.SD,2) , nanmean(out.SDP,2), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})  %,30, out.FR(idx_sig)/max(out.FR) , 'filled')
                set(gca,'ylim',[0,max(nanmean(out.SD,2))]);
                set(gca,'xlim',[0,max(nanmean(out.SD,2))]);
                x = linspace(0,max(nanmean(out.SD,2))); y = linspace(0,max(nanmean(out.SD,2))); hold on;
                plot(x,y, 'k');
                xlabel('average Firing rate within analysis window','fontsize',14 );
                ylabel('mean FR of jittered data','fontsize',14 );
                axis square;
                box on;
                title([  Ana_TargetBrainArea{i_BrArea}]);
            end
        end
    end
    
    filename= ['MeanSurrogate_MeanFr'];
    
    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
    
    
    %% Example for
    if true
        for i_BrArea = 1:length(fieldnames(Out)); %%% 3
            for i_tsk = 1: numel(TaskTyp)
                hf = figure('Name',sprintf(Ana_TargetBrainArea{i_BrArea}),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
                
                out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
                idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
                
                ha1 = subplot(2,4,[1]);
                box on;
                hold on;
                if ~isempty(out.quantSNR)
                    scatter(out.quantSNR(idx_sig) , out.FR_ModIndex_PcS(idx_sig), 'filled', 'MarkerFaceColor', condition_colors{i_tsk});
                    scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_PcS(~idx_sig), 'filled', 'MarkerFaceColor', condition_colors{i_tsk}/2);
                    xlabel('Signal-to-Noise','fontsize',14 );
                    ylabel('Modulation index (%)','fontsize',14 );
                    axis square;
                
                    xf = [min(out.quantSNR), max(out.quantSNR)];
                    [p,S] = polyfit(out.quantSNR(~isnan(out.FR_ModIndex_PcS)),out.FR_ModIndex_PcS(~isnan(out.FR_ModIndex_PcS)),1); %
                    [y_fit,delta] = polyval(p,out.quantSNR(~isnan(out.FR_ModIndex_PcS)),S);
                    [coef, pval] = corr(out.quantSNR,out.FR_ModIndex_PcS, 'rows','complete') ;
                    plot(out.quantSNR(~isnan(out.FR_ModIndex_PcS)), y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
                    title({[Ana_TargetBrainArea{i_BrArea} , ' ' ,(TaskTyp{i_tsk}) ], ...
                        ['coef, p ', num2str([roundto(coef,2), roundto(pval,4)])]})
                    legend({'Significant Units', 'Non-Significant Units', 'Linear Fit (Overall)'}, 'Location', 'Best')
                end
                
                ha1 = subplot(2,4,[2]);
                box on;
                hold on;
                scatter(out.FR(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor', condition_colors{i_tsk})
                scatter(out.FR(~idx_sig) , out.FR_ModIndex_SubtrSDP(~idx_sig), 'filled', 'MarkerFaceColor', condition_colors{i_tsk}/2)
                ylabel('Modulation index (%)','fontsize',14 );
                xlabel('mean FR','fontsize',14 );
                axis square;
                hold on;
                legend({'Significant Units', 'Non-Significant Units'}, 'Location', 'Best')
                
                ha1 = subplot(2,4,[3]);
                box on;
                hold on;
                scatter(out.sig_n_bins(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor', condition_colors{i_tsk})
                scatter(out.sig_n_bins(~idx_sig) , out.FR_ModIndex_SubtrSDP(~idx_sig), 'filled', 'MarkerFaceColor', condition_colors{i_tsk}/2)
                ylabel('Modulation index (%)','fontsize',14 );
                xlabel('sig. Nr. bins','fontsize',14 );
                axis square;
                hold off;
                legend({'Significant Units', 'Non-Significant Units'}, 'Location', 'Best')
                
                ha1 = subplot(2,4,[5:6]);
                box on
                UnitSig_Rest = out.FR_ModIndex_SubtrSDP >30 & idx_sig;
                UnitNotSign_Rest = out.FR_ModIndex_PcS > 30 & ~idx_sig;
                
                text(-400,-15, [out.unit_ID(UnitSig_Rest)],'Color', condition_colors{i_tsk});
                
                if sum(UnitSig_Rest)
                    line((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000 , out.SDsubstractedSDP_normalized(UnitSig_Rest,:), 'color', condition_colors{i_tsk}, 'LineWidth', 1);
                    hold on;
                end
                hold on;
                text(300,20, [out.unit_ID(UnitNotSign_Rest)],'Color','k');
                xlabel('Time from R-peak, ms')
                ylabel('% signal change')
                title('Units with Modulation strength > 30%'); %axis square;
                
                ha1 = subplot(2,4,[7:8]);
                box on
                text(300,20, [out.unit_ID(UnitNotSign_Rest)],'Color','k');
                if sum(UnitNotSign_Rest)
                    line((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000 , out.SDsubstractedSDP_normalized(UnitNotSign_Rest,:), 'color', condition_colors{i_tsk},'LineWidth', 4);
                end
                hold on;
                xlabel('Time from R-peak, ms')
                ylabel('% signal change')
                title('Non-Sign units with Modulation strength > 30%'); hold off;
                if savePlot;
                    export_fig([basepath_to_save,filesep ,['Check_ModulationIndex_',Ana_TargetBrainArea{i_BrArea},'_',(TaskTyp{i_tsk}) ]], '-pdf'); %,'-transparent'
                    close all;
                end
            end
        end
    end
    
    if Graph_SelectionCriterion
        %%  normalizing the Firing rate calcuations
        for i_BrArea = 1: length(fieldnames(Out))
            hf = figure('Name',sprintf(Ana_TargetBrainArea{i_BrArea}),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
            for i_tsk = 1: numel(TaskTyp)
                
                out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
                
                idx_sig =  ~isnan(out.sig_FR_diff);
                
                ha1 = subplot(2,2,[1:2]); %
                hold on
                box on
                % bar plot how many are significant & positive and negative FR?
                SDmean_SEM = nanstd(out.SDsubstractedSDP(idx_sig, :))/ sqrt(length(nanmean(out.SDsubstractedSDP(idx_sig, :)))) ;
                shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,nanmean(out.SDsubstractedSDP(idx_sig, :)) ,SDmean_SEM, lineProps{i_tsk}, 1);
                ylabel('Signal Change, spikes/s')
                title(['Pop:  (all significant) Cal Per Unit:SD-SDP ' (Ana_TargetBrainArea{i_BrArea}) ' units'],'interpreter','none');
                legend(TaskTyp)
                
                ha1 = subplot(2,2,[3:4]); %
                hold on
                box on
                SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_sig, :),[], 1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_sig, :), 1))) ;
                shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_sig, :), 1) ,SDmean_SEM ,lineProps{i_tsk},1);
                title(['Population:  (all significant)' (Ana_TargetBrainArea{i_BrArea}) ' units'],'interpreter','none');
                ylabel('% signal change','fontsize',14 );
                xlabel('Time relative to R-peak (ms)','fontsize',14 );
                title(['Pop:  (all significant) Cal Per Unit:(SD-SDP)/SDP ' (Ana_TargetBrainArea{i_BrArea}) ' units'],'interpreter','none');
                legend(TaskTyp)
            end
            if savePlot
                filename= ['Normalization_FR_ECG_triggered_spike_' ,(Ana_TargetBrainArea{i_BrArea})];
                export_fig([basepath_to_save filesep filename ], '-pdf'); %
                close all;
            end
        end
        
        %% Plot the dataset - check up
        for i_BrArea = 1: length(fieldnames(Out))
            for i_tsk = 1: numel(TaskTyp)
                out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
                
                hf = figure('Name',sprintf([(Ana_TargetBrainArea{i_BrArea}) '  ' TaskTyp{i_tsk}]),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
                rows_plot = 2;
                colums_plot = 7;
                
                for I = 1: 5
                    % How many different groups
                    if I == 1
                        ha1 = subplot(rows_plot,colums_plot,[4:5]); %
                        
                        idx1 = ([out.sig_n_bins] > 0 );
                        idx2 = ([out.FR_ModIndex_SubtrSDP] > 5 );
                        
                        Y1 = out.SD( idx1 & idx2,:)  ;
                        Y2 = out.SDP( idx1 & idx2,:)  ;
                        
                        if sum(idx2) ~= size(Y1,1)
                            disp('incorrect Selection')
                        end
                        title(['FR_ModIndex_SubtrSDP > 5 '],'interpreter','none');
                        
                    elseif I == 2
                        ha1 = subplot(rows_plot,colums_plot,[6:7]); %
                        idx1 = ([out.sig_n_bins] > 0 );
                        idx_ex = ([out.FR_ModIndex_SubtrSDP] < 5 );  %Smaller
                        Y1 = out.SD( idx1 & idx_ex,:)  ;
                        Y2 = out.SDP( idx1 & idx_ex,:)  ;
                        title(['FR_ModIndex_SubtrSDP < 5 '],'interpreter','none');      %& Nr. of bin: 0 - 10 bins for sign. Intervals
                        
                    elseif I == 3
                        ha1 = subplot(rows_plot,colums_plot,[8:9]); %
                        
                        idx1 = ([out.sig_n_bins] > 0 );
                        idx2 = ([out.sig_n_bins] <= 4 );
                        
                        Y1 = out.SD( idx1 & idx2,:)  ;
                        Y2 = out.SDP( idx1 & idx2,:)  ;
                        
                        %Y1 = out.SD( out.CmbDataSig_SD(:, end-2) < 5 & out.CmbDataSig_SD(:, end-2) > 0,:);
                        %Y2 = out.SDP( out.CmbDataSig_SD(:, end-2) < 5 & out.CmbDataSig_SD(:, end-2) > 0,:);
                        title(['Nr. of bin: 1 - 4 bins for sign. Intervals'],'interpreter','none');
                        
                    elseif I == 4
                        idx1 = ([out.sig_n_bins] >= 5  );
                        idx2 = ([out.sig_n_bins] < 7 );
                        
                        Y1 = out.SD( idx1 & idx2,:)  ;
                        Y2 = out.SDP( idx1 & idx2,:)  ;
                        
                        
                        ha2 = subplot(rows_plot,colums_plot,[10:11]); %
                        title(['Nr. of bin: 5 - 7 bins for sign. Intervals'],'interpreter','none');
                        
                    elseif I == 5
                        
                        idx1 = ([out.sig_n_bins] > 6 );
                        
                        Y1 = out.SD( idx1 ,:)  ;
                        Y2 = out.SDP( idx1 ,:)  ;
                        ha3 = subplot(rows_plot,colums_plot,[12:13]); %
                        title(['Nr. of bin: above 6 bins for sign. Intervals'],'interpreter','none');
                        
                    end
                    
                    if ~isempty(Y1)
                        cmap=jet(size(Y1,1));
                        colororder(cmap)
                        line((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000 , Y1 - Y2, 'linewidth', 3)
                        set(gca,'ylim',[-9, 9]);
                        ylabel('Firing rate (%)','fontsize',12 );
                        xlabel('Time relative to ECG peak (ms)','fontsize',12);
                    end
                    hold off;
                end
                
                idx_ex = ([out.sig_n_bins] <= 4 );
                
                ha4 = subplot(rows_plot,colums_plot,[3]); %
                scatter(out.sig_n_bins , out.FR_ModIndex_SubtrSDP, 'filled', 'MarkerFaceColor', condition_colors{i_tsk}); hold on;
                scatter(out.sig_n_bins(idx_ex) , out.FR_ModIndex_SubtrSDP(idx_ex), 'filled', 'MarkerFaceColor',[0 0 0])
                ylabel('FR Modulation (spike/s)','fontsize',14 );
                xlabel('Nr. of bins of sig. Interval','fontsize',14 );
                hold off;  axis square;
                
                ha4 = subplot(rows_plot,colums_plot,[1:2]); %
                scatter(out.sig_n_bins , out.quantSNR, 'filled', 'MarkerFaceColor', condition_colors{i_tsk}); hold on;
                scatter(out.sig_n_bins(idx_ex) , out.quantSNR(idx_ex), 'filled', 'MarkerFaceColor',[0 0 0])
                ylabel('Signal-to-Noise Ratio','fontsize',14 );
                xlabel('Nr. of bins of sig. Interval','fontsize',14 );
                hold off;  axis square;
                
                ha4 = subplot(rows_plot,colums_plot,[14]); %
                scatter(out.sig_time , out.FR_ModIndex_SubtrSDP, 'filled', 'MarkerFaceColor', condition_colors{i_tsk}); hold on;
                scatter(out.sig_time(idx_ex) , out.FR_ModIndex_SubtrSDP(idx_ex), 'filled', 'MarkerFaceColor',[0 0 0])
                ylabel('FR Modulation (spike/s)','fontsize',10);
                xlabel('TimePoint of sig. highest diff in FR','fontsize',10 );
                hold off; axis square;
                
                
                if savePlot
                    filename= ['ModulationInRelationToBins_ECG_triggered_spike_' (TaskTyp{i_tsk}) , '_' ,(Ana_TargetBrainArea{i_BrArea})];
                    export_fig([basepath_to_save filesep filename ], '-pdf'); %
                    close all;
                end
            end
        end
    end
    
    %% Plot the averages
    for i_BrArea = 1: length(fieldnames(Out))
        hf = figure('Name',sprintf(Ana_TargetBrainArea{i_BrArea}),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
        for i_tsk = 1: numel(TaskTyp)
            
            out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            % from all nicht NAN units - how many were significant?
            Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
            Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
            
            idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
            %idx_sig = double(idx_sig);
            
            
            hold on
            ha1 = subplot(2,4,[1:2]);
            % bar plot how many are significant & positive and negative FR?
            if sum(idx_sig)
                line((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000 , out.SDsubstractedSDP_normalized(idx_sig,:), 'color', condition_colors{i_tsk}, 'LineWidth', 1);
                hold on;
            end
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_sig, :), [], 1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_sig, :), 1))) ;
            shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_sig, :), 1) ,SDmean_SEM,lineProps{i_tsk},1);
            title({[TaskTyp{i_tsk} ': units = ' ,num2str(sum(idx_sig)), ' of ' ,num2str(sum(Idx_Units_NonNaN)) ], ...
                ['Population:  (all significant)' (Ana_TargetBrainArea{i_BrArea}) ' units']},'interpreter','none');
            ylabel('normalized Firing rate (%)','fontsize',14 );
            xlabel('Time relative to R-peak (ms)','fontsize',14 );
            hold off;
            
            %Bar graph
            ha1 = subplot(2,4,3);% pie plot how many
            
            if i_tsk == 1
                Pc_SignFR_rest = ([sum(out.sig_sign(idx_sig) == 1  ), sum(out.sig_sign(idx_sig) == -1), (sum(~idx_sig) -Idx_Units_NaN)  ] / sum(Idx_Units_NonNaN)) *100;
                Pc_SignFR_task = ([sum(out.sig_sign(idx_sig) == 1  ), sum(out.sig_sign(idx_sig) == -1), (sum(~idx_sig) -Idx_Units_NaN)  ] / sum(Idx_Units_NonNaN)) *100;
                
            else
                Pc_SignFR_task = ([sum(out.sig_sign(idx_sig) == 1  ), sum(out.sig_sign(idx_sig) == -1), (sum(~idx_sig) -Idx_Units_NaN)  ] / sum(Idx_Units_NonNaN)) *100;
            end
            barpairs =  [Pc_SignFR_rest; Pc_SignFR_task];
            b = bar(barpairs,'stacked', 'Facecolor','flat' );
            b(3).FaceColor = [1 1 1];
            legend(b, {'increase FR', 'decrease FR', 'non-significant'}, 'Location', 'Best')
            set(gca,'XTickLabel',TaskTyp,'fontsize',10);
            
            % display only significant units showing a increase in FR
            ha1 = subplot(2,4,[5]); %
            idx_SigDec = (out.sig_sign == -1);
            idx_SigInc = (out.sig_sign == 1);
            text(-400,1, [TaskTyp{i_tsk} ': units = ' ,num2str(sum(idx_SigInc & idx_sig)) ],'Color',condition_colors{i_tsk})
            
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig,:), 1))) ;
            shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig,:),1), SDmean_SEM ,lineProps{i_tsk},1);
            set(gca,'ylim',[-10, 10]);
            title(['units showing a sig. INCREASE in FR'],'interpreter','none');
            ylabel('normalized Firing rate (%)','fontsize',14 );
            xlabel('Time relative to R-peak (ms)','fontsize',14 );
            % display only significant units showing a decrease in FR
            ha1 = subplot(2,4,[6]);
            text(-400,1, [TaskTyp{i_tsk} ': units = ' ,num2str(sum(idx_SigInc & idx_sig)) ],'Color',condition_colors{i_tsk})
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig,:), 1))) ;
            shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig,:),1), SDmean_SEM ,lineProps{i_tsk},1);
            set(gca,'ylim',[-10, 10]);
            title(['units showing a sig. DECREASE in FR'],'interpreter','none');
            ylabel('normalized Firing rate (spike/s)','fontsize',14 );
            xlabel('Time relative to R-peak (ms)','fontsize',14 );
            
            
            %% MODULATION INDEX
            
            ha1 = subplot(2,4,[4]);
            hold on;
            if ~isempty(out.quantSNR)
                scatter(out.quantSNR(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
                ylabel('Modulation index (%)','fontsize',14 );
                xlabel('Signal-to-Noise','fontsize',14 );
                axis square;
                hold on;
            
                xf = [min(out.quantSNR), max(out.quantSNR)];
                [p,S] = polyfit(out.quantSNR,out.FR_ModIndex_PcS,1); %
                [y_fit,delta] = polyval(p,out.quantSNR,S);
                [coef, pval] = corr(out.quantSNR,out.FR_ModIndex_PcS, 'rows','complete') ;
                scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_PcS(~idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
                plot(out.quantSNR, y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
                text(8,35, ['coef, p ', num2str([roundto(coef,2),roundto(pval,3)])], 'Color',condition_colors{i_tsk})
            end
            
            if i_tsk == 1
                % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
            else
                % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
            end
            
            
            %%
            Dat = [];
            ModIndex = out.FR_ModIndex_SubtrSDP(idx_sig);
            SNR =  out.quantSNR(idx_sig);
            %Dat = table(ModIndex, SNR);
            
            ha1 = subplot(2,4,[7]);        hold on;
            
            if sum(idx_sig)
                
                scatter(SNR, ModIndex, 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
                ylabel('Modulation index(%)','fontsize',14 );
                xlabel('Signal-to-Noise','fontsize',14 );        axis square;
                hold on;
                
                xf = [min(out.quantSNR(idx_sig)), max(out.quantSNR(idx_sig))];
                [p,S] = polyfit(SNR,ModIndex,1); %
                [y_fit,delta] = polyval(p,SNR,S);
                [coef, pval] = corr(SNR,ModIndex) ;
                plot(SNR, y_fit,'LineWidth', 2, 'Color',condition_colors{i_tsk});
                text(8,35, ['coef, p ', num2str([roundto(coef,2),roundto(pval,3)])], 'Color', condition_colors{i_tsk})
                if i_tsk == 1
                    % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
                else
                    % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
                end
                
            end
            hold off;
            
            ha1 = subplot(2,4,[8]); %
            scatter(out.FR(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
            
            %
            %         xf = [min(out.quantSNR(idx_sig)), max(out.quantSNR(idx_sig))];
            %         [p,S] = polyfit(SNR,ModIndex,1); %
            %         [y_fit,delta] = polyval(p,SNR,S);
            %         [coef, pval] = corr(SNR,ModIndex) ;
            %         if i_tsk == 1
            %
            %             plot(SNR, y_fit,'LineWidth', 2, 'Color', 'b');
            %             text(8,35, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', 'b')
            %             % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
            %         else
            %             plot(SNR, y_fit,'LineWidth', 2, 'Color', 'r');
            %             text(8,30, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', 'r')
            %             % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
            %         end
            
            ylabel('Modulation index(%)','fontsize',14 );
            xlabel('mean firing rate','fontsize',14);
            title(['all significant units'],'interpreter','none');
            
            axis square;
            
        end
        filename= ['Average_ECG_triggered_spike_' (Ana_TargetBrainArea{i_BrArea})];
        
        if savePlot;
            export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
            close all;
        end
        
    end
    
    
    %% MODULATION STRENGTH - for the percentage signal change
    hf = figure('Name',sprintf('ModulationIndex'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for i_BrArea = 1: length(fieldnames(Out))
        for i_tsk = 1: numel(TaskTyp)
            
            out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            % from all nicht NAN units - how many were significant?
            Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
            Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
            
            idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;

            %% Modulation Index
            ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);
            hold on;
            if ~isempty(out.quantSNR)
                % all not significant & NaN units
                scatter(out.quantSNR(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
                ylabel('Modulation index (%)','fontsize',14 );
                xlabel('Signal-to-Noise','fontsize',14 );        axis square; box on;
                title(Ana_TargetBrainArea{i_BrArea})
                set(gca,'xlim',[0, 25]);
                set(gca,'ylim',[0, 80]);
                hold on;
            
                xf = [min(out.quantSNR), max(out.quantSNR)];
                [p,S] = polyfit(out.quantSNR(~isnan(out.FR_ModIndex_PcS)),out.FR_ModIndex_PcS(~isnan(out.FR_ModIndex_PcS)),1); %
                [y_fit,delta] = polyval(p,out.quantSNR(~isnan(out.FR_ModIndex_PcS)),S);
                [coef, pval] = corr(out.quantSNR,out.FR_ModIndex_PcS, 'rows','complete') ;
                scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_PcS(~idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk}/2)
                plot(out.quantSNR(~isnan(out.FR_ModIndex_PcS)), y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
                text(8,35, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', condition_colors{i_tsk})
            end
            if i_tsk == 1
                % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
            else
                % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
            end
            
            Dat = [];
            ModIndex = out.FR_ModIndex_SubtrSDP(idx_sig);
            SNR =  out.quantSNR(idx_sig);
            %Dat = table(ModIndex, SNR);
            
            ha1 = subplot(2,length(fieldnames(Out)),i_BrArea +length(fieldnames(Out)));        hold on;
            
            scatter(SNR, ModIndex, 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
            ylabel('Modulation index(%)','fontsize',14 );
            xlabel('Signal-to-Noise','fontsize',14 );        axis square;box on;
            title(Ana_TargetBrainArea{i_BrArea});
            set(gca,'xlim',[0, 25]);
            set(gca,'ylim',[0, 80]);
            
            hold on;
            
            if sum(idx_sig)
                
                xf = [min(out.quantSNR(idx_sig)), max(out.quantSNR(idx_sig))];
                [p,S] = polyfit(SNR,ModIndex,1); %
                [y_fit,delta] = polyval(p,SNR,S);
                [coef, pval] = corr(SNR,ModIndex) ;
                plot(SNR, y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
                text(8,35, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', condition_colors{i_tsk})
                if i_tsk == 1
                    % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
                else
                    % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
                end
            end
        end
    end
    filename= ['ModulationIndex_Npsc'];
    
    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
    
    %% MODULATION STRENGTH - for SDP subtraction
    hf = figure('Name',sprintf('ModulationIndex'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for i_BrArea = 1: length(fieldnames(Out))
        for i_tsk = 1: numel(TaskTyp)
            
            out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            % from all nicht NAN units - how many were significant?
            Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
            Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
            
            idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
            
            %% Modulation Index
            ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);
            hold on;
            % all not significant & NaN units
            if ~isempty(out.quantSNR)
                scatter(out.quantSNR(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
                ylabel('Modulation index (spike/s)','fontsize',14 );
                xlabel('Signal-to-Noise ratio (mV)','fontsize',14 );
                axis square; box on;
                title(Ana_TargetBrainArea{i_BrArea})
                set(gca,'xlim',[0, 25]);
                set(gca,'ylim',[0, 15]);
                hold on;
                
                xf = [min(out.quantSNR), max(out.quantSNR)];
                [p,S] = polyfit(out.quantSNR(~isnan(out.FR_ModIndex_SubtrSDP)),out.FR_ModIndex_SubtrSDP(~isnan(out.FR_ModIndex_SubtrSDP)),1); %
                [y_fit,delta] = polyval(p,out.quantSNR(~isnan(out.FR_ModIndex_SubtrSDP)),S);
                [coef, pval] = corr(out.quantSNR,out.FR_ModIndex_SubtrSDP, 'rows','complete') ;
                scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_SubtrSDP(~idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk}/2)
                plot(out.quantSNR(~isnan(out.FR_ModIndex_SubtrSDP)), y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
                text(8,10, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', condition_colors{i_tsk})
            end
            if i_tsk == 1
                % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
            else
                % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
            end
            
            ha1 = subplot(2,length(fieldnames(Out)),i_BrArea +length(fieldnames(Out)));
            hold on;
            
            scatter(out.quantSNR(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
            ylabel('Modulation index(spike/s)','fontsize',14 );
            xlabel('Signal-to-Noise ratio (mV)','fontsize',14 );
            axis square;box on;
            title(Ana_TargetBrainArea{i_BrArea});
            set(gca,'xlim',[0, 25]);
            set(gca,'ylim',[0, 15]);
            hold on;
            
            if sum(idx_sig)
                
                xf = [min(out.quantSNR(idx_sig)), max(out.quantSNR(idx_sig))];
                [p,S] = polyfit(out.quantSNR(idx_sig),out.FR_ModIndex_SubtrSDP(idx_sig),1); %
                [y_fit,delta] = polyval(p,out.quantSNR(idx_sig),S);
                
                [coef, pval] = corr(out.quantSNR(idx_sig),out.FR_ModIndex_SubtrSDP(idx_sig)) ;
                plot(out.quantSNR(idx_sig), y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
                text(8,10, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', condition_colors{i_tsk},'fontsize',6)
                if i_tsk == 1
                    % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
                else
                    % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
                end
            end
        end
    end
    filename= ['ModulationIndex_NSubtr'];
    
    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
    
    %% MODULATION STRENGTH - compare both Normalizations ... SDP subtraction
    hf = figure('Name',sprintf('ModulationIndex'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for i_BrArea = 1: length(fieldnames(Out))
        for i_tsk = 1: numel(TaskTyp)
            
            out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            % from all nicht NAN units - how many were significant?
            Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
            Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
            
            idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
            
            %% Modulation Index
            ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);
            hold on;
            % all not significant & NaN units
            if ~isempty(out.quantSNR)
            scatter(out.FR_ModIndex_PcS(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
            ylabel('Modulation index (spike/s)','fontsize',14 );
            xlabel('Modulation index(%pSc)','fontsize',14 );
            axis square; box on;
            title(Ana_TargetBrainArea{i_BrArea})
            set(gca,'xlim',[0, 80]);
            set(gca,'ylim',[0, 15]);
            hold on;
            
            xf = [min(out.quantSNR), max(out.quantSNR)];
            [p,S] = polyfit(out.FR_ModIndex_PcS(~isnan(out.FR_ModIndex_SubtrSDP)),out.FR_ModIndex_SubtrSDP(~isnan(out.FR_ModIndex_SubtrSDP)),1); %
            [y_fit,delta] = polyval(p,out.FR_ModIndex_PcS(~isnan(out.FR_ModIndex_SubtrSDP)),S);
            [coef, pval] = corr(out.FR_ModIndex_PcS,out.FR_ModIndex_SubtrSDP, 'rows','complete') ;
            scatter(out.FR_ModIndex_PcS(~idx_sig) , out.FR_ModIndex_SubtrSDP(~idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk}/2)
            plot(out.FR_ModIndex_PcS(~isnan(out.FR_ModIndex_SubtrSDP)), y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
            text(8,10, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', condition_colors{i_tsk},'fontsize',14)
            end
            if i_tsk == 1
                % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
            else
                % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
            end
            
            ha1 = subplot(2,length(fieldnames(Out)),i_BrArea +length(fieldnames(Out)));
            hold on;
            scatter(out.FR_ModIndex_PcS(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
            ylabel('Modulation index(spike/s)','fontsize',14 );
            xlabel('Modulation index(%pSc)','fontsize',14 );        axis square;box on;
            title(Ana_TargetBrainArea{i_BrArea});
            set(gca,'xlim',[0, 80]);
            
            set(gca,'ylim',[0, 15]);
            
            
            hold on;
            
            if sum(idx_sig)
                xf = [min(out.FR_ModIndex_PcS(idx_sig)), max(out.quantSNR(idx_sig))];
                [p,S] = polyfit(out.FR_ModIndex_PcS(idx_sig),out.FR_ModIndex_SubtrSDP(idx_sig),1); %
                [y_fit,delta] = polyval(p,out.FR_ModIndex_PcS(idx_sig),S);
                
                [coef, pval] = corr(out.FR_ModIndex_PcS(idx_sig),out.FR_ModIndex_SubtrSDP(idx_sig)) ;
                plot(out.FR_ModIndex_PcS(idx_sig), y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
                text(8,10, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', condition_colors{i_tsk},'fontsize',14)
                if i_tsk == 1
                    % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
                else
                    % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
                end
            end
        end
    end
    filename= ['ModulationIndex_Cmp_Npsc_NSubtr'];
    
    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
    
    %% MODULATION STRENGTH - FR as explanation - compare both Normalizations ... SDP subtraction
    hf = figure('Name',sprintf('ModulationIndex'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for i_BrArea = 1: length(fieldnames(Out))
        for i_tsk = 1: numel(TaskTyp)
            
            out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            % from all nicht NAN units - how many were significant?
            Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
            Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
            
            idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
            
            %% Modulation Index
            if i_tsk == 1
                ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);        hold on;
                
                
                % all not significant & NaN units
                if ~isequal(idx_sig, zeros(length(idx_sig),1))
                    % scatter(out.FR_ModIndex_PcS(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',Color)
                    scatter(out.FR_ModIndex_PcS(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig),30, out.FR_perECGTriggeredAverage(idx_sig)/max(out.FR_perECGTriggeredAverage(idx_sig)) , 'filled')
                    %scatter(out.FR_ModIndex_PcS(~idx_sig) , out.FR_ModIndex_SubtrSDP(~idx_sig),30, out.FR(~idx_sig)/max(out.FR) , 'filled')
                    
                    colorbar();
                    colormap jet
                    ylabel('Modulation index (spike/s)','fontsize',14 );
                    xlabel('Modulation index(%pSc)','fontsize',14 );        axis square; box on;
                    title([ (TaskTyp{i_tsk}), Ana_TargetBrainArea{i_BrArea}]);
                    if strcmp(Ana_TargetBrainArea{i_BrArea}, 'VPL_R') || strcmp(Ana_TargetBrainArea{i_BrArea}, 'dPul_R') || strcmp(Ana_TargetBrainArea{i_BrArea}, 'MD')
                        set(gca,'xlim',[0, 60]);
                    end
                    set(gca,'ylim',[0, 15]);
                end
            end
            if i_tsk == 2
                ha1 = subplot(2,length(fieldnames(Out)),(i_BrArea + length(fieldnames(Out))));        hold on;
                if ~isequal(idx_sig, zeros(length(idx_sig),1))
                    % all not significant & NaN units
                    % scatter(out.FR_ModIndex_PcS(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',Color)
                    scatter(out.FR_ModIndex_PcS(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig),30, out.FR_perECGTriggeredAverage(idx_sig)/max(out.FR_perECGTriggeredAverage(idx_sig)) , 'filled')
                    % scatter(out.FR_ModIndex_PcS(~idx_sig) , out.FR_ModIndex_SubtrSDP(~idx_sig),30, out.FR(~idx_sig)/max(out.FR) , 'filled')
                    
                    colorbar
                    colormap jet
                    ylabel('Modulation index (spike/s)','fontsize',14 );
                    xlabel('Modulation index(%pSc)','fontsize',14 );        axis square; box on;
                    title([ (TaskTyp{i_tsk}),Ana_TargetBrainArea{i_BrArea}]);
                    % set(gca,'xlim',[0, 25]);
                    set(gca,'ylim',[0, 15]);
                    colorbar();
                    colormap jet
                    if strcmp(Ana_TargetBrainArea{i_BrArea}, 'VPL_R') || strcmp(Ana_TargetBrainArea{i_BrArea}, 'dPul_R') || strcmp(Ana_TargetBrainArea{i_BrArea}, 'MD')
                        set(gca,'xlim',[0, 60]);
                    end
                end
            end
        end
    end
    filename= ['ModulationIndex_Compare_NpSc_Subtr_WithFR_ONLYSignificantUnits'];
    
    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
    
    %% MODULATION STRENGTH - FR as explanation - compare both Normalizations ... SDP subtraction
    hf = figure('Name',sprintf('ModulationIndex'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for i_BrArea = 1: length(fieldnames(Out))
        for i_tsk = 1: numel(TaskTyp)
            
            out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            % from all nicht NAN units - how many were significant?
            Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
            Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
            
            idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
            
            %% Modulation Index
            if i_tsk == 1
                ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);        hold on;
                
                if ~isequal(idx_sig, zeros(length(idx_sig),1))
                    % all not significant & NaN units
                    % scatter(out.FR_ModIndex_PcS(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',Color)
                    scatter(out.quantSNR(idx_sig),out.FR_ModIndex_PcS(idx_sig) , 30, out.FR_perECGTriggeredAverage(idx_sig)/max(out.FR_perECGTriggeredAverage(idx_sig)) , 'filled')
                    %scatter(out.FR_ModIndex_PcS(~idx_sig) , out.FR_ModIndex_SubtrSDP(~idx_sig),30, out.FR(~idx_sig)/max(out.FR) , 'filled')
                    
                    colorbar();
                    colormap jet
                    ylabel('Modulation index (%sc)','fontsize',14 );
                    xlabel('signal-to-noise ratio','fontsize',14 );        axis square; box on;
                    title([ (TaskTyp{i_tsk}), Ana_TargetBrainArea{i_BrArea}]);
                    set(gca,'ylim',[0, 80]);
                    set(gca,'xlim',[0, 25]);
                end
            end
            if i_tsk == 2
                ha1 = subplot(2,length(fieldnames(Out)),(i_BrArea + length(fieldnames(Out))));        hold on;
                if ~isequal(idx_sig, zeros(length(idx_sig),1))
                    % all not significant & NaN units
                    % scatter(out.FR_ModIndex_PcS(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',Color)
                    scatter(out.quantSNR(idx_sig),out.FR_ModIndex_PcS(idx_sig) , 30, out.FR_perECGTriggeredAverage(idx_sig)/max(out.FR_perECGTriggeredAverage(idx_sig)) , 'filled')
                    % scatter(out.FR_ModIndex_PcS(~idx_sig) , out.FR_ModIndex_SubtrSDP(~idx_sig),30, out.FR(~idx_sig)/max(out.FR) , 'filled')
                    
                    colorbar
                    colormap jet
                    ylabel('Modulation index (%sc)','fontsize',14 );
                    xlabel('signal-to-noise ratio','fontsize',14 );        axis square; box on;
                    title([ (TaskTyp{i_tsk}),Ana_TargetBrainArea{i_BrArea}]);
                    colorbar();
                    colormap jet
                    set(gca,'ylim',[0, 80]);
                    set(gca,'xlim',[0, 25]);
                end
            end
        end
    end
    filename= ['ModulationIndex_NpSc_SNR_FR_ONLYSignificantUnits'];
    
    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
    %% MODULATION STRENGTH pSC- FR
    hf = figure('Name',sprintf('ModulationIndex_FiringRate'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for i_BrArea = 1: length(fieldnames(Out))
        for i_tsk = 1: numel(TaskTyp)
            
            out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            % from all nicht NAN units - how many were significant?
            Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
            Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
            
            idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
            
            %% Modulation Index
            ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);        hold on;
            % all not significant & NaN units
            if ~isempty(out.quantSNR)
            scatter(out.FR(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
            ylabel('Modulation index (%)','fontsize',14 );
            xlabel('average Firing rate','fontsize',14 );        axis square; box on;
            title(Ana_TargetBrainArea{i_BrArea})
            set(gca,'xlim',[0, 200]);
            set(gca,'ylim',[0, 80]);
            hold on;
            
            xf = [min(out.FR), max(out.quantSNR)];
            [p,S] = polyfit(out.FR(~isnan(out.FR_ModIndex_PcS)),out.FR_ModIndex_PcS(~isnan(out.FR_ModIndex_PcS)),1); %
            [y_fit,delta] = polyval(p,out.FR(~isnan(out.FR_ModIndex_PcS)),S);
            [coef, pval] = corr(out.FR,out.FR_ModIndex_PcS, 'rows','complete') ;
            scatter(out.FR(~idx_sig) , out.FR_ModIndex_PcS(~idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk}/2)
            plot(out.FR(~isnan(out.FR_ModIndex_PcS)), y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
            text(8,35, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', condition_colors{i_tsk})
            
            if i_tsk == 1
                % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
            else
                % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
            end
            end
            ha1 = subplot(2,length(fieldnames(Out)),i_BrArea +length(fieldnames(Out)));        hold on;
            
            scatter(out.FR(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
            ylabel('Modulation index(%)','fontsize',14 );
            xlabel('average Firing rate','fontsize',14 );
            axis square;box on;
            title(Ana_TargetBrainArea{i_BrArea});
            set(gca,'xlim',[0, 200]);
            set(gca,'ylim',[0, 80]);
            
            hold on;
            
            if sum(idx_sig)
                
                xf = [min(out.FR(idx_sig)), max(out.quantSNR(idx_sig))];
                [p,S] = polyfit(out.FR(idx_sig),out.FR_ModIndex_SubtrSDP(idx_sig),1); %
                [y_fit,delta] = polyval(p,out.FR(idx_sig),S);
                [coef, pval] = corr(out.FR(idx_sig),out.FR_ModIndex_SubtrSDP(idx_sig)) ;
                plot(out.FR(idx_sig), y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
                text(8,35, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', condition_colors{i_tsk})
                if i_tsk == 1
                    % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
                else
                    % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
                end
            end
        end
    end
    filename= ['ModulationIndex_Npsc_FiringRate'];
    
    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
    
    %% Correlation between SNR - FR
    hf = figure('Name',sprintf('Correlation_SNR_FiringRate'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for i_BrArea = 1: length(fieldnames(Out))
        for i_tsk = 1: numel(TaskTyp)
            
            out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            % from all nicht NAN units - how many were significant?
            Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
            Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
            
            idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
            
            %% Modulation Index
            ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);        hold on;
            % all not significant & NaN units
            scatter( out.quantSNR(idx_sig),out.FR_perECGTriggeredAverage(idx_sig) , 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
            xlabel('Signal-to-Noise ratio','fontsize',14 );
            ylabel('average Firing rate','fontsize',14 );        axis square; box on;
            title(Ana_TargetBrainArea{i_BrArea})
            set(gca,'ylim',[0, 200]);
            set(gca,'xlim',[0, 25]);
            hold on;
            
            if sum(idx_sig)
                xf = [min(out.quantSNR), max(out.quantSNR)];
                [p,S] = polyfit(out.quantSNR(~isnan(out.FR_perECGTriggeredAverage(:, 1))),out.FR_perECGTriggeredAverage(~isnan(out.FR_perECGTriggeredAverage(:, 1))),1); %
                [y_fit,delta] = polyval(p,out.FR_perECGTriggeredAverage(~isnan(out.FR_perECGTriggeredAverage)),S);
                [coef, pval] = corr(out.FR_perECGTriggeredAverage(~isnan(out.FR_perECGTriggeredAverage)) ,out.quantSNR(~isnan(out.FR_perECGTriggeredAverage)), 'rows','complete') ;
                scatter( out.quantSNR(~idx_sig),out.FR_perECGTriggeredAverage(~idx_sig) , 'filled', 'MarkerFaceColor',condition_colors{i_tsk}/2)
                plot(out.FR_perECGTriggeredAverage(~isnan(out.FR_perECGTriggeredAverage)), y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
                text(8,15, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', condition_colors{i_tsk})
                if i_tsk == 1
                    % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
                else
                    % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
                end
            end
            ha1 = subplot(2,length(fieldnames(Out)),i_BrArea +length(fieldnames(Out)));        hold on;
            
            scatter( out.quantSNR(idx_sig),out.FR_perECGTriggeredAverage(idx_sig) , 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
            xlabel('Signal-to-Noise ratio','fontsize',14 );
            ylabel('average Firing rate','fontsize',14 );        axis square;box on;
            title(Ana_TargetBrainArea{i_BrArea});
            set(gca,'ylim',[0, 200]);
            set(gca,'xlim',[0, 25]);
            hold on;
            
            if sum(idx_sig)
                
                xf = [min(out.quantSNR(idx_sig)), max(out.quantSNR(idx_sig))];
                [p,S] = polyfit(out.quantSNR(idx_sig),out.FR_perECGTriggeredAverage(idx_sig),1); %
                [y_fit,delta] = polyval(p,out.FR_perECGTriggeredAverage(idx_sig),S);
                [coef, pval] = corr(out.quantSNR(idx_sig),out.FR_perECGTriggeredAverage(idx_sig)) ;
                plot(out.FR_perECGTriggeredAverage(idx_sig), y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
                text(8,15, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', condition_colors{i_tsk})
                if i_tsk == 1
                    % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
                else
                    % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
                end
            end
        end
    end
    filename= ['SNR_FiringRate'];
    
    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
    
    %% MODULATION STRENGTH - substractiveNormalization - FR
    % Shows that higher average firing rate is related to a higher
    % modulation index
    hf = figure('Name',sprintf('ModulationIndex_FiringRate'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for i_BrArea = 1: length(fieldnames(Out))
        for i_tsk = 1: numel(TaskTyp)
            
            out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            % from all nicht NAN units - how many were significant?
            Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
            Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
            
            idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
            
            %% Modulation Index
            ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);        hold on;
            % all not significant & NaN units
            if ~isempty(out.quantSNR)
            scatter(out.FR(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
            ylabel('Modulation index (%)','fontsize',14 );
            xlabel('average Firing rate','fontsize',14 );        axis square; box on;
            title(Ana_TargetBrainArea{i_BrArea})
            set(gca,'xlim',[0, 200]);
            set(gca,'ylim',[0, 15]);
            
            hold on;
            
            xf = [min(out.FR), max(out.FR)];
            [p,S] = polyfit(out.FR(~isnan(out.FR_ModIndex_SubtrSDP)),out.FR_ModIndex_SubtrSDP(~isnan(out.FR_ModIndex_SubtrSDP)),1); %
            [y_fit,delta] = polyval(p,out.FR(~isnan(out.FR_ModIndex_SubtrSDP)),S);
            [coef, pval] = corr(out.FR,out.FR_ModIndex_SubtrSDP, 'rows','complete') ;
            scatter(out.FR(~idx_sig) , out.FR_ModIndex_SubtrSDP(~idx_sig), 'filled', 'MarkerFaceColor', condition_colors{i_tsk}/2)
            plot(out.FR(~isnan(out.FR_ModIndex_SubtrSDP)), y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
            text(8,10, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', condition_colors{i_tsk})
            if i_tsk == 1
                % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
            else
                % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
            end
            end
            
            ha1 = subplot(2,length(fieldnames(Out)),i_BrArea +length(fieldnames(Out)));        hold on;
            
            scatter(out.FR(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
            ylabel('Modulation index(%)','fontsize',14 );
            xlabel('average Firing rate','fontsize',14 );        axis square;box on;
            title(Ana_TargetBrainArea{i_BrArea});
            set(gca,'xlim',[0, 200]);        set(gca,'ylim',[0, 15]);
            
            hold on;
            
            if sum(idx_sig)
                
                xf = [min(out.FR(idx_sig)), max(out.FR(idx_sig))];
                [p,S] = polyfit(out.FR(idx_sig),out.FR_ModIndex_SubtrSDP(idx_sig),1); %
                [y_fit,delta] = polyval(p,out.FR(idx_sig),S);
                [coef, pval] = corr(out.FR(idx_sig),out.FR_ModIndex_SubtrSDP(idx_sig)) ;
                plot(out.FR(idx_sig), y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
                text(8,10, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', condition_colors{i_tsk})
                if i_tsk == 1
                    % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
                else
                    % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
                end
            end
        end
    end
    filename= ['ModulationIndex_NSubt_FR'];
    
    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
    
    
    %% MODULATION STRENGTH -  Bin size
    hf = figure('Name',sprintf('ModulationIndex_BinSize'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for i_BrArea = 1: length(fieldnames(Out))
        for i_tsk = 1: numel(TaskTyp)
            
            out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            % from all nicht NAN units - how many were significant?
            Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
            Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
            
            idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
            
            %% Modulation Index
            ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);        hold on;
            % all not significant & NaN units
            if ~isempty(out.quantSNR)
            scatter(out.sig_n_bins(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
            ylabel('Modulation index (%)','fontsize',14 );
            xlabel('Bin size of sig. Interval','fontsize',14 );        axis square; box on;
            title(Ana_TargetBrainArea{i_BrArea})
            set(gca,'xlim',[0, 30]);
            set(gca,'ylim',[0, 60]);
            
            hold on;
            
            xf = [min(out.sig_n_bins), max(out.quantSNR)];
            [p,S] = polyfit(out.sig_n_bins(~isnan(out.FR_ModIndex_PcS)),out.FR_ModIndex_PcS(~isnan(out.FR_ModIndex_PcS)),1); %
            [y_fit,delta] = polyval(p,out.sig_n_bins(~isnan(out.FR_ModIndex_PcS)),S);
            [coef, pval] = corr(out.sig_n_bins,out.FR_ModIndex_PcS, 'rows','complete') ;
            scatter(out.sig_n_bins(~idx_sig) , out.FR_ModIndex_PcS(~idx_sig), 'filled', 'MarkerFaceColor', condition_colors{i_tsk}/2)
            plot(out.sig_n_bins(~isnan(out.FR_ModIndex_PcS)), y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
            text(8,35, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', condition_colors{i_tsk})
            if i_tsk == 1
                % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
            else
                % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
            end
            end
            ha1 = subplot(2,length(fieldnames(Out)),i_BrArea +length(fieldnames(Out)));        hold on;
            
            scatter(out.sig_n_bins(idx_sig) , out.FR_ModIndex_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',condition_colors{i_tsk})
            ylabel('Modulation index(%)','fontsize',14 );
            xlabel('Bin size of sig. Interval','fontsize',14 );        axis square;box on;
            title(Ana_TargetBrainArea{i_BrArea});
            set(gca,'xlim',[0, 30]);
            set(gca,'ylim',[0, 80]);
            
            hold on;
            
            if sum(idx_sig)
                
                xf = [min(out.sig_n_bins(idx_sig)), max(out.quantSNR(idx_sig))];
                [p,S] = polyfit(out.sig_n_bins(idx_sig),out.FR_ModIndex_SubtrSDP(idx_sig),1); %
                [y_fit,delta] = polyval(p,out.sig_n_bins(idx_sig),S);
                [coef, pval] = corr(out.sig_n_bins(idx_sig),out.FR_ModIndex_SubtrSDP(idx_sig)) ;
                plot(out.sig_n_bins(idx_sig), y_fit,'LineWidth', 2, 'Color', condition_colors{i_tsk});
                text(8,35, ['coef, p ', num2str([roundto(coef,2), roundto(pval,3)])], 'Color', condition_colors{i_tsk})
                if i_tsk == 1
                    % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
                else
                    % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
                end
            end
        end
    end
    filename= ['ModulationIndex_Nspc_Binsize'];
    
    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
    
    
    %% overview about all units - change to include the
    %             ThreeTiming = {'T<-50', '-50>T<50', 'T>50'} ;
    %
    %         hf = figure('Name',sprintf('BarPlot'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    %         for i_BrArea = 1: length(fieldnames(Out))
    %                 for i_tsk = 1: numel(TaskTyp)
    %                     for i_Time = 1: length(ThreeTiming)
    %
    %
    %                         out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
    %                         % from all nicht NAN units - how many were significant?
    %                         Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
    %                         Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
    %
    %                         idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
    %                 idx_SigDec      = (out.sig_sign == -1);
    %                 idx_SigInc      = (out.sig_sign == 1);
    %                 idx_SigTime_BeforeMinus50 = (out.sig_time < -50 );
    %                 idx_SigTime_Around0       = (out.sig_time > -50 ) & (out.sig_time < 50) ;
    %                 idx_SigTime_After50       = (out.sig_time > 50);
    %                 idx_Time = [];
    %                 switch i_Time
    %                     case 1 %{ 'BeforeMinus50'}
    %                         idx_Time_Before = (idx_SigTime_BeforeMinus50 & idx_SigInc)+  (idx_SigTime_After50 & idx_SigDec) ;
    %
    %                     case 2 %'Around0'
    %                         idx_Time_During = idx_SigTime_Around0 ;
    %
    %                     case 3 %'After50'
    %                         idx_Time_End = (idx_SigTime_BeforeMinus50 & idx_SigDec)+  (idx_SigTime_After50 & idx_SigInc) ;
    %
    %                 end
    %
    %                         if i_tsk == 1
    %                             Pc_SignFR_rest(i_BrArea,:) = ([sum(out.sig_sign(idx_sig & idx_Time) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] / sum(Idx_Units_NonNaN)) *100;
    %                             Nb_SignFR_rest(i_BrArea,:) = ([sum(out.sig_sign(idx_sig) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] ) ;
    %                             Pc_SignFR_rest2(i_BrArea,:) = round(([sum(idx_sig), (sum(~idx_sig) -Idx_Units_NaN)] / sum(Idx_Units_NonNaN)) *100);
    %                             % sum(Nb_SignFR_rest, 2)
    %                         else
    %                             Pc_SignFR_task(i_BrArea,:) = ([sum(out.sig_sign(idx_sig) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] / sum(Idx_Units_NonNaN)) *100;
    %                             Nb_SignFR_task(i_BrArea,:) = ([sum(out.sig_sign(idx_sig) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] ) ;
    %                             Pc_SignFR_task2(i_BrArea,:) = round(([sum(idx_sig),(sum(~idx_sig) -Idx_Units_NaN),] / sum(Idx_Units_NonNaN)) *100);
    %                             % sum(Nb_SignFR_task, 2)
    %                         end
    %
    %                     end
    %                 end
    %             end
    %         ha1 = subplot(1,2,1);% pie plot how many
    %         barpairs =  [Pc_SignFR_rest];
    %         b = bar(barpairs,'stacked', 'Facecolor','flat' );
    %         title('Rest: non-sig.yellow,iFR-blue,dFR-green','interpreter','none');
    %         set(gca,'XTickLabel',fieldnames(Out),'fontsize',10);
    %
    %         ha1 = subplot(1,2,2);% pie plot how many
    %         barpairs =  [Pc_SignFR_task];
    %         b = bar(barpairs,'stacked', 'Facecolor','flat' );
    %         title('Task: non-sig.yellow,iFR-blue,dFR-green','interpreter','none');
    %         set(gca,'XTickLabel',fieldnames(Out),'fontsize',10);
    %
    %
    %
    %
    %     filename= ['Pc_CardiacRelatedUnits'];
    %
    %     if savePlot;
    %         export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
    %         close all;
    %     end
    %% Decrease and Increase grouped for brain region
    hf = figure('Name',sprintf('CardiacRelated_Change_FR'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    
    % Color_BrainArea = [[0 0 0];  colors(7,:);     colors(13,:); colors(21,:)  ];  %[0 0.9 0.4] %[0 0.6 0] [0.8 0.3 0.1]
    Color_BrainArea = [colors(7,:);     colors(13,:); colors(21,:)  ];
    Color_BrainArea = distinguishable_colors(6);
    Color_BrainArea = [1.0000         0         0;
    1.0000    0.5300         0;
    0.5500    0.9000    0.1000;
    0.0500    0.6500    0.3000;
    0.0500    0.6500    0.7000;
    0.5500    0.2000    0.7500];
    
    for i_tsk = 1: numel(TaskTyp)
        for i_BrArea = 1: length(fieldnames(Out))
            
            out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            % from all nicht NAN units - how many were significant?
            Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
            Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
            
            idx_sig         =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
            idx_SigDec      = (out.sig_sign == -1);
            idx_SigInc      = (out.sig_sign == 1);
            idx_SigTime_BeforeMinus50 = (out.sig_time < -50 );
            idx_SigTime_Around0       = (out.sig_time > -50 ) & (out.sig_time < 50) ;
            idx_SigTime_After50       = (out.sig_time > 50);
            
            
            
            hold on
            if i_tsk == 1
                lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',4};
                subId = 0;
                
            else
                lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',4}; %{'r-o','markerfacecolor','r'}
                subId = 1;
            end
            
            ha1 = subplot(2,4,1 +subId); %
            lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
            text(-400,-1* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), ' ' TaskTyp{i_tsk} ': units = ' ,num2str(sum(idx_SigInc & idx_sig)) ],'Color',Color_BrainArea(i_BrArea,:))
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig,:),1))) ;
            shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig,:),1), SDmean_SEM ,lineProps,1);
            set(gca,'ylim',[-10, 10]);
            title(['units showing a sig. INCREASE in FR'],'interpreter','none');
            ylabel('normalized Firing rate (%)','fontsize',14 );
            xlabel('Time relative to R-peak (ms)','fontsize',14 );axis square; box on;
            set(gca, 'XTick', (ecg_bna_cfg.analyse_states{1,3}*1000):100:(ecg_bna_cfg.analyse_states{1,4}*1000))
            
            
            % display only significant units showing a decrease in FR
            
            ha1 = subplot(2,4,5 +subId); %
            if i_tsk == 1
                lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                text(-400,1* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), 'Rest: units = ' ,num2str(sum(idx_SigDec & idx_sig)) ],'Color',Color_BrainArea(i_BrArea,:))
            else
                lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                text(-400,1* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), 'Task: units = ' ,num2str(sum(idx_SigDec & idx_sig)) ],'Color',Color_BrainArea(i_BrArea,:))
            end
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig,:),1))) ;
            shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig,:),1), SDmean_SEM ,lineProps,1);
            set(gca,'ylim',[-10, 10]);
            title([(TaskTyp{i_tsk}), 'units showing a sig. DECREASE in FR'],'interpreter','none');
            ylabel('normalized Firing rate (spike/s)','fontsize',14 );
            xlabel('Time relative to R-peak (ms)','fontsize',14 );axis square; box on;
            set(gca, 'XTick', (ecg_bna_cfg.analyse_states{1,3}*1000):100:(ecg_bna_cfg.analyse_states{1,4}*1000))
            
            
            
            if i_tsk == 1
                ha4 = subplot(2,4,7); %
                scatter(out.sig_time(idx_SigDec & idx_sig) , out.FR_ModIndex_SubtrSDP(idx_SigDec & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
                ylabel('Modulation strength (% pSc)','fontsize',10);
                xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
                axis square; box on;
                title([(TaskTyp{i_tsk}),'Decrease' ])
                
                ha4 = subplot(2,4,3); %
                scatter(out.sig_time(idx_SigInc & idx_sig) , out.FR_ModIndex_SubtrSDP(idx_SigInc & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
                ylabel('Modulation strength (% pSc)','fontsize',10);
                xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
                axis square; box on;
                title([(TaskTyp{i_tsk}),'Increase' ])
            else
                ha4 = subplot(2,4,8); %
                scatter(out.sig_time(idx_SigDec & idx_sig) , out.FR_ModIndex_SubtrSDP(idx_SigDec & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
                ylabel('Modulation strength (% pSc)','fontsize',10);
                xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
                axis square; box on;
                title([(TaskTyp{i_tsk}),'Decrease' ])
                
                ha4 = subplot(2,4,4); %
                scatter(out.sig_time(idx_SigInc & idx_sig) , out.FR_ModIndex_SubtrSDP(idx_SigInc & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
                ylabel('Modulation strength (% pSc)','fontsize',10);
                xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
                axis square; box on;
                title([(TaskTyp{i_tsk}),'Increase' ])
            end
            
        end
    end
    filename= ['Suppression_Enhancement'];
    
    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
    
    
    
    %% Decrease and Increase grouped for brain region
    hf = figure('Name',sprintf('CardiacRelated_ChangeFR_Time'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    Table_NrUnits = [];
    %     Color_BrainArea = [[0 0 0];  colors(7,:);     colors(13,:); colors(21,:)  ];  %[0 0.9 0.4] %[0 0.6 0] [0.8 0.3 0.1]
    %     Color_BrainArea = [colors(7,:);     colors(13,:); colors(21,:)  ];  %[0 0.9 0.4] %[0 0.6 0] [0.8 0.3 0.1]
    ThreeTiming = {'T<0',  'T>0'} ; %{'T<-50', '-50>T<50', 'T>50'} ;
    c_UPpanel = 1;
    c_Lowpanel = 7;
    for i_tsk = 1: numel(TaskTyp)
        for i_Time = 1: length(ThreeTiming)
            
            for i_BrArea = 1: length(fieldnames(Out))
                
                out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
                % from all nicht NAN units - how many were significant?
                Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP_normalized(:,end));
                Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP_normalized(:,end)));
                
                idx_sig         =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
                idx_SigDec      = (out.sig_sign == -1);
                idx_SigInc      = (out.sig_sign == 1);
                idx_SigTime_BeforeMinus50 = (out.sig_time < 0 );
                idx_SigTime_After50       = (out.sig_time > 0);
                
                %                 idx_SigTime_BeforeMinus50 = (out.sig_time < -50 );
                %                 idx_SigTime_Around0       = (out.sig_time > -50 ) & (out.sig_time < 50) ;
                %                 idx_SigTime_After50       = (out.sig_time > 50);
                A = (ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4});
                
                switch i_Time
                    case 1 %{ 'BeforeMinus50'}
                        idx_Time = idx_SigTime_BeforeMinus50;
                        WindowIdx = find(A < 0);
                        
                        
                    case 2 %'After50'
                        idx_Time = idx_SigTime_After50 ;
                        WindowIdx = find(A > 0);
                    case 1 %'Around0'
                        idx_Time = idx_SigTime_Around0 ;
                        WindowIdx = find(A == 0);
                        
                end
                if i_tsk == 1; subId = 0; else  subId = 6;  end
                ha1 = subplot(2,6,i_Time +subId); %
                
                if i_tsk == 1
                    lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                    %   text(-400,-4* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), ' R: n = ' ,num2str(sum(idx_SigInc & idx_sig & idx_Time)) ],'Color',Color_BrainArea(i_BrArea,:))
                    subId = 0;
                else
                    lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                    %   text(-400,-4* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), 'T: n = ' ,num2str(sum(idx_SigInc & idx_sig & idx_Time)) ],'Color',Color_BrainArea(i_BrArea,:))
                    subId = 6;
                    
                end
                
                if size(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,:),1) < 2 &&  ~size(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,:),1) == 0
                    SDmean_SEM = nan(1, length(nanmean(out.SDsubstractedSDP_normalized)));
                    shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,:), SDmean_SEM ,lineProps,1);
                    MI_groups =  max(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,WindowIdx)) - min(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,:)) ;
                    
                else
                    SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,:)))) ;
                    shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,:)), SDmean_SEM ,lineProps,1);
                    MI_groups =  max(nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,WindowIdx))) - min(nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,:))) ;
                    
                    
                end
                set(gca,'ylim',[-15, 15]);
                title([(TaskTyp{i_tsk}), ' sig.INCREASE ', ThreeTiming(i_Time)],'interpreter','none');
                ylabel('normalized Firing rate (% pSc)','fontsize',14 );
                xlabel('Time relative to R-peak (ms)','fontsize',14 );axis square; box on;
                % vline(50); vline(-50); vline(250); vline(-250);
                % set(gca, 'XTick', (ecg_bna_cfg.analyse_states{1,3}*1000):50:(ecg_bna_cfg.analyse_states{1,4}*1000))
                
                
                %% Calculate the Modulatioin index
                
                %                 Table_Units = table(Ana_TargetBrainArea(i_BrArea),TaskTyp(i_tsk),{'sigINCREASE'}, ThreeTiming(i_Time),  sum(idx_SigInc & idx_sig & idx_Time), roundto(MI_groups,2) );
                %                 Table_NrUnits = [Table_NrUnits; Table_Units];
                
                
                %%
                
                c_UPpanel = c_UPpanel +1;
                % DECREASE -
                
                
                ha1 = subplot(2,6,i_Time +3 +subId); %
                if i_tsk == 1
                    lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                    %   text(-400,4* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), 'R: n = ' ,num2str(sum(idx_SigDec & idx_sig & idx_Time)) ],'Color',Color_BrainArea(i_BrArea,:))
                else
                    lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                    %  text(-400,4* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), 'T: n = ' ,num2str(sum(idx_SigDec & idx_sig & idx_Time)) ],'Color',Color_BrainArea(i_BrArea,:))
                end
                if size(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,:),1) < 2 &&  ~size(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,:),1) == 0
                    SDmean_SEM = nan(1, length(nanmean(out.SDsubstractedSDP_normalized)));
                    shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,:), SDmean_SEM ,lineProps,1);
                    MI_groups =   max(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,:)) - min(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,WindowIdx));
                    
                else
                    
                    SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,:)))) ;
                    shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,:)), SDmean_SEM ,lineProps,1);
                    MI_groups =   max(nanmean(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,:))) - min(nanmean(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,WindowIdx)));
                    
                end
                
                
                
                %                 Table_Units = table(Ana_TargetBrainArea(i_BrArea),TaskTyp(i_tsk),{'sigDECREASE'}, ThreeTiming(i_Time),  sum(idx_SigDec & idx_sig & idx_Time), roundto(MI_groups,2)  );
                %                 Table_NrUnits = [Table_NrUnits; Table_Units];
                
                set(gca,'ylim',[-15, 15]);
                title([(TaskTyp{i_tsk}), ' sig.DECREASE ', ThreeTiming(i_Time)],'interpreter','none');
                ylabel('normalized Firing rate (% pSc)','fontsize',14 );
                xlabel('Time relative to R-peak (ms)','fontsize',14 );axis square; box on;
                %vline(50); vline(-50); vline(250); vline(-250);
                
                %set(gca, 'XTick', (ecg_bna_cfg.analyse_states{1,3}*1000):100:(ecg_bna_cfg.analyse_states{1,4}*1000))
                
                c_Lowpanel = c_Lowpanel +1;
                
                %
                %         if i_tsk == 1
                %             ha4 = subplot(2,4,7); %
                %             scatter(out.sig_time(idx_SigDec & idx_sig) , out.FR_ModIndex_SubtrSDP(idx_SigDec & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
                %             ylabel('Modulation strength (% pSc)','fontsize',10);
                %             xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
                %             axis square; box on;
                %             title([(TaskTyp{i_tsk}),'Decrease' ])
                %
                %             lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                %             SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_BeforeMinus50 & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_BeforeMinus50 & idx_sig,:)))) ;
                %             shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_BeforeMinus50 & idx_sig,:)), SDmean_SEM ,lineProps,1);
                %             set(gca,'ylim',[-10, 10]);
                %
                %             lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                %             SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_Around0 & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_Around0 & idx_sig,:)))) ;
                %             shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_Around0 & idx_sig,:)), SDmean_SEM ,lineProps,1);
                %             set(gca,'ylim',[-10, 10]);
                %
                %             lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                %             SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_After50 & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_After50 & idx_sig,:)))) ;
                %             shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_After50 & idx_sig,:)), SDmean_SEM ,lineProps,1);
                %             set(gca,'ylim',[-10, 10]);
                %
                %
                %
                %
                %
                %
                %
                %
                %             ha4 = subplot(2,4,3); %
                %             scatter(out.sig_time(idx_SigInc & idx_sig) , out.FR_ModIndex_SubtrSDP(idx_SigInc & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
                %             ylabel('Modulation strength (% pSc)','fontsize',10);
                %             xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
                %             axis square; box on;
                %             title([(TaskTyp{i_tsk}),'Increase' ])
                %         else
                %             ha4 = subplot(2,4,4); %
                %             scatter(out.sig_time(idx_SigDec & idx_sig) , out.FR_ModIndex_SubtrSDP(idx_SigDec & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
                %             ylabel('Modulation strength (% pSc)','fontsize',10);
                %             xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
                %             axis square; box on;
                %             title([(TaskTyp{i_tsk}),'Decrease' ])
                %
                %             ha4 = subplot(2,4,8); %
                %             scatter(out.sig_time(idx_SigInc & idx_sig) , out.FR_ModIndex_SubtrSDP(idx_SigInc & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
                %             ylabel('Modulation strength (% pSc)','fontsize',10);
                %             xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
                %             axis square; box on;
                %             title([(TaskTyp{i_tsk}),'Increase' ])
                %         end
                
            end
        end
    end
    
    filename= ['Suppression_Enhancement_SeparatedForTime_NoLINES'];
    
    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
    if saveTable;
        filename = [basepath_to_save , 'Table_NbUnits_IncreaseDecreaseAtSpecificTime.xlsx' ] ;
        writetable(Table_NrUnits,filename,'Sheet',1,  'Range' ,'A1' )
        
    end
    
    %% GENERATE A TABLE FOR ALL THE NUMBER OF UNITS PER CATEGORY
    
    %%
    hf = figure('Name',sprintf('CardiacRelated_ChangeFR_Time'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    
    %     Color_BrainArea = [[0 0 0];  colors(7,:);     colors(13,:); colors(21,:)  ];  %[0 0.9 0.4] %[0 0.6 0] [0.8 0.3 0.1]
    ThreeTiming = {'T<-50', '-50>T<50', 'T>50'} ;
    c_UPpanel = 1;
    c_Lowpanel = 4;
    for i_tsk = 1: numel(TaskTyp)
        for i_Time = 1: length(ThreeTiming)
            
            for i_BrArea = 1: length(fieldnames(Out))
                
                out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
                % from all nicht NAN units - how many were significant?
                Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
                Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
                
                idx_sig         =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
                idx_SigDec      = (out.sig_sign == -1);
                idx_SigInc      = (out.sig_sign == 1);
                idx_SigTime_BeforeMinus50 = (out.sig_time < -50 );
                idx_SigTime_Around0       = (out.sig_time > -50 ) & (out.sig_time < 50) ;
                idx_SigTime_After50       = (out.sig_time > 50);
                idx_Time = [];
                switch i_Time
                    case 1 %{ 'BeforeMinus50'}
                        idx_Time = (idx_SigTime_BeforeMinus50 & idx_SigInc)+  (idx_SigTime_After50 & idx_SigDec) ;
                        
                    case 2 %'Around0'
                        idx_Time = idx_SigTime_Around0 ;
                        
                    case 3 %'After50'
                        idx_Time = (idx_SigTime_BeforeMinus50 & idx_SigDec)+  (idx_SigTime_After50 & idx_SigInc) ;
                        
                end
                if i_tsk == 1; subId = 0; else  subId = 3;  end
                ha1 = subplot(2,3,(i_Time +subId)); %
                
                if i_tsk == 1
                    lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                    text(-400,-4* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), ' R: n = ' ,num2str(sum( idx_sig & idx_Time)) ],'Color',Color_BrainArea(i_BrArea,:))
                    subId = 0;
                else
                    lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                    text(-400,-4* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), 'T: n = ' ,num2str(sum( idx_sig & idx_Time)) ],'Color',Color_BrainArea(i_BrArea,:))
                    subId = 3;
                    
                end
                
                if size(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:),1) < 2 &&  ~size(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:),1) == 0
                    SDmean_SEM = nan(1, length(nanmean(out.SDsubstractedSDP_normalized)));
                    shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:), SDmean_SEM ,lineProps,1);
                    
                else
                    SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:)))) ;
                    shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:)), SDmean_SEM ,lineProps,1);
                end
                set(gca,'ylim',[-15, 15]);
                title([(TaskTyp{i_tsk}), ThreeTiming(i_Time)],'interpreter','none');
                ylabel('normalized Firing rate (% pSc)','fontsize',14 );
                xlabel('Time relative to R-peak (ms)','fontsize',14 );axis square; box on;
                vline(50); vline(-50); vline(250); vline(-250);
                % set(gca, 'XTick', (ecg_bna_cfg.analyse_states{1,3}*1000):50:(ecg_bna_cfg.analyse_states{1,4}*1000))
                
                %                 c_UPpanel = c_UPpanel +1;
                %                 % DECREASE -
                %
                %
                %                 ha1 = subplot(2,3,(i_Time +3) ); %
                %                 if i_tsk == 1
                %                     lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                %                     text(-400,4* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), 'R: n = ' ,num2str(sum( idx_sig & idx_Time)) ],'Color',Color_BrainArea(i_BrArea,:))
                %                 else
                %                     lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                %                     text(-400,4* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), 'T: n = ' ,num2str(sum( idx_sig & idx_Time)) ],'Color',Color_BrainArea(i_BrArea,:))
                %                 end
                %                 if size(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:),1) < 2 &&  ~size(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:),1) == 0
                %                     SDmean_SEM = nan(1, length(nanmean(out.SDsubstractedSDP_normalized)));
                %                     shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:), SDmean_SEM ,lineProps,1);
                %
                %                 else
                %
                %                     SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:)))) ;
                %                     shadedErrorBar((ecg_bna_cfg.analyse_states{1,3}:ecg_bna_cfg.PSTH_binwidth:ecg_bna_cfg.analyse_states{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:)), SDmean_SEM ,lineProps,1);
                %                 end
                %                 set(gca,'ylim',[-15, 15]);
                %                 title([(TaskTyp{i_tsk}), ThreeTiming(i_Time)],'interpreter','none');
                %                 ylabel('normalized Firing rate (% pSc)','fontsize',14 );
                %                 xlabel('Time relative to R-peak (ms)','fontsize',14 );axis square; box on;
                %                 vline(50); vline(-50); vline(250); vline(-250);
                %
                %                 %set(gca, 'XTick', (ecg_bna_cfg.analyse_states{1,3}*1000):100:(ecg_bna_cfg.analyse_states{1,4}*1000))
                %
                %                 c_Lowpanel = c_Lowpanel +1;
                
            end
        end
    end
    
    filename= ['Suppression_Enhancement_SeparatedForTime_GroupedAccordingTo'];
    
    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
    %% Graph for units having rest & task
    
end
end

function rounded=roundto(unrounded,n)
factor= 10^n;
rounded=round(unrounded*factor)/factor;
end
