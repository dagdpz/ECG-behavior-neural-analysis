function ecg_bna_avg_spike_histogram(SPK_PSTH,session_info)
% Here comes some sort of across population plot i assume?

savePlot = 1;
OnlyUnits_withRestANDTask = 0;
Graph_SelectionCriterion = 0; 

ECG_event=-1;
keys.PSTH_WINDOWS={'ECG',ECG_event,-0.5,0.5};
keys.PSTH_binwidth=0.01;
keys.kernel_type='gaussian';
keys.gaussian_kernel=0.02;

basepath_to_save=[session_info(1).SPK_fldr filesep 'Population'];

TaskTyp = {'Rest', 'Task'};

TargetBrainArea = SPK_PSTH{1}.Rest.target(1);
for i = 1: length(SPK_PSTH)
    BrainArea = [SPK_PSTH{i}.Rest.target, SPK_PSTH{i}.Task.target];
    for i_brain = 1 : length(BrainArea)
        if ~ismember(TargetBrainArea, BrainArea(i_brain))
            TargetBrainArea = [TargetBrainArea, BrainArea(i_brain)];
        end
    end
end
TargetBrainArea = sort(TargetBrainArea); 

%% Create function to concatenate the variables per TargetBrainArea

% initate empty structure
for i_BrArea = 1: numel(TargetBrainArea)
    for i_tsk = 1: numel(TaskTyp)
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD             = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD_SEM         = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDP            = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_period     = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR             = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).quantSNR       = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_n_bins     = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_time       = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_sign       = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_FR_diff    = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).unit_ID        = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).target         = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).NrEvents       = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).NrTrials       = [];
        
    end
end

% sort the data according to brain areas
for f_brain = 1: length(TargetBrainArea)
    for i = 1: length(SPK_PSTH)
        for i_tsk = 1: numel(TaskTyp)
            O = [SPK_PSTH{1,i}.(TaskTyp{i_tsk})];
            % index to select only specific units from a brain areas
            idx_brain = ismember(O.target, TargetBrainArea{f_brain});
            
            
            Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SD           = [Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SD;                O.SD(idx_brain,:)];
            Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SDP          = [Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SDP;               O.SDP(idx_brain,:)];
            Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SD_SEM       = [Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SD_SEM;            O.SD_SEM(idx_brain,:)];
            Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_period   = [Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_period;        O.sig(idx_brain,:)];
            Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_FR_diff  = [Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_FR_diff;       O.sig_FR_diff(idx_brain)];
            Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_n_bins   = [Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_n_bins;        O.sig_n_bins(idx_brain)];
            Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_sign     = [Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_sign;          O.sig_sign(idx_brain)];
            Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_time     = [Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_time;          O.sig_time(idx_brain)];
            Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).quantSNR     = [Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).quantSNR;          O.quantSNR(idx_brain)];
            Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).FR           = [Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).FR;                O.FR(idx_brain)];
            Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).unit_ID      = [Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).unit_ID;           O.unit_ID(idx_brain)'];
            Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).target       = [Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).target;            O.target(idx_brain)'];
            Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).NrEvents     = [Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).NrEvents;          O.NrEvents(idx_brain)];
            Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).NrTrials     = [Out.(TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).NrTrials;          O.NrTrials(idx_brain)];
            
            
            
        end
    end
end



%% Crieriua for the dataset
if OnlyUnits_withRestANDTask
    for i_BrArea = 1: length(fieldnames(Out))
        Idx_Units_RestTask = []; Idx_NoSpikes_T = [];  Idx_NoSpikes_R = [];
        
        %  Task = NAN and Rest = NAN - should not exist?
        Idx_NoSpikes_R = find(any(isnan([Out.(TargetBrainArea{i_BrArea}).(TaskTyp{1}).SD]),2));
        Idx_NoSpikes_T = find(any(isnan([Out.(TargetBrainArea{i_BrArea}).(TaskTyp{2}).SD]),2));
        Idx_NoSpikes_RestTask = intersect( Idx_NoSpikes_T,Idx_NoSpikes_R);
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{1}).unit_ID(Idx_NoSpikes_RestTask,:)       %What happenend here?
        
        %     'Bac_20210806_01'
        %     'Bac_20210806_02'
        %     'Bac_20210806_03'
        %     'Bac_20210806_04'
        %     'Bac_20210806_05'
        %     'Bac_20210906_04'
        %     'Bac_20210906_05'
        %     'Bac_20210906_06'
        
        Idx_Spikes_R = find(~any(isnan([Out.(TargetBrainArea{i_BrArea}).(TaskTyp{1}).SD]),2));
        Idx_Spikes_T = find(~any(isnan([Out.(TargetBrainArea{i_BrArea}).(TaskTyp{2}).SD]),2));
        Idx_Units_RestTask = intersect(Idx_Spikes_R, Idx_Spikes_T);
        
        disp([(TargetBrainArea{i_BrArea}),': ', num2str(length(Idx_Units_RestTask))])
        for i_tsk = 1: numel(TaskTyp)
            if i_tsk == 1;  Idx = Idx_NoSpikes_R ; else Idx = Idx_NoSpikes_T ; end
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD(Idx, :)           = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD_SEM(Idx, :)       = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDP(Idx, :)          = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_period(Idx, :)   = [] ;
            
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR(Idx)              = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).quantSNR(Idx)        = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_n_bins(Idx)      = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_time(Idx)        = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_sign(Idx)        = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_FR_diff(Idx)     = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).unit_ID(Idx)         = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).target(Idx)          = [] ;
        end
        
        Idx_Spikes_R = find(~any(isnan([Out.(TargetBrainArea{i_BrArea}).(TaskTyp{1}).SD]),2));
        Idx_Spikes_T = find(~any(isnan([Out.(TargetBrainArea{i_BrArea}).(TaskTyp{2}).SD]),2));
        Idx_Units_RestTask = intersect(Idx_Spikes_R, Idx_Spikes_T);
        disp([(TargetBrainArea{i_BrArea}),': ', num2str(length(Idx_Units_RestTask))])
        
        
    end
end
%% Selection criteria
% criterium 1 spikes per second
Tab_ExcludedUnits = [];
for i_BrArea = 1: length(fieldnames(Out))
    for i_tsk = 1: numel(TaskTyp)
        Criterium_SpkPerSec = 1;
        InVal_idx1 = find([Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR] <= Criterium_SpkPerSec);
        
        Criterium_NrCardiacCycles = 120*5;
        InVal_idx2 = find([Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).NrEvents] <= Criterium_NrCardiacCycles);
        
        InVal_idx = [InVal_idx1(~ismember(InVal_idx1, InVal_idx2)) ; InVal_idx2] ;
        
        
        if ~isempty(InVal_idx)
            InVal_unit_ID = Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).unit_ID(InVal_idx);
            if i_tsk == 1 ; InVal_unit_ID_Rest = InVal_unit_ID;
                idx_IdentRest_Task = 0;
                InVal_Nr = repmat(length(InVal_unit_ID), length(InVal_idx), 1);
                
            else
                idx_IdentRest_Task =  ismember(InVal_unit_ID, InVal_unit_ID_Rest);
                InVal_Nr = repmat(length(InVal_unit_ID(~idx_IdentRest_Task)), length(InVal_idx), 1);
            end
            
            
            TaskType = repmat((TaskTyp(i_tsk)), length(InVal_idx), 1);
            if strcmp((TargetBrainArea{i_BrArea}), 'VPL_R')||   strcmp((TargetBrainArea{i_BrArea}), 'mdT_L') ||   strcmp((TargetBrainArea{i_BrArea}), 'mdT_R')
                BrainArea = repmat({[ '_', (TargetBrainArea{i_BrArea})]}, length(InVal_idx), 1);
            else
                BrainArea = repmat((TargetBrainArea(i_BrArea)), length(InVal_idx), 1);
            end
            Criterium_SpkPerSec = repmat(Criterium_SpkPerSec, length(InVal_idx), 1);
            Criterium_NrCardiacCycles = repmat(Criterium_NrCardiacCycles, length(InVal_idx), 1);
            Nr_InVal_idx2 = repmat(length(InVal_idx2), length(InVal_idx), 1);
            Nr_InVal_idx1 = repmat(length(InVal_idx1), length(InVal_idx), 1);
            
            ExcludedUnits = table(Criterium_SpkPerSec,Nr_InVal_idx1, Criterium_NrCardiacCycles, Nr_InVal_idx2, TaskType,BrainArea, InVal_unit_ID,InVal_idx,InVal_Nr  );
            Tab_ExcludedUnits = [Tab_ExcludedUnits;ExcludedUnits ];
            
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD(InVal_idx, :)           = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD_SEM(InVal_idx, :)       = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDP(InVal_idx, :)          = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_period(InVal_idx, :)   = [] ;
            
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR(InVal_idx)              = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).quantSNR(InVal_idx)        = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_n_bins(InVal_idx)      = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_time(InVal_idx)        = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_sign(InVal_idx)        = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_FR_diff(InVal_idx)     = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).unit_ID(InVal_idx)         = [] ;
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).target(InVal_idx)          = [] ;
            
        end
    end
end

% sum the neurons task + rest for each brain area

%%  Calculations
for i_BrArea = 1: length(fieldnames(Out))
    for i_tsk = 1: numel(TaskTyp)
        out = [Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
        
        idx_sig =  ~isnan(out.sig_FR_diff);
        
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP                     =   out.SD - out.SDP;
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized          =   ((out.SD - out.SDP) ./ out.SDP *100);
        
        % mean(SD - SDP)
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean          =   nanmean(out.SD - out.SDP);
        % standard error of the SDmean
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SEM      =  nanstd(out.SD - out.SDP)/ sqrt(length(Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean )) ;
        
        % Calculate the amount of modulation in FR from the significant interval
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_AllUnits   = zeros(size(Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized,1),1);
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_AllUnits2  = NaN(size(Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized,1),1);
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_Modulation          = NaN(size(Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized,1),1);
        
        % the window of analysis is restricted
        A = (keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4}); 
        A1 = find(A == -0.25)+1; 
        A2 = find(A == 0.25)-1;
        WindowIdx = A1:A2; 
        for i = 1: size(Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized,1)
            Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_AllUnits2(i) = max(Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized(i,WindowIdx)) -  min(Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized(i,WindowIdx));
            
            if any(logical(out.sig_period(i,:)))
%                 switch out.sig_sign(i)
%                     case 1
%                         Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_Modulation(i) = ...
%                          max(Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized(i,logical(out.sig_period(i,:)))) -  ...
%                          min(Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized(i,WindowIdx));
%                     case -1
%                         Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_Modulation(i) = ...
%                             max(Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized(i,WindowIdx)) -  ...
%                             min(Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized(i,logical(out.sig_period(i,:))));
%                 end
%                 
%                 a=Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_Modulation(i) == Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_AllUnits2(i);
%                 
                Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_Modulation(i) = Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_AllUnits2(i);
                Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_AllUnits(i) = Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_Modulation(i);
            end
        end
        
        %% Which units were recorded only during rest?
        % Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).nonSig =  isnan(Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP(:,end))
        
        %   [ Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).nonSig  , Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR]
    end
end

%% Example for 
% i_BrArea = 1; 
%  hf = figure('Name',sprintf(TargetBrainArea{i_BrArea}),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
%     for i_tsk = 1: numel(TaskTyp)
%         
%         out = [Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
%         idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
%         ha1 = subplot(2,4,[1:2]); %
%         if i_tsk == 1
%             lineProps={'color','b','linewidth',4};
%         else
%             lineProps={'color','r','linewidth',4};
%         end
%         
%         
%         % bar plot how many are significant & positive and negative FR?
%         if i_tsk == 1
%             text(-400,-10, ['Rest: units = ' ,num2str(sum(idx_sig)), ' of ' ,num2str(sum(Idx_Units_NonNaN)) ],'Color','b');
%             for i = find(idx_sig)
%                 line((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000 , out.SDsubstractedSDP_normalized(i,:), 'color',[0 0 1 0.2],'LineWidth', 1); hold on;
%             end
%         else
%             text(-400,-12, ['Task: units = ' ,num2str(sum(idx_sig)), ' of ' ,num2str(sum(Idx_Units_NonNaN))  ],'Color','r')
%             
%             for i = find(idx_sig)
%                 line((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000 , out.SDsubstractedSDP_normalized(i,:), 'color',[1 0 0 0.2] ,'LineWidth', 1);hold on;
%             end
%         end
%     end
% 


if Graph_SelectionCriterion
    %%  normalizing the Firing rate calcuations
    for i_BrArea = 1: length(fieldnames(Out))
        hf = figure('Name',sprintf(TargetBrainArea{i_BrArea}),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
        for i_tsk = 1: numel(TaskTyp)
            
            out = [Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            
            idx_sig =  ~isnan(out.sig_FR_diff);
            hold on
            ha1 = subplot(2,2,[1:2]); %
            if i_tsk == 1
                lineProps={'color','b','linewidth',4};
            else
                lineProps={'color','r','linewidth',4};
            end
            
            % bar plot how many are significant & positive and negative FR?
            SDmean_SEM = nanstd(out.SDsubstractedSDP(idx_sig, :))/ sqrt(length(nanmean(out.SDsubstractedSDP(idx_sig, :)))) ;
            shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.SDsubstractedSDP(idx_sig, :)) ,SDmean_SEM ,lineProps,1);
            title(['Pop:  (all significant) Cal Per Unit:SD-SDP' (TargetBrainArea{i_BrArea}) ' units'],'interpreter','none');
            
            ha1 = subplot(2,2,[3:4]); %
            
            if i_tsk == 1
                lineProps={'color','b','linewidth',4};
            else
                lineProps={'color','r','linewidth',4};
            end
            
            
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_sig, :))/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_sig, :)))) ;
            shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_sig, :)) ,SDmean_SEM ,lineProps,1);
            title(['Population:  (all significant)' (TargetBrainArea{i_BrArea}) ' units'],'interpreter','none');
            ylabel('normalized Firing rate (%)','fontsize',14 );
            xlabel('Time relative to R-peak (ms)','fontsize',14 );
            title(['Pop:  (all significant) Cal Per Unit:(SD-SDP)/SDP' (TargetBrainArea{i_BrArea}) ' units'],'interpreter','none');
            
            hold off;
        end
        
        if savePlot
            filename= ['Normalization_FR_ECG_triggered_spike_' (TaskTyp{i_tsk}) , '_' ,(TargetBrainArea{i_BrArea})];
            export_fig([basepath_to_save filesep filename ], '-pdf'); %
            close all;
        end
    end
    
    
    %% Plot the dataset - check up
    for i_BrArea = 1: length(fieldnames(Out))
        for i_tsk = 1: numel(TaskTyp)
            out = [Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            
            hf = figure('Name',sprintf([(TargetBrainArea{i_BrArea}) '  ' TaskTyp{i_tsk}]),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
            rows_plot = 2;
            colums_plot = 7;
            
            for I = 1: 5
                % How many different groups
                if I == 1
                    ha1 = subplot(rows_plot,colums_plot,[4:5]); %
                    
                    idx1 = ([out.sig_n_bins] > 0 );
                    idx2 = ([out.FR_Modulation] > 5 );
                    
                    Y1 = out.SD( idx1 & idx2,:)  ;
                    Y2 = out.SDP( idx1 & idx2,:)  ;
                    
                    if ~length(idx2) == size(Y1,1)
                        disp('incorrect Selection')
                    end
                    title(['FR_Modulation > 5 '],'interpreter','none');
                    
                elseif I == 2
                    ha1 = subplot(rows_plot,colums_plot,[6:7]); %
                    idx1 = ([out.sig_n_bins] > 0 );
                    idx_ex = ([out.FR_Modulation] < 5 );  %Smaller
                    Y1 = out.SD( idx1 & idx_ex,:)  ;
                    Y2 = out.SDP( idx1 & idx_ex,:)  ;
                    title(['FR_Modulation < 5 '],'interpreter','none');      %& Nr. of bin: 0 - 10 bins for sign. Intervals
                    
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
                
                colormap jet;
                cmap=colormap;
                for i = 1: size(Y1(:, 1:end),1)
                    if size(Y1(:, 1:end),1) > 63
                        cmap = [colormap; colormap; colormap; colormap ; colormap ;];
                        t = 1;
                    elseif size(Y1(:, 1:end),1) > 20
                        t = 1;
                        
                    elseif    size(Y1(:, 1:end),1) > 8
                        t = 3;
                    else
                        t = 7;
                    end
                    line((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000 , ( Y1(i,:) - Y2(i,:)), 'color',cmap(i*t,:),'linewidth',3);
                    hold on;
                end
                set(gca,'ylim',[-9, 9]);
                ylabel('Firing rate (%)','fontsize',12 );
                xlabel('Time relative to ECg peak (ms)','fontsize',12);
                
                hold off;
            end
            
            if i_tsk == 1
                Color = [1 0 0];
                
            else
                Color = [0 0 1];
                
            end
            
            idx_ex = ([out.sig_n_bins] <= 4 );
            
            ha4 = subplot(rows_plot,colums_plot,[3]); %
            scatter(out.sig_n_bins , out.FR_Modulation, 'filled', 'MarkerFaceColor',Color); hold on;
            scatter(out.sig_n_bins(idx_ex) , out.FR_Modulation(idx_ex), 'filled', 'MarkerFaceColor',[0 0 0])
            ylabel('FR Modulation (spike/s)','fontsize',14 );
            xlabel('Nr. of bins of sig. Interval','fontsize',14 );
            hold off;  axis square;
            
            ha4 = subplot(rows_plot,colums_plot,[1:2]); %
            scatter(out.sig_n_bins , out.quantSNR, 'filled', 'MarkerFaceColor',Color); hold on;
            scatter(out.sig_n_bins(idx_ex) , out.quantSNR(idx_ex), 'filled', 'MarkerFaceColor',[0 0 0])
            ylabel('Signal-to-Noise Ratio','fontsize',14 );
            xlabel('Nr. of bins of sig. Interval','fontsize',14 );
            hold off;  axis square;
            
            ha4 = subplot(rows_plot,colums_plot,[14]); %
            scatter(out.sig_time , out.FR_Modulation, 'filled', 'MarkerFaceColor',Color); hold on;
            scatter(out.sig_time(idx_ex) , out.FR_Modulation(idx_ex), 'filled', 'MarkerFaceColor',[0 0 0])
            ylabel('FR Modulation (spike/s)','fontsize',10);
            xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
            hold off; axis square;
            
            
            if savePlot
                filename= ['ModulationInRelationToBins_ECG_triggered_spike_' (TaskTyp{i_tsk}) , '_' ,(TargetBrainArea{i_BrArea})];
                export_fig([basepath_to_save filesep filename ], '-pdf'); %
                close all;
            end
        end
    end
end

%% Plot the averages
for i_BrArea = 1: length(fieldnames(Out))
    hf = figure('Name',sprintf(TargetBrainArea{i_BrArea}),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    for i_tsk = 1: numel(TaskTyp)
        
        out = [Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
        % from all nicht NAN units - how many were significant?
        Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
        Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
        
        idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
        %idx_sig = double(idx_sig);
        
        
        hold on
        ha1 = subplot(2,4,[1:2]); %
        if i_tsk == 1
            lineProps={'color','b','linewidth',4};
        else
            lineProps={'color','r','linewidth',4};
        end
        
        
        % bar plot how many are significant & positive and negative FR?
        if i_tsk == 1
            text(-400,-10, ['Rest: units = ' ,num2str(sum(idx_sig)), ' of ' ,num2str(sum(Idx_Units_NonNaN)) ],'Color','b');
            for i = find(idx_sig)
                line((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000 , out.SDsubstractedSDP_normalized(i,:), 'color',[0 0 1 0.2],'LineWidth', 1); hold on;
            end
        else
            text(-400,-12, ['Task: units = ' ,num2str(sum(idx_sig)), ' of ' ,num2str(sum(Idx_Units_NonNaN))  ],'Color','r')
            
            for i = find(idx_sig)
                line((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000 , out.SDsubstractedSDP_normalized(i,:), 'color',[1 0 0 0.2] ,'LineWidth', 1);hold on;
            end
        end
        SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_sig, :))/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_sig, :)))) ;
        shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_sig, :)) ,SDmean_SEM ,lineProps,1);
        title(['Population:  (all significant)' (TargetBrainArea{i_BrArea}) ' units'],'interpreter','none');
        ylabel('normalized Firing rate (%)','fontsize',14 );
        xlabel('Time relative to R-peak (ms)','fontsize',14 );
        hold off;
        
        %Bar graph
        ha1 = subplot(2,4,3);% pie plot how many
        
        if i_tsk == 1
            Pc_SignFR_rest = ([(sum(~idx_sig) -Idx_Units_NaN), sum(out.sig_sign(idx_sig) == -1), sum(out.sig_sign(idx_sig) == 1) ] / sum(Idx_Units_NonNaN)) *100;
            Pc_SignFR_task = ([(sum(~idx_sig)-Idx_Units_NaN), sum(out.sig_sign(idx_sig) == -1), sum(out.sig_sign(idx_sig) == 1) ] / sum(Idx_Units_NonNaN)) *100;
            
        else
            Pc_SignFR_task = ([sum(~idx_sig)-Idx_Units_NaN, sum(out.sig_sign(idx_sig) == -1), sum(out.sig_sign(idx_sig) == 1) ] / sum(Idx_Units_NonNaN)) *100;
        end
        barpairs =  [Pc_SignFR_rest; Pc_SignFR_task];
        b = bar(barpairs,'stacked', 'Facecolor','flat' );
        title('non-sig.blue,iFR-yellow,dFR-green','interpreter','none');
        set(gca,'XTickLabel',{'Rest' , 'Task'},'fontsize',10);
        
        
        % display only significant units showing a increase in FR
        ha1 = subplot(2,4,[5]); %
        idx_SigDec = (out.sig_sign == -1);
        idx_SigInc = (out.sig_sign == 1);
        if i_tsk == 1
            lineProps={'color','b','linewidth',3};
            text(-400,1, ['Rest: units = ' ,num2str(sum(idx_SigInc & idx_sig)) ],'Color','b')
            
        else
            lineProps={'color','r','linewidth',3};
            text(-400,1.5, ['Task: units = ' ,num2str(sum(idx_SigInc & idx_sig)) ],'Color','r')
            
        end
        
        SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig,:)))) ;
        shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig,:)), SDmean_SEM ,lineProps,1);
        set(gca,'ylim',[-5, 10]);
        title(['units showing a sig. INCREASE in FR'],'interpreter','none');
        ylabel('normalized Firing rate (%)','fontsize',14 );
        xlabel('Time relative to R-peak (ms)','fontsize',14 );
        % display only significant units showing a decrease in FR
        ha1 = subplot(2,4,[6]); %
        if i_tsk == 1
            lineProps={'color','b','linewidth',3};
            text(-400,1, ['Rest: units = ' ,num2str(sum(idx_SigDec & idx_sig)) ],'Color','b')
            Color = [0 0 1];
        else
            lineProps={'color','r','linewidth',3};
            text(-400,1.5, ['Task: units = ' ,num2str(sum(idx_SigDec & idx_sig)) ],'Color','r')
            Color = [1 0 0];
        end
        if size(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig,:),1) < 2
            
        else
            
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig,:)))) ;
            
            shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig,:)), SDmean_SEM ,lineProps,1);
        end
        set(gca,'ylim',[-5, 10]);
        title(['units showing a sig. DECREASE in FR'],'interpreter','none');
        ylabel('normalized Firing rate (spike/s)','fontsize',14 );
        xlabel('Time relative to R-peak (ms)','fontsize',14 );
        
        
        %% MODULATION INDEX
        
        ha1 = subplot(2,4,[4]);        hold on;
        % all not significant & NaN units
        % scatter(out.quantSNR(~idx_sig) , out.FR_Modulation(~idx_sig) == 0, 'filled', 'MarkerFaceColor','k')
        scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_AllUnits2(~idx_sig), 'filled', 'MarkerFaceColor','k')
        %  scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_AllUnits(~idx_sig) == 0, 'filled', 'MarkerFaceColor','g')
        
        
        
        scatter(out.quantSNR(idx_sig) , out.FR_Modulation(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('Signal-to-Noise','fontsize',14 );        axis square;
        hold on;
        
        xf = [min(out.quantSNR), max(out.quantSNR)];
        [p,S] = polyfit(out.quantSNR,out.FR_ModIndex_AllUnits2,1); %
        [y_fit,delta] = polyval(p,out.quantSNR,S);
        [coef, pval] = corr(out.quantSNR,out.FR_ModIndex_AllUnits2, 'rows','complete') ;
        if i_tsk == 1
            
            plot(out.quantSNR, y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,35, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            plot(out.quantSNR, y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,30, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r')
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
        end
        
        
        %%
        Dat = [];
        ModIndex = out.FR_Modulation(idx_sig);
        SNR =  out.quantSNR(idx_sig);
        Dat = table(ModIndex, SNR);
        
        ha1 = subplot(2,4,[7]);        hold on;
        
        scatter(out.quantSNR(idx_sig) , out.FR_Modulation(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index(%)','fontsize',14 );
        xlabel('Signal-to-Noise','fontsize',14 );        axis square;
        hold on;
        
        xf = [min(out.quantSNR(idx_sig)), max(out.quantSNR(idx_sig))];
        [p,S] = polyfit(SNR,ModIndex,1); %
        [y_fit,delta] = polyval(p,SNR,S);
        [coef, pval] = corr(SNR,ModIndex) ;
        if i_tsk == 1
            
            plot(SNR, y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,35, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            plot(SNR, y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,30, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r')
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
        end
        
        
        %         mdl = fitlm(Dat);
        %         coefs = mdl.Coefficients.Estimate; % 2x1 [intercept; slope]
        %         pValue = mdl.Coefficients.pValue(2);
        %         Rsquared = mdl.Rsquared;
        %         a =  refline(coefs(2),coefs(1) );
        %         if i_tsk == 1
        %         set(a,'LineWidth', 2, 'Color', 'r');
        %         else
        %         set(a,'LineWidth', 2, 'Color', 'b');
        %         end
        
        
        hold off;
        
        ha1 = subplot(2,4,[8]); %
        scatter(out.FR(idx_sig) , out.FR_Modulation(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index(%)','fontsize',14 );
        xlabel('mean firing rate','fontsize',14);
        title(['all significant units'],'interpreter','none');
        
        axis square;
        
    end
    filename= ['Average_ECG_triggered_spike_' (TargetBrainArea{i_BrArea})];
    
    if savePlot;
        % print(gcf,['Y:\Projects\Pulv_distractor_spatial_choice\ECG_triggered_spikes\ver1\per_unit\Population' filesep, filename,'.pdf'],'-dpdf','-r0');
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
    
end


hf = figure('Name',sprintf('ModulationIndex'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
for i_BrArea = 1: length(fieldnames(Out))
    for i_tsk = 1: numel(TaskTyp)
        
        out = [Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
        % from all nicht NAN units - how many were significant?
        Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
        Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
        
        idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
        if i_tsk == 1

         Color = [0 0 1];
        else
        Color = [1 0 0];

        end
        
        
        %% Modulation Index
        ha1 = subplot(2,5,i_BrArea);        hold on;
        % all not significant & NaN units
        scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_AllUnits2(~idx_sig), 'filled', 'MarkerFaceColor','k')   
        scatter(out.quantSNR(idx_sig) , out.FR_Modulation(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('Signal-to-Noise','fontsize',14 );      %  axis square;
        hold on;
        
        xf = [min(out.quantSNR), max(out.quantSNR)];
        [p,S] = polyfit(out.quantSNR(~isnan(out.FR_ModIndex_AllUnits2)),out.FR_ModIndex_AllUnits2(~isnan(out.FR_ModIndex_AllUnits2)),1); %
        [y_fit,delta] = polyval(p,out.quantSNR(~isnan(out.FR_ModIndex_AllUnits2)),S);
        [coef, pval] = corr(out.quantSNR,out.FR_ModIndex_AllUnits2, 'rows','complete') ;
        if i_tsk == 1
            plot(out.quantSNR(~isnan(out.FR_ModIndex_AllUnits2)), y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,35, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            plot(out.quantSNR(~isnan(out.FR_ModIndex_AllUnits2)), y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,30, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r')
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
        end
        
        Dat = [];
        ModIndex = out.FR_Modulation(idx_sig);
        SNR =  out.quantSNR(idx_sig);
        Dat = table(ModIndex, SNR);
        
        ha1 = subplot(2,5,i_BrArea +length(fieldnames(Out)));        hold on;
        
        scatter(out.quantSNR(idx_sig) , out.FR_Modulation(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index(%)','fontsize',14 );
        xlabel('Signal-to-Noise','fontsize',14 );       % axis square;
        title(TargetBrainArea{i_BrArea})
        hold on;
        
        xf = [min(out.quantSNR(idx_sig)), max(out.quantSNR(idx_sig))];
        [p,S] = polyfit(SNR,ModIndex,1); %
        [y_fit,delta] = polyval(p,SNR,S);
        [coef, pval] = corr(SNR,ModIndex) ;
        if i_tsk == 1
            
            plot(SNR, y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,35, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            plot(SNR, y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,30, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r')
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
        end
    end
end
    filename= ['ModulationIndex'];

    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end

%% overview about all units
hf = figure('Name',sprintf('BarPlot'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
for i_BrArea = 1: length(fieldnames(Out))
    for i_tsk = 1: numel(TaskTyp)
        
        out = [Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
        % from all nicht NAN units - how many were significant?
        Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
        Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
        
        idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
        %idx_sig = double(idx_sig);
        
        
        if i_tsk == 1
            Pc_SignFR_rest(i_BrArea,:) = ([sum(out.sig_sign(idx_sig) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] / sum(Idx_Units_NonNaN)) *100;            
        else
            Pc_SignFR_task(i_BrArea,:) = ([sum(~idx_sig)-Idx_Units_NaN, sum(out.sig_sign(idx_sig) == -1), sum(out.sig_sign(idx_sig) == 1) ] / sum(Idx_Units_NonNaN)) *100;
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
        
    filename= ['Pc_CardiacRelatedUnits'];

    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end 

%% Graph for units having rest & task


