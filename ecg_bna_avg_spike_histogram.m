function ecg_bna_avg_spike_histogram(SPK_PSTH,session_info)
% Here comes some sort of across population plot i assume?

savePlot = 1;
OnlyUnits_withRestANDTask = 0;
Graph_SelectionCriterion = 0; 
    colors = distinguishable_colors(25); 

    
ECG_event=-1;
keys.PSTH_WINDOWS={'ECG',ECG_event,-0.5,0.5};
keys.PSTH_binwidth=0.01;
keys.kernel_type='gaussian';
keys.gaussian_kernel=0.02;

basepath_to_save=[session_info(1).SPK_fldr filesep 'Population'];
if ~exist(basepath_to_save,'dir')
    mkdir(basepath_to_save);
end

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
Ana_TargetBrainArea = TargetBrainArea; 

%% Create function to concatenate the variables per TargetBrainArea
Out = []; 
% initate empty structure
for i_BrArea = 1: numel(TargetBrainArea)
    if strcmp((TargetBrainArea{i_BrArea}), 'mdT_R') || strcmp((TargetBrainArea{i_BrArea}), 'mdT_L')
        Ana_TargetBrainArea{i_BrArea} = 'MD'; 
    end
    for i_tsk = 1: numel(TaskTyp)
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD             = [];
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD_SEM         = [];
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDP            = [];
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_period     = [];
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR             = [];
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).quantSNR       = [];
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_n_bins     = [];
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_time       = [];
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_sign       = [];
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_FR_diff    = [];
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).unit_ID        = [];
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).target         = [];
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).NrEvents       = [];
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).NrTrials       = [];
        
    end
end

% sort the data according to brain areas
for f_brain = 1: length(TargetBrainArea)

    for i = 1: length(SPK_PSTH)
        for i_tsk = 1: numel(TaskTyp)
            O = [SPK_PSTH{1,i}.(TaskTyp{i_tsk})];
            % index to select only specific units from a brain areas
            idx_brain = ismember(O.target, TargetBrainArea{f_brain})';

            
            
            Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SD           = [Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SD;                O.SD(idx_brain,:)];
            Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SDP          = [Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SDP;               O.SDP(idx_brain,:)];
            Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SD_SEM       = [Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).SD_SEM;            O.SD_SEM(idx_brain,:)];
            Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_period   = [Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_period;        O.sig(idx_brain,:)];
            Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_FR_diff  = [Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_FR_diff;       O.sig_FR_diff(idx_brain)];
            Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_n_bins   = [Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_n_bins;        O.sig_n_bins(idx_brain)];
            Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_sign     = [Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_sign;          O.sig_sign(idx_brain)];
            Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_time     = [Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).sig_time;          O.sig_time(idx_brain)];
            Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).quantSNR     = [Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).quantSNR;          O.quantSNR(idx_brain)];
            Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).FR           = [Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).FR;                O.FR(idx_brain)];
            Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).unit_ID      = [Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).unit_ID;           O.unit_ID(idx_brain)'];
            Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).target       = [Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).target;            O.target(idx_brain)'];
            Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).NrEvents     = [Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).NrEvents;          O.NrEvents(idx_brain)];
            Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).NrTrials     = [Out.(Ana_TargetBrainArea{f_brain}).(TaskTyp{i_tsk}).NrTrials;          O.NrTrials(idx_brain)];
            
            
            
        end
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
        
        %     'Bac_20210806_01'
        %     'Bac_20210806_02'
        %     'Bac_20210806_03'
        %     'Bac_20210806_04'
        %     'Bac_20210806_05'
        %     'Bac_20210906_04'
        %     'Bac_20210906_05'
        %     'Bac_20210906_06'
        
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
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_n_bins(Idx)      = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_time(Idx)        = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_sign(Idx)        = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_FR_diff(Idx)     = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).unit_ID(Idx)         = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).target(Idx)          = [] ;
        end
        
        Idx_Spikes_R = find(~any(isnan([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{1}).SD]),2));
        Idx_Spikes_T = find(~any(isnan([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{2}).SD]),2));
        Idx_Units_RestTask = intersect(Idx_Spikes_R, Idx_Spikes_T);
        disp([(Ana_TargetBrainArea{i_BrArea}),': ', num2str(length(Idx_Units_RestTask))])
        
        
    end
end
%% Selection criteria
% criterium 1 spikes per second
Tab_ExcludedUnits = [];
for i_BrArea = 1: length(fieldnames(Out))
    for i_tsk = 1: numel(TaskTyp)
        Criterium_SpkPerSec = 2;
        InVal_idx1 = find([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR] <= Criterium_SpkPerSec);
        
        Criterium_NrCardiacCycles = 120*5;
        InVal_idx2 = find([Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).NrEvents] <= Criterium_NrCardiacCycles);
        
        InVal_idx = [InVal_idx1(~ismember(InVal_idx1, InVal_idx2)) ; InVal_idx2] ;
        
        
        if ~isempty(InVal_idx)
            InVal_unit_ID = Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).unit_ID(InVal_idx);
            if i_tsk == 1 ; InVal_unit_ID_Rest = InVal_unit_ID;
                idx_IdentRest_Task = 0;
                InVal_Nr = repmat(length(InVal_unit_ID), length(InVal_idx), 1);
                
            else
                idx_IdentRest_Task =  ismember(InVal_unit_ID, InVal_unit_ID_Rest);
                InVal_Nr = repmat(length(InVal_unit_ID(~idx_IdentRest_Task)), length(InVal_idx), 1);
            end
            
            
            TaskType = repmat((TaskTyp(i_tsk)), length(InVal_idx), 1);
            if strcmp((Ana_TargetBrainArea{i_BrArea}), 'VPL_R')||   strcmp((Ana_TargetBrainArea{i_BrArea}), 'mdT_L') ||   strcmp((Ana_TargetBrainArea{i_BrArea}), 'mdT_R')
                BrainArea = repmat({[ '_', (Ana_TargetBrainArea{i_BrArea})]}, length(InVal_idx), 1);
            else
                BrainArea = repmat((Ana_TargetBrainArea(i_BrArea)), length(InVal_idx), 1);
            end
            Criterium_SpkPerSec = repmat(Criterium_SpkPerSec, length(InVal_idx), 1);
            Criterium_NrCardiacCycles = repmat(Criterium_NrCardiacCycles, length(InVal_idx), 1);
            Nr_InVal_idx2 = repmat(length(InVal_idx2), length(InVal_idx), 1);
            Nr_InVal_idx1 = repmat(length(InVal_idx1), length(InVal_idx), 1);
            
            ExcludedUnits = table(Criterium_SpkPerSec,Nr_InVal_idx1, Criterium_NrCardiacCycles, Nr_InVal_idx2, TaskType,BrainArea, InVal_unit_ID,InVal_idx,InVal_Nr  );
            Tab_ExcludedUnits = [Tab_ExcludedUnits;ExcludedUnits ];
            
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD(InVal_idx, :)           = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD_SEM(InVal_idx, :)       = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDP(InVal_idx, :)          = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_period(InVal_idx, :)   = [] ;
            
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR(InVal_idx)              = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).quantSNR(InVal_idx)        = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_n_bins(InVal_idx)      = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_time(InVal_idx)        = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_sign(InVal_idx)        = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_FR_diff(InVal_idx)     = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).unit_ID(InVal_idx)         = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).target(InVal_idx)          = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).NrEvents(InVal_idx)        = [] ;
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).NrTrials(InVal_idx)        = [] ;
            
        end
    end
end

save([basepath_to_save,filesep , 'Table_excludedUnits' ],'Tab_ExcludedUnits');
filename = [basepath_to_save,filesep , 'Table_excludedUnits.xlsx' ];  
writetable(Tab_ExcludedUnits,filename,'Sheet',1,  'Range' ,'A1' )
disp(['SAVED   ', basepath_to_save,filesep , 'Table_excludedUnits' ])


%%  Calculations
for i_BrArea = 1: length(fieldnames(Out))
    for i_tsk = 1: numel(TaskTyp)
        out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
        
        idx_sig =  ~isnan(out.sig_FR_diff);
        
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP                     =   out.SD - out.SDP;
        
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized          =   ((out.SD - out.SDP) ./ out.SDP *100);
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_perECGTriggeredAverage          = nanmean(out.SD,2); 
        % mean(SD - SDP)
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean          =   nanmean(out.SD - out.SDP);
        % standard error of the SDmean
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SEM      =  nanstd(out.SD - out.SDP)/ sqrt(length(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean )) ;
        
        % Calculate the amount of modulation in FR from the significant interval
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_AllUnits   = zeros(size(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized,1),1);
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_AllUnits_PcS  = NaN(size(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized,1),1);
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_Modulation          = NaN(size(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized,1),1);
       
        Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_AllUnits_SubtrSDP  = NaN(size(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP,1),1);
        
        % the window of analysis is restricted
        A = (keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4}); 
        A1 = find(A == -0.25)+1; 
        A2 = find(A == 0.25)-1;
        WindowIdx = A1:A2; 
        for i = 1: size(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized,1)
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_AllUnits_PcS(i) = max(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized(i,WindowIdx)) -  min(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized(i,WindowIdx));
            Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_AllUnits_SubtrSDP(i) = max(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP(i,WindowIdx)) -  min(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP(i,WindowIdx));
            
            if any(logical(out.sig_period(i,:)))
%                 switch out.sig_sign(i)
%                     case 1
%                         Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_Modulation(i) = ...
%                          max(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized(i,logical(out.sig_period(i,:)))) -  ...
%                          min(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized(i,WindowIdx));
%                     case -1
%                         Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_Modulation(i) = ...
%                             max(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized(i,WindowIdx)) -  ...
%                             min(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP_normalized(i,logical(out.sig_period(i,:))));
%                 end
%                 
%                 a=Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_Modulation(i) == Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_AllUnits2(i);
%                 
                Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_Modulation(i) = Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_AllUnits_PcS(i);
                Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_ModIndex_AllUnits(i) = Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR_Modulation(i);
            end
        end
        
        %% Which units were recorded only during rest?
        % Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).nonSig =  isnan(Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDsubstractedSDP(:,end))
        
        %   [ Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).nonSig  , Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR]
    end
end




 %% How much does the surrogate and mean firing divergate
 hf = figure('Name',sprintf('Surrogate_MeanFr'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
 for i_BrArea = 1: length(fieldnames(Out))
     for i_tsk = 1: numel(TaskTyp)
         out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
         
        if i_tsk == 1 ; %: numel(TaskTyp)
            Color  = [0 0 1]; 
        else
            Color  = [1 0 0]; 
        end
         ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);        hold on;
         scatter(out.FR , nanmean(out.SDP,2), 'filled', 'MarkerFaceColor',Color)  %,30, out.FR(idx_sig)/max(out.FR) , 'filled')
         set(gca,'ylim',[0,max(out.FR)]);
         set(gca,'xlim',[0,max(out.FR)]);
         x = linspace(0,max(out.FR)); y = linspace(0,max(out.FR)); hold on; 
            plot(x,y, 'k');
         ylabel('mean surrogate','fontsize',14 );
         xlabel('average Firing rate','fontsize',14 );        axis square; box on;
         title([ (TaskTyp{i_tsk}), Ana_TargetBrainArea{i_BrArea}]);
         
         
          ha1 = subplot(2,length(fieldnames(Out)), (i_BrArea + length(fieldnames(Out))));        hold on;
         scatter(nanmean(out.SD,2) , nanmean(out.SDP,2), 'filled', 'MarkerFaceColor',Color)  %,30, out.FR(idx_sig)/max(out.FR) , 'filled')
         set(gca,'ylim',[0,max(nanmean(out.SD,2))]);
         set(gca,'xlim',[0,max(nanmean(out.SD,2))]);
         x = linspace(0,max(nanmean(out.SD,2))); y = linspace(0,max(nanmean(out.SD,2))); hold on; 
            plot(x,y, 'k');
         ylabel('mean surrogate','fontsize',14 );
         xlabel('average Firing rate of ECG-triggeerd Average','fontsize',14 );        axis square; box on;
         title([  Ana_TargetBrainArea{i_BrArea}]);
         
     end
 end

    filename= ['MeanSurrogate_MeanFr'];

    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
 
 
%% Example for 
i_BrArea = 3; 
    for i_tsk = 1: numel(TaskTyp)   
  hf = figure('Name',sprintf(Ana_TargetBrainArea{i_BrArea}),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');

        if i_tsk == 1 ; %: numel(TaskTyp)
            Color  = [0 0 1]; 
        else
            Color  = [1 0 0]; 
        end
        
        out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
        idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
        
        ha1 = subplot(2,4,[1]);        hold on;

        scatter(out.quantSNR(idx_sig) , out.FR_Modulation(idx_sig), 'filled', 'MarkerFaceColor',Color)        
        if i_tsk == 1 ; 
            scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_AllUnits_PcS(~idx_sig), 'filled', 'MarkerFaceColor',[0 0 0])
        else
            scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_AllUnits_PcS(~idx_sig), 'filled', 'MarkerFaceColor',[0.5 0 0.2])

        end

        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('Signal-to-Noise','fontsize',14 );        axis square;
        hold on;
        title([Ana_TargetBrainArea{i_BrArea} , ' ' ,(TaskTyp{i_tsk}) ])
        
        xf = [min(out.quantSNR), max(out.quantSNR)];
        [p,S] = polyfit(out.quantSNR(~isnan(out.FR_ModIndex_AllUnits_PcS)),out.FR_ModIndex_AllUnits_PcS(~isnan(out.FR_ModIndex_AllUnits_PcS)),1); %
        [y_fit,delta] = polyval(p,out.quantSNR(~isnan(out.FR_ModIndex_AllUnits_PcS)),S);
        [coef, pval] = corr(out.quantSNR,out.FR_ModIndex_AllUnits_PcS, 'rows','complete') ;
          if i_tsk == 1 ; 
        plot(out.quantSNR(~isnan(out.FR_ModIndex_AllUnits_PcS)), y_fit,'LineWidth', 2, 'Color', 'b');
        text(8,25, ['coef, p ', num2str([round(coef,2), round(pval,4)])], 'Color', 'b')
          else
        plot(out.quantSNR(~isnan(out.FR_ModIndex_AllUnits_PcS)), y_fit,'LineWidth', 2, 'Color', 'r');
        text(8,30, ['coef, p ', num2str([round(coef,2), round(pval,4)])], 'Color', 'r')
          end

        ha1 = subplot(2,4,[2]);        hold on;
        scatter(out.FR(idx_sig) , out.FR_Modulation(idx_sig), 'filled', 'MarkerFaceColor',Color)
          if i_tsk == 1 ; 
        scatter(out.FR(~idx_sig) , out.FR_ModIndex_AllUnits_PcS(~idx_sig), 'filled', 'MarkerFaceColor',[0 0 0])
          else
        scatter(out.FR(~idx_sig) , out.FR_ModIndex_AllUnits_PcS(~idx_sig), 'filled', 'MarkerFaceColor',[0.5 0 0.2])

          end
        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('mean FR','fontsize',14 );        axis square;
        hold on;
        
        
        ha1 = subplot(2,4,[3]);        hold on;
        scatter(out.sig_n_bins(idx_sig) , out.FR_Modulation(idx_sig), 'filled', 'MarkerFaceColor',Color)
         if i_tsk == 1 ; 
        scatter(out.sig_n_bins(~idx_sig) , out.FR_ModIndex_AllUnits_PcS(~idx_sig), 'filled', 'MarkerFaceColor',[0 0 0])
         else
        scatter(out.sig_n_bins(~idx_sig) , out.FR_ModIndex_AllUnits_PcS(~idx_sig), 'filled', 'MarkerFaceColor',[0.5 0 0.2])

         end
        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('sig. Nr. bins','fontsize',14 );        axis square;
        hold off;
        
        
        ha1 = subplot(2,4,[5:6]); %
        UnitSig_Rest = intersect(find(out.FR_Modulation >30) , find(idx_sig)); 
        UnitNotSign_Rest = intersect(find(out.FR_ModIndex_AllUnits_PcS >30) , find(~idx_sig)); 

        
            text(-400,-15, [out.unit_ID(UnitSig_Rest)],'Color',Color);
            lineProps={'color','b','linewidth',4};
            for i = UnitSig_Rest
                if ~ isempty(i)
                line((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000 , out.SDsubstractedSDP_normalized(i,:), 'color',Color,'LineWidth', 1); hold on;
                end
            end
            hold on; 
            
            colormap bone;
            cmap=colormap;                idx = 1;
            
            lineProps={'color','k','linewidth',4};
            text(300,20, [out.unit_ID(UnitNotSign_Rest)],'Color','k');
            col = cmap([1;30],:);
            for i = 1: length(UnitNotSign_Rest)
                idx_unit = UnitNotSign_Rest(i);
                line((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000 , out.SDsubstractedSDP_normalized(idx_unit,:), 'color',cmap(i,:),'LineWidth', 4); hold on;
            end
            title('Units with Modulation strength > 30%'); %axis square;
            
            ha1 = subplot(2,4,[7:8]); %
            text(300,20, [out.unit_ID(UnitNotSign_Rest)],'Color','k');
            for i = 1: length(UnitNotSign_Rest)
                idx_unit = UnitNotSign_Rest(i);
                line((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000 , out.SDsubstractedSDP_normalized(idx_unit,:), 'color',cmap(i,:),'LineWidth', 4); hold on;
            end
            
            title('Non-Sign units with Modulation strength > 30%'); hold off; 
    if savePlot;
        export_fig([basepath_to_save,filesep ,['Check_ModulationIndex_',Ana_TargetBrainArea{i_BrArea},'_',(TaskTyp{i_tsk}) ]], '-pdf'); %,'-transparent'
        close all;
    end
    end

if Graph_SelectionCriterion
    %%  normalizing the Firing rate calcuations
    for i_BrArea = 1: length(fieldnames(Out))
        hf = figure('Name',sprintf(Ana_TargetBrainArea{i_BrArea}),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
        for i_tsk = 1: numel(TaskTyp)
            
            out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            
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
            title(['Pop:  (all significant) Cal Per Unit:SD-SDP' (Ana_TargetBrainArea{i_BrArea}) ' units'],'interpreter','none');
            
            ha1 = subplot(2,2,[3:4]); %
            
            if i_tsk == 1
                lineProps={'color','b','linewidth',4};
            else
                lineProps={'color','r','linewidth',4};
            end
            
            
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_sig, :))/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_sig, :)))) ;
            shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_sig, :)) ,SDmean_SEM ,lineProps,1);
            title(['Population:  (all significant)' (Ana_TargetBrainArea{i_BrArea}) ' units'],'interpreter','none');
            ylabel('normalized Firing rate (%)','fontsize',14 );
            xlabel('Time relative to R-peak (ms)','fontsize',14 );
            title(['Pop:  (all significant) Cal Per Unit:(SD-SDP)/SDP' (Ana_TargetBrainArea{i_BrArea}) ' units'],'interpreter','none');
            
            hold off;
        end
        
        if savePlot
            filename= ['Normalization_FR_ECG_triggered_spike_' (TaskTyp{i_tsk}) , '_' ,(Ana_TargetBrainArea{i_BrArea})];
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
        title(['Population:  (all significant)' (Ana_TargetBrainArea{i_BrArea}) ' units'],'interpreter','none');
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
        title('iFR-blue,dFR-green, non-sig.yellow,','interpreter','none');
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
        set(gca,'ylim',[-10, 10]);
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
        set(gca,'ylim',[-10, 10]);
        title(['units showing a sig. DECREASE in FR'],'interpreter','none');
        ylabel('normalized Firing rate (spike/s)','fontsize',14 );
        xlabel('Time relative to R-peak (ms)','fontsize',14 );
        
        
        %% MODULATION INDEX
        
        ha1 = subplot(2,4,[4]);        hold on;
        % all not significant & NaN units
        % scatter(out.quantSNR(~idx_sig) , out.FR_Modulation(~idx_sig) == 0, 'filled', 'MarkerFaceColor','k')
        %  scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_AllUnits(~idx_sig) == 0, 'filled', 'MarkerFaceColor','g')
        
        
        
        scatter(out.quantSNR(idx_sig) , out.FR_Modulation(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('Signal-to-Noise','fontsize',14 );        axis square;
        hold on;
        
        xf = [min(out.quantSNR), max(out.quantSNR)];
        [p,S] = polyfit(out.quantSNR,out.FR_ModIndex_AllUnits_PcS,1); %
        [y_fit,delta] = polyval(p,out.quantSNR,S);
        [coef, pval] = corr(out.quantSNR,out.FR_ModIndex_AllUnits_PcS, 'rows','complete') ;
        if i_tsk == 1
           scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_AllUnits_PcS(~idx_sig), 'filled', 'MarkerFaceColor',[0.4 0.3 0.99])

            plot(out.quantSNR, y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,35, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_AllUnits_PcS(~idx_sig), 'filled', 'MarkerFaceColor',[0 0.4 1])

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
        
%         
%         xf = [min(out.quantSNR(idx_sig)), max(out.quantSNR(idx_sig))];
%         [p,S] = polyfit(SNR,ModIndex,1); %
%         [y_fit,delta] = polyval(p,SNR,S);
%         [coef, pval] = corr(SNR,ModIndex) ;
%         if i_tsk == 1
%             
%             plot(SNR, y_fit,'LineWidth', 2, 'Color', 'b');
%             text(8,35, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
%             % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
%         else
%             plot(SNR, y_fit,'LineWidth', 2, 'Color', 'r');
%             text(8,30, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r')
%             % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
%         end
        
        ylabel('Modulation index(%)','fontsize',14 );
        xlabel('mean firing rate','fontsize',14);
        title(['all significant units'],'interpreter','none');
        
        axis square;
        
    end
    filename= ['Average_ECG_triggered_spike_' (Ana_TargetBrainArea{i_BrArea})];
    
    if savePlot;
        % print(gcf,['Y:\Projects\Pulv_distractor_spatial_choice\ECG_triggered_spikes\ver1\per_unit\Population' filesep, filename,'.pdf'],'-dpdf','-r0');
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
        if i_tsk == 1

        Color = [0 0 1];
        else
        Color = [1 0 0];

        end
        
        
        %% Modulation Index
        ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);        hold on;
        % all not significant & NaN units
        scatter(out.quantSNR(idx_sig) , out.FR_Modulation(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('Signal-to-Noise','fontsize',14 );        axis square; box on; 
        title(Ana_TargetBrainArea{i_BrArea})
        set(gca,'xlim',[0, 25]);
        set(gca,'ylim',[0, 80]);
        hold on;
        
        xf = [min(out.quantSNR), max(out.quantSNR)];
        [p,S] = polyfit(out.quantSNR(~isnan(out.FR_ModIndex_AllUnits_PcS)),out.FR_ModIndex_AllUnits_PcS(~isnan(out.FR_ModIndex_AllUnits_PcS)),1); %
        [y_fit,delta] = polyval(p,out.quantSNR(~isnan(out.FR_ModIndex_AllUnits_PcS)),S);
        [coef, pval] = corr(out.quantSNR,out.FR_ModIndex_AllUnits_PcS, 'rows','complete') ;
        if i_tsk == 1
            scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_AllUnits_PcS(~idx_sig), 'filled', 'MarkerFaceColor',colors(11,:))   

            plot(out.quantSNR(~isnan(out.FR_ModIndex_AllUnits_PcS)), y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,35, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_AllUnits_PcS(~idx_sig), 'filled', 'MarkerFaceColor',colors(16,:))   

            plot(out.quantSNR(~isnan(out.FR_ModIndex_AllUnits_PcS)), y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,30, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r')
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
        end
        
        Dat = [];
        ModIndex = out.FR_Modulation(idx_sig);
        SNR =  out.quantSNR(idx_sig);
        Dat = table(ModIndex, SNR);
        
        ha1 = subplot(2,length(fieldnames(Out)),i_BrArea +length(fieldnames(Out)));        hold on;
        
        scatter(out.quantSNR(idx_sig) , out.FR_Modulation(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index(%)','fontsize',14 );
        xlabel('Signal-to-Noise','fontsize',14 );        axis square;box on; 
        title(Ana_TargetBrainArea{i_BrArea}); 
        set(gca,'xlim',[0, 25]);
        set(gca,'ylim',[0, 80]);

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
        if i_tsk == 1

        Color = [0 0 1];
        else
        Color = [1 0 0];

        end
        
        
        %% Modulation Index
        ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);        hold on;
        % all not significant & NaN units
        scatter(out.quantSNR(idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index (spike/s)','fontsize',14 );
        xlabel('Signal-to-Noise ratio (mV)','fontsize',14 );        axis square; box on; 
        title(Ana_TargetBrainArea{i_BrArea})
        set(gca,'xlim',[0, 25]);
        set(gca,'ylim',[0, 15]);


        hold on;
        
        xf = [min(out.quantSNR), max(out.quantSNR)];
        [p,S] = polyfit(out.quantSNR(~isnan(out.FR_ModIndex_AllUnits_SubtrSDP)),out.FR_ModIndex_AllUnits_SubtrSDP(~isnan(out.FR_ModIndex_AllUnits_SubtrSDP)),1); %
        [y_fit,delta] = polyval(p,out.quantSNR(~isnan(out.FR_ModIndex_AllUnits_SubtrSDP)),S);
        [coef, pval] = corr(out.quantSNR,out.FR_ModIndex_AllUnits_SubtrSDP, 'rows','complete') ;
        if i_tsk == 1
            scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(~idx_sig), 'filled', 'MarkerFaceColor',colors(11,:))   

            plot(out.quantSNR(~isnan(out.FR_ModIndex_AllUnits_SubtrSDP)), y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,10, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            scatter(out.quantSNR(~idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(~idx_sig), 'filled', 'MarkerFaceColor',colors(16,:))   

            plot(out.quantSNR(~isnan(out.FR_ModIndex_AllUnits_SubtrSDP)), y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,8, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r')
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
        end

        
        ha1 = subplot(2,length(fieldnames(Out)),i_BrArea +length(fieldnames(Out)));        hold on;
        
        scatter(out.quantSNR(idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index(spike/s)','fontsize',14 );
        xlabel('Signal-to-Noise ratio (mV)','fontsize',14 );        axis square;box on; 
        title(Ana_TargetBrainArea{i_BrArea}); 
        set(gca,'xlim',[0, 25]);

        set(gca,'ylim',[0, 15]);
 
  
        hold on;
        
        xf = [min(out.quantSNR(idx_sig)), max(out.quantSNR(idx_sig))];
        [p,S] = polyfit(out.quantSNR(idx_sig),out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig),1); %
        [y_fit,delta] = polyval(p,out.quantSNR(idx_sig),S);
        
        [coef, pval] = corr(out.quantSNR(idx_sig),out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig)) ;
        if i_tsk == 1
            
            plot(out.quantSNR(idx_sig), y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,10, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b','fontsize',6)
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            plot(out.quantSNR(idx_sig), y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,8, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r' ,'fontsize',6)
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
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
        if i_tsk == 1

        Color = [0 0 1];
        else
        Color = [1 0 0];

        end
        
        
        %% Modulation Index
        ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);        hold on;
        % all not significant & NaN units
        scatter(out.FR_ModIndex_AllUnits_PcS(idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',Color)              
        ylabel('Modulation index (spike/s)','fontsize',14 );
        xlabel('Modulation index(%pSc)','fontsize',14 );        axis square; box on; 
        title(Ana_TargetBrainArea{i_BrArea})
        set(gca,'xlim',[0, 80]);
        set(gca,'ylim',[0, 15]);


        hold on;
        
        xf = [min(out.quantSNR), max(out.quantSNR)];
        [p,S] = polyfit(out.FR_ModIndex_AllUnits_PcS(~isnan(out.FR_ModIndex_AllUnits_SubtrSDP)),out.FR_ModIndex_AllUnits_SubtrSDP(~isnan(out.FR_ModIndex_AllUnits_SubtrSDP)),1); %
        [y_fit,delta] = polyval(p,out.FR_ModIndex_AllUnits_PcS(~isnan(out.FR_ModIndex_AllUnits_SubtrSDP)),S);
        [coef, pval] = corr(out.FR_ModIndex_AllUnits_PcS,out.FR_ModIndex_AllUnits_SubtrSDP, 'rows','complete') ;
        if i_tsk == 1
            scatter(out.FR_ModIndex_AllUnits_PcS(~idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(~idx_sig), 'filled', 'MarkerFaceColor',colors(11,:))   

            plot(out.FR_ModIndex_AllUnits_PcS(~isnan(out.FR_ModIndex_AllUnits_SubtrSDP)), y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,10, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b','fontsize',14)
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            scatter(out.FR_ModIndex_AllUnits_PcS(~idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(~idx_sig), 'filled', 'MarkerFaceColor',colors(16,:))   

            plot(out.FR_ModIndex_AllUnits_PcS(~isnan(out.FR_ModIndex_AllUnits_SubtrSDP)), y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,8, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r','fontsize',14)
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
        end

        
        ha1 = subplot(2,length(fieldnames(Out)),i_BrArea +length(fieldnames(Out)));        hold on;        
        scatter(out.FR_ModIndex_AllUnits_PcS(idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index(spike/s)','fontsize',14 );
        xlabel('Modulation index(%pSc)','fontsize',14 );        axis square;box on; 
        title(Ana_TargetBrainArea{i_BrArea}); 
        set(gca,'xlim',[0, 80]);

        set(gca,'ylim',[0, 15]);
 
  
        hold on;
        
        xf = [min(out.FR_ModIndex_AllUnits_PcS(idx_sig)), max(out.quantSNR(idx_sig))];
        [p,S] = polyfit(out.FR_ModIndex_AllUnits_PcS(idx_sig),out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig),1); %
        [y_fit,delta] = polyval(p,out.FR_ModIndex_AllUnits_PcS(idx_sig),S);
        
        [coef, pval] = corr(out.FR_ModIndex_AllUnits_PcS(idx_sig),out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig)) ;
        if i_tsk == 1
            
            plot(out.FR_ModIndex_AllUnits_PcS(idx_sig), y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,10, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b','fontsize',14)
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            plot(out.FR_ModIndex_AllUnits_PcS(idx_sig), y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,8, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r' ,'fontsize',14)
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
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
       % scatter(out.FR_ModIndex_AllUnits_PcS(idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',Color)                
        scatter(out.FR_ModIndex_AllUnits_PcS(idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig),30, out.FR_perECGTriggeredAverage(idx_sig)/max(out.FR_perECGTriggeredAverage(idx_sig)) , 'filled')
        %scatter(out.FR_ModIndex_AllUnits_PcS(~idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(~idx_sig),30, out.FR(~idx_sig)/max(out.FR) , 'filled')   

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
        if i_tsk == 2
        ha1 = subplot(2,length(fieldnames(Out)),(i_BrArea + 4));        hold on;
        % all not significant & NaN units
       % scatter(out.FR_ModIndex_AllUnits_PcS(idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',Color)                
        scatter(out.FR_ModIndex_AllUnits_PcS(idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig),30, out.FR_perECGTriggeredAverage(idx_sig)/max(out.FR_perECGTriggeredAverage(idx_sig)) , 'filled')
       % scatter(out.FR_ModIndex_AllUnits_PcS(~idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(~idx_sig),30, out.FR(~idx_sig)/max(out.FR) , 'filled')   

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
    filename= ['ModulationIndex_Compare_NpSc_Subtr_WithFR_ONLYSignificantUnits'];

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
        if i_tsk == 1

        Color = [0 0 1];
        else
        Color = [1 0 0];

        end
        
        
        %% Modulation Index
        ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);        hold on;
        % all not significant & NaN units
        scatter(out.FR(idx_sig) , out.FR_Modulation(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('average Firing rate','fontsize',14 );        axis square; box on; 
        title(Ana_TargetBrainArea{i_BrArea})
        set(gca,'xlim',[0, 200]);
        
        set(gca,'ylim',[0, 80]);
       

        hold on;
        
        xf = [min(out.FR), max(out.quantSNR)];
        [p,S] = polyfit(out.FR(~isnan(out.FR_ModIndex_AllUnits_PcS)),out.FR_ModIndex_AllUnits_PcS(~isnan(out.FR_ModIndex_AllUnits_PcS)),1); %
        [y_fit,delta] = polyval(p,out.FR(~isnan(out.FR_ModIndex_AllUnits_PcS)),S);
        [coef, pval] = corr(out.FR,out.FR_ModIndex_AllUnits_PcS, 'rows','complete') ;
        if i_tsk == 1
            scatter(out.FR(~idx_sig) , out.FR_ModIndex_AllUnits_PcS(~idx_sig), 'filled', 'MarkerFaceColor',colors(11,:))   

            plot(out.FR(~isnan(out.FR_ModIndex_AllUnits_PcS)), y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,35, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            scatter(out.FR(~idx_sig) , out.FR_ModIndex_AllUnits_PcS(~idx_sig), 'filled', 'MarkerFaceColor',colors(16,:))   

            plot(out.FR(~isnan(out.FR_ModIndex_AllUnits_PcS)), y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,30, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r')
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
        end
   ha1 = subplot(2,length(fieldnames(Out)),i_BrArea +length(fieldnames(Out)));        hold on;
        
        scatter(out.FR(idx_sig) , out.FR_Modulation(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index(%)','fontsize',14 );
        xlabel('average Firing rate','fontsize',14 );        axis square;box on; 
        title(Ana_TargetBrainArea{i_BrArea}); 
        set(gca,'xlim',[0, 200]);
        set(gca,'ylim',[0, 80]);

        hold on;
        
        xf = [min(out.FR(idx_sig)), max(out.quantSNR(idx_sig))];
        [p,S] = polyfit(out.FR(idx_sig),out.FR_Modulation(idx_sig),1); %
        [y_fit,delta] = polyval(p,out.FR(idx_sig),S);
        [coef, pval] = corr(out.FR(idx_sig),out.FR_Modulation(idx_sig)) ;
        if i_tsk == 1
            
            plot(out.FR(idx_sig), y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,35, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            plot(out.FR(idx_sig), y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,30, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r')
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
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
        if i_tsk == 1

        Color = [0 0 1];
        else
        Color = [1 0 0];

        end
        
        
        %% Modulation Index
        ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);        hold on;
        % all not significant & NaN units
        scatter(out.FR(idx_sig) , out.quantSNR(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Signal-to-Noise ratio','fontsize',14 );
        xlabel('average Firing rate','fontsize',14 );        axis square; box on; 
        title(Ana_TargetBrainArea{i_BrArea})
        set(gca,'xlim',[0, 200]);
        set(gca,'ylim',[0, 25]);
        hold on;
        
        xf = [min(out.FR), max(out.quantSNR)];
        [p,S] = polyfit(out.FR(~isnan(out.quantSNR)),out.quantSNR(~isnan(out.quantSNR)),1); %
        [y_fit,delta] = polyval(p,out.FR(~isnan(out.quantSNR)),S);
        [coef, pval] = corr(out.FR,out.quantSNR, 'rows','complete') ;
        if i_tsk == 1
            scatter(out.FR(~idx_sig) , out.quantSNR(~idx_sig), 'filled', 'MarkerFaceColor',colors(11,:))   

            plot(out.FR(~isnan(out.quantSNR)), y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,15, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            scatter(out.FR(~idx_sig) , out.quantSNR(~idx_sig), 'filled', 'MarkerFaceColor',colors(16,:))   

            plot(out.FR(~isnan(out.quantSNR)), y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,20, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r')
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
        end
   ha1 = subplot(2,length(fieldnames(Out)),i_BrArea +length(fieldnames(Out)));        hold on;
        
        scatter(out.FR(idx_sig) , out.quantSNR(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Signal-to-Noise ratio','fontsize',14 );
        xlabel('average Firing rate','fontsize',14 );        axis square;box on; 
        title(Ana_TargetBrainArea{i_BrArea}); 
        set(gca,'xlim',[0, 200]);

        set(gca,'ylim',[0, 25]);

        hold on;
        
        xf = [min(out.FR(idx_sig)), max(out.quantSNR(idx_sig))];
        [p,S] = polyfit(out.FR(idx_sig),out.quantSNR(idx_sig),1); %
        [y_fit,delta] = polyval(p,out.FR(idx_sig),S);
        [coef, pval] = corr(out.FR(idx_sig),out.quantSNR(idx_sig)) ;
        if i_tsk == 1
            
            plot(out.FR(idx_sig), y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,15, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            plot(out.FR(idx_sig), y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,20, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r')
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
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
        if i_tsk == 1

        Color = [0 0 1];
        else
        Color = [1 0 0];

        end
        
        
        %% Modulation Index
        ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);        hold on;
        % all not significant & NaN units
        scatter(out.FR(idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('average Firing rate','fontsize',14 );        axis square; box on; 
        title(Ana_TargetBrainArea{i_BrArea})
        set(gca,'xlim',[0, 200]);
        set(gca,'ylim',[0, 15]);



        hold on;
        
        xf = [min(out.FR), max(out.FR)];
        [p,S] = polyfit(out.FR(~isnan(out.FR_ModIndex_AllUnits_SubtrSDP)),out.FR_ModIndex_AllUnits_SubtrSDP(~isnan(out.FR_ModIndex_AllUnits_SubtrSDP)),1); %
        [y_fit,delta] = polyval(p,out.FR(~isnan(out.FR_ModIndex_AllUnits_SubtrSDP)),S);
        [coef, pval] = corr(out.FR,out.FR_ModIndex_AllUnits_SubtrSDP, 'rows','complete') ;
        if i_tsk == 1
            scatter(out.FR(~idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(~idx_sig), 'filled', 'MarkerFaceColor',colors(11,:))   

            plot(out.FR(~isnan(out.FR_ModIndex_AllUnits_SubtrSDP)), y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,10, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            scatter(out.FR(~idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(~idx_sig), 'filled', 'MarkerFaceColor',colors(16,:))   

            plot(out.FR(~isnan(out.FR_ModIndex_AllUnits_SubtrSDP)), y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,12, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r')
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
        end
        
        
        ha1 = subplot(2,length(fieldnames(Out)),i_BrArea +length(fieldnames(Out)));        hold on;
        
        scatter(out.FR(idx_sig) , out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index(%)','fontsize',14 );
        xlabel('average Firing rate','fontsize',14 );        axis square;box on; 
        title(Ana_TargetBrainArea{i_BrArea}); 
        set(gca,'xlim',[0, 200]);        set(gca,'ylim',[0, 15]);

        hold on;
        
        xf = [min(out.FR(idx_sig)), max(out.FR(idx_sig))];
        [p,S] = polyfit(out.FR(idx_sig),out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig),1); %
        [y_fit,delta] = polyval(p,out.FR(idx_sig),S);
        [coef, pval] = corr(out.FR(idx_sig),out.FR_ModIndex_AllUnits_SubtrSDP(idx_sig)) ;
        if i_tsk == 1
            
            plot(out.FR(idx_sig), y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,10, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            plot(out.FR(idx_sig), y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,12, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r')
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
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
        if i_tsk == 1

        Color = [0 0 1];
        else
        Color = [1 0 0];

        end
        
        
        %% Modulation Index
        ha1 = subplot(2,length(fieldnames(Out)),i_BrArea);        hold on;
        % all not significant & NaN units
        scatter(out.sig_n_bins(idx_sig) , out.FR_Modulation(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index (%)','fontsize',14 );
        xlabel('Bin size of sig. Interval','fontsize',14 );        axis square; box on; 
        title(Ana_TargetBrainArea{i_BrArea})
        set(gca,'xlim',[0, 30]);
        set(gca,'ylim',[0, 60]);


        hold on;
        
        xf = [min(out.sig_n_bins), max(out.quantSNR)];
        [p,S] = polyfit(out.sig_n_bins(~isnan(out.FR_ModIndex_AllUnits_PcS)),out.FR_ModIndex_AllUnits_PcS(~isnan(out.FR_ModIndex_AllUnits_PcS)),1); %
        [y_fit,delta] = polyval(p,out.sig_n_bins(~isnan(out.FR_ModIndex_AllUnits_PcS)),S);
        [coef, pval] = corr(out.sig_n_bins,out.FR_ModIndex_AllUnits_PcS, 'rows','complete') ;
        if i_tsk == 1
            scatter(out.sig_n_bins(~idx_sig) , out.FR_ModIndex_AllUnits_PcS(~idx_sig), 'filled', 'MarkerFaceColor',colors(11,:))   

            plot(out.sig_n_bins(~isnan(out.FR_ModIndex_AllUnits_PcS)), y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,35, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            scatter(out.sig_n_bins(~idx_sig) , out.FR_ModIndex_AllUnits_PcS(~idx_sig), 'filled', 'MarkerFaceColor',colors(16,:))   

            plot(out.sig_n_bins(~isnan(out.FR_ModIndex_AllUnits_PcS)), y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,30, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r')
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
        end
   ha1 = subplot(2,length(fieldnames(Out)),i_BrArea +length(fieldnames(Out)));        hold on;
        
        scatter(out.sig_n_bins(idx_sig) , out.FR_Modulation(idx_sig), 'filled', 'MarkerFaceColor',Color)
        ylabel('Modulation index(%)','fontsize',14 );
        xlabel('Bin size of sig. Interval','fontsize',14 );        axis square;box on; 
        title(Ana_TargetBrainArea{i_BrArea}); 
        set(gca,'xlim',[0, 30]);
        set(gca,'ylim',[0, 80]);

        hold on;
        
        xf = [min(out.sig_n_bins(idx_sig)), max(out.quantSNR(idx_sig))];
        [p,S] = polyfit(out.sig_n_bins(idx_sig),out.FR_Modulation(idx_sig),1); %
        [y_fit,delta] = polyval(p,out.sig_n_bins(idx_sig),S);
        [coef, pval] = corr(out.sig_n_bins(idx_sig),out.FR_Modulation(idx_sig)) ;
        if i_tsk == 1
            
            plot(out.sig_n_bins(idx_sig), y_fit,'LineWidth', 2, 'Color', 'b');
            text(8,35, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'b')
            % plot(SNR,y_fit+2*delta,'r--',SNR,y_fit-2*delta,'r--')
        else
            plot(out.sig_n_bins(idx_sig), y_fit,'LineWidth', 2, 'Color', 'r');
            text(8,30, ['coef, p ', num2str([round(coef,2), round(pval,3)])], 'Color', 'r')
            % plot(SNR,y_fit+2*delta,'b-',SNR,y_fit-2*delta,'b-')
        end
    end
end
    filename= ['ModulationIndex_Nspc_Binsize'];

    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end
%% overview about all units
hf = figure('Name',sprintf('BarPlot'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
for i_BrArea = 1: length(fieldnames(Out))
    for i_tsk = 1: numel(TaskTyp)
            
            out = [Out.(Ana_TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            Idx_Units_NonNaN = ~isnan(out.SDsubstractedSDP(:,end));
            Idx_Units_NaN =  sum(isnan(out.SDsubstractedSDP(:,end)));
            idx_sig =  ~isnan(out.sig_FR_diff) & (out.sig_n_bins > 4) ;
            
            if i_tsk == 1
                Pc_SignFR_rest(i_BrArea,:) = ([sum(out.sig_sign(idx_sig) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] / sum(Idx_Units_NonNaN)) *100;
                Nb_SignFR_rest(i_BrArea,:) = ([sum(out.sig_sign(idx_sig) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] ) ;
                Pc_SignFR_rest2(i_BrArea,:) = round(([sum(idx_sig), (sum(~idx_sig) -Idx_Units_NaN)] / sum(Idx_Units_NonNaN)) *100);
                % sum(Nb_SignFR_rest, 2)
            else
                Pc_SignFR_task(i_BrArea,:) = ([sum(out.sig_sign(idx_sig) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] / sum(Idx_Units_NonNaN)) *100;
                Nb_SignFR_task(i_BrArea,:) = ([sum(out.sig_sign(idx_sig) == 1) ,sum(out.sig_sign(idx_sig) == -1),(sum(~idx_sig) -Idx_Units_NaN),] ) ;
                Pc_SignFR_task2(i_BrArea,:) = round(([sum(idx_sig),(sum(~idx_sig) -Idx_Units_NaN),] / sum(Idx_Units_NonNaN)) *100);
                % sum(Nb_SignFR_task, 2)
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

Color_BrainArea = [[0 0 0];  colors(7,:);     colors(13,:); colors(21,:)  ];  %[0 0.9 0.4] %[0 0.6 0] [0.8 0.3 0.1]

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

        if i_tsk == 1
            lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
            text(-400,-1* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), ' Rest: units = ' ,num2str(sum(idx_SigInc & idx_sig)) ],'Color',Color_BrainArea(i_BrArea,:))
            
        else
            lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
            text(-400,-1* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), 'Task: units = ' ,num2str(sum(idx_SigInc & idx_sig)) ],'Color',Color_BrainArea(i_BrArea,:))
            
        end
        
        
        
        SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig,:)))) ;
        shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig,:)), SDmean_SEM ,lineProps,1);
        set(gca,'ylim',[-10, 10]);
        title(['units showing a sig. INCREASE in FR'],'interpreter','none');
        ylabel('normalized Firing rate (%)','fontsize',14 );
        xlabel('Time relative to R-peak (ms)','fontsize',14 );axis square; box on;
        set(gca, 'XTick', (keys.PSTH_WINDOWS{1,3}*1000):100:(keys.PSTH_WINDOWS{1,4}*1000))
        
        
        % display only significant units showing a decrease in FR
        
        ha1 = subplot(2,4,5 +subId); %
        if i_tsk == 1
            lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
            text(-400,1* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), 'Rest: units = ' ,num2str(sum(idx_SigDec & idx_sig)) ],'Color',Color_BrainArea(i_BrArea,:))
        else
            lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
            text(-400,1* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), 'Task: units = ' ,num2str(sum(idx_SigDec & idx_sig)) ],'Color',Color_BrainArea(i_BrArea,:))
        end
        if size(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig,:),1) < 2
            
        else
            
            SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig,:)))) ;
            shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig,:)), SDmean_SEM ,lineProps,1);
        end
        set(gca,'ylim',[-10, 10]);
        title([(TaskTyp{i_tsk}), 'units showing a sig. DECREASE in FR'],'interpreter','none');
        ylabel('normalized Firing rate (spike/s)','fontsize',14 );
        xlabel('Time relative to R-peak (ms)','fontsize',14 );axis square; box on;
        set(gca, 'XTick', (keys.PSTH_WINDOWS{1,3}*1000):100:(keys.PSTH_WINDOWS{1,4}*1000))
        
        
        
        if i_tsk == 1
            ha4 = subplot(2,4,7); %
            scatter(out.sig_time(idx_SigDec & idx_sig) , out.FR_Modulation(idx_SigDec & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
            ylabel('Modulation strength (% pSc)','fontsize',10);
            xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
            axis square; box on;
            title([(TaskTyp{i_tsk}),'Decrease' ])
             
            ha4 = subplot(2,4,3); %
            scatter(out.sig_time(idx_SigInc & idx_sig) , out.FR_Modulation(idx_SigInc & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
            ylabel('Modulation strength (% pSc)','fontsize',10);
            xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
            axis square; box on;
            title([(TaskTyp{i_tsk}),'Increase' ])
        else
            ha4 = subplot(2,4,8); %
            scatter(out.sig_time(idx_SigDec & idx_sig) , out.FR_Modulation(idx_SigDec & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
            ylabel('Modulation strength (% pSc)','fontsize',10);
            xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
            axis square; box on;
            title([(TaskTyp{i_tsk}),'Decrease' ])
            
            ha4 = subplot(2,4,4); %
            scatter(out.sig_time(idx_SigInc & idx_sig) , out.FR_Modulation(idx_SigInc & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
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
    
    Color_BrainArea = [[0 0 0];  colors(7,:);     colors(13,:); colors(21,:)  ];  %[0 0.9 0.4] %[0 0.6 0] [0.8 0.3 0.1]
    ThreeTiming = {'T<-50', '-50>T<50', 'T>50'} ;
    c_UPpanel = 1;
    c_Lowpanel = 7;
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
                
                switch i_Time
                    case 1 %{ 'BeforeMinus50'}
                        idx_Time = idx_SigTime_BeforeMinus50;
                    case 2 %'Around0'
                        idx_Time = idx_SigTime_Around0 ;
                        
                    case 3 %'After50'
                        idx_Time = idx_SigTime_After50 ;
                        
                end
                if i_tsk == 1; subId = 0; else  subId = 6;  end
                ha1 = subplot(2,6,i_Time +subId); %

                if i_tsk == 1
                    lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                    text(-400,-4* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), ' R: n = ' ,num2str(sum(idx_SigInc & idx_sig & idx_Time)) ],'Color',Color_BrainArea(i_BrArea,:))
                    subId = 0;
                else
                    lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                    text(-400,-4* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), 'T: n = ' ,num2str(sum(idx_SigInc & idx_sig & idx_Time)) ],'Color',Color_BrainArea(i_BrArea,:))
                    subId = 6;
                    
                end
                
                if size(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,:),1) < 2 &&  ~size(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,:),1) == 0
                    SDmean_SEM = nan(1, length(nanmean(out.SDsubstractedSDP_normalized)));
                    shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,:), SDmean_SEM ,lineProps,1);
                    
                else
                    SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,:)))) ;
                    shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_sig & idx_Time,:)), SDmean_SEM ,lineProps,1);
                end
                set(gca,'ylim',[-15, 15]);
                title([(TaskTyp{i_tsk}), ' sig.INCREASE ', ThreeTiming(i_Time)],'interpreter','none');
                ylabel('normalized Firing rate (% pSc)','fontsize',14 );
                xlabel('Time relative to R-peak (ms)','fontsize',14 );axis square; box on;
                vline(50); vline(-50); vline(250); vline(-250);
               % set(gca, 'XTick', (keys.PSTH_WINDOWS{1,3}*1000):50:(keys.PSTH_WINDOWS{1,4}*1000))
                
                c_UPpanel = c_UPpanel +1;
                % DECREASE -
                
                
                ha1 = subplot(2,6,i_Time +3 +subId); %
                if i_tsk == 1
                    lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                    text(-400,4* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), 'R: n = ' ,num2str(sum(idx_SigDec & idx_sig & idx_Time)) ],'Color',Color_BrainArea(i_BrArea,:))
                else
                    lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                    text(-400,4* i_BrArea, [(Ana_TargetBrainArea{i_BrArea}), 'T: n = ' ,num2str(sum(idx_SigDec & idx_sig & idx_Time)) ],'Color',Color_BrainArea(i_BrArea,:))
                end
                if size(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,:),1) < 2 &&  ~size(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,:),1) == 0
                    SDmean_SEM = nan(1, length(nanmean(out.SDsubstractedSDP_normalized)));
                    shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,:), SDmean_SEM ,lineProps,1);
                    
                else
                    
                    SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,:)))) ;
                    shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigDec & idx_sig & idx_Time,:)), SDmean_SEM ,lineProps,1);
                end
                set(gca,'ylim',[-15, 15]);
                title([(TaskTyp{i_tsk}), ' sig.DECREASE ', ThreeTiming(i_Time)],'interpreter','none');
                ylabel('normalized Firing rate (% pSc)','fontsize',14 );
                xlabel('Time relative to R-peak (ms)','fontsize',14 );axis square; box on;
                vline(50); vline(-50); vline(250); vline(-250);

                %set(gca, 'XTick', (keys.PSTH_WINDOWS{1,3}*1000):100:(keys.PSTH_WINDOWS{1,4}*1000))
                
                c_Lowpanel = c_Lowpanel +1;
                
                %
                %         if i_tsk == 1
                %             ha4 = subplot(2,4,7); %
                %             scatter(out.sig_time(idx_SigDec & idx_sig) , out.FR_Modulation(idx_SigDec & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
                %             ylabel('Modulation strength (% pSc)','fontsize',10);
                %             xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
                %             axis square; box on;
                %             title([(TaskTyp{i_tsk}),'Decrease' ])
                %
                %             lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                %             SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_BeforeMinus50 & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_BeforeMinus50 & idx_sig,:)))) ;
                %             shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_BeforeMinus50 & idx_sig,:)), SDmean_SEM ,lineProps,1);
                %             set(gca,'ylim',[-10, 10]);
                %
                %             lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                %             SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_Around0 & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_Around0 & idx_sig,:)))) ;
                %             shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_Around0 & idx_sig,:)), SDmean_SEM ,lineProps,1);
                %             set(gca,'ylim',[-10, 10]);
                %
                %             lineProps={'color',Color_BrainArea(i_BrArea,:),'linewidth',3};
                %             SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_After50 & idx_sig,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_After50 & idx_sig,:)))) ;
                %             shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized(idx_SigInc & idx_SigTime_After50 & idx_sig,:)), SDmean_SEM ,lineProps,1);
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
                %             scatter(out.sig_time(idx_SigInc & idx_sig) , out.FR_Modulation(idx_SigInc & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
                %             ylabel('Modulation strength (% pSc)','fontsize',10);
                %             xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
                %             axis square; box on;
                %             title([(TaskTyp{i_tsk}),'Increase' ])
                %         else
                %             ha4 = subplot(2,4,4); %
                %             scatter(out.sig_time(idx_SigDec & idx_sig) , out.FR_Modulation(idx_SigDec & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
                %             ylabel('Modulation strength (% pSc)','fontsize',10);
                %             xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
                %             axis square; box on;
                %             title([(TaskTyp{i_tsk}),'Decrease' ])
                %
                %             ha4 = subplot(2,4,8); %
                %             scatter(out.sig_time(idx_SigInc & idx_sig) , out.FR_Modulation(idx_SigInc & idx_sig), 'filled', 'MarkerFaceColor',Color_BrainArea(i_BrArea,:)); hold on;
                %             ylabel('Modulation strength (% pSc)','fontsize',10);
                %             xlabel('TinePoint of sig. highest diff in FR','fontsize',10 );
                %             axis square; box on;
                %             title([(TaskTyp{i_tsk}),'Increase' ])
                %         end
                
            end
        end
      end
    
   filename= ['Suppression_Enhancement_SeparatedForTime'];

    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end 

    %% 
    hf = figure('Name',sprintf('CardiacRelated_ChangeFR_Time'),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    
    Color_BrainArea = [[0 0 0];  colors(7,:);     colors(13,:); colors(21,:)  ];  %[0 0.9 0.4] %[0 0.6 0] [0.8 0.3 0.1]
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
                    shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:), SDmean_SEM ,lineProps,1);
                    
                else
                    SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:)))) ;
                    shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:)), SDmean_SEM ,lineProps,1);
                end
                set(gca,'ylim',[-15, 15]);
                title([(TaskTyp{i_tsk}), ThreeTiming(i_Time)],'interpreter','none');
                ylabel('normalized Firing rate (% pSc)','fontsize',14 );
                xlabel('Time relative to R-peak (ms)','fontsize',14 );axis square; box on;
                vline(50); vline(-50); vline(250); vline(-250);
               % set(gca, 'XTick', (keys.PSTH_WINDOWS{1,3}*1000):50:(keys.PSTH_WINDOWS{1,4}*1000))
                
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
%                     shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:), SDmean_SEM ,lineProps,1);
%                     
%                 else
%                     
%                     SDmean_SEM = nanstd(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:),0,1)/ sqrt(length(nanmean(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:)))) ;
%                     shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.SDsubstractedSDP_normalized( idx_sig & idx_Time,:)), SDmean_SEM ,lineProps,1);
%                 end
%                 set(gca,'ylim',[-15, 15]);
%                 title([(TaskTyp{i_tsk}), ThreeTiming(i_Time)],'interpreter','none');
%                 ylabel('normalized Firing rate (% pSc)','fontsize',14 );
%                 xlabel('Time relative to R-peak (ms)','fontsize',14 );axis square; box on;
%                 vline(50); vline(-50); vline(250); vline(-250);
% 
%                 %set(gca, 'XTick', (keys.PSTH_WINDOWS{1,3}*1000):100:(keys.PSTH_WINDOWS{1,4}*1000))
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
    

