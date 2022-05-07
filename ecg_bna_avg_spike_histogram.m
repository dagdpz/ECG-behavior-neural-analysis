function ecg_bna_avg_spike_histogram(SPK_PSTH,session_info)
% Here comes some sort of across population plot i assume?


ECG_event=-1;
keys.PSTH_WINDOWS={'ECG',ECG_event,-0.5,0.5};
keys.PSTH_binwidth=0.01;
keys.kernel_type='gaussian';
keys.gaussian_kernel=0.02;

basepath_to_save=[session_info(1).SPK_fldr filesep 'Population'];

TaskTyp = {'Rest', 'Task'};
figure;
TargetBrainArea = fieldnames(SPK_PSTH{1}); 

for i = 1: length(SPK_PSTH)
    BrainArea = fieldnames(SPK_PSTH{i});
    if ~ismember(TargetBrainArea, BrainArea)
        TargetBrainArea = [TargetBrainArea, BrainArea];
    end
end

%% Create function to concatenate the variables per TargetBrainArea

% initate empty structure
for i_BrArea = 1: numel(TargetBrainArea)
    for i_tsk = 1: numel(TaskTyp)
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD_SEM  = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDP = [];
        
    end
end
for i = 1: length(SPK_PSTH)
   BrainArea = fieldnames(SPK_PSTH{i});
    for i_BrArea = 1: numel(BrainArea)
        for i_tsk = 1: numel(TaskTyp)
            O = [SPK_PSTH{1,i}.(BrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            
            NrUnit = size(Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD, 1) ; 
            
            SD      = [O.SD,    repmat(NrUnit,  size(O.SD,1),1)+ (1:size(O.SD,1))',    repmat(i_tsk,  size(O.SD,1),1) ,    repmat(i_BrArea,  size(O.SD,1),1) , O.sig_FR_diff, O.sig_sign , O.sig_n_bins];
            SDP     = [O.SDP,   repmat(NrUnit,  size(O.SD,1),1)+ (1:size(O.SD,1))',    repmat(i_tsk,  size(O.SDP,1),1) ,   repmat(i_BrArea,  size(O.SDP,1),1) , O.sig_FR_diff, O.sig_sign , O.sig_n_bins];
            SD_SEM  = [O.SD_SEM,repmat(NrUnit,  size(O.SD,1),1)+ (1:size(O.SD,1))',     repmat(i_tsk,  size(O.SD_SEM,1),1), repmat(i_BrArea,  size(O.SD_SEM,1),1) , O.sig_FR_diff, O.sig_sign , O.sig_n_bins];

            Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD = [Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD; SD]; 
            Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDP = [Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDP; SDP]; 
            Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD_SEM = [Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD_SEM; SD_SEM]; 
        end
    end
end
%%
for i_BrArea = 1: length(Out)
    
    title(['Mean_', (TargetBrainArea{i_BrArea})],'interpreter','none');
    for i_tsk = 1: numel(TaskTyp)
        O = [Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
        TaskType(i_tsk).SDmean          =   nanmean(O.SD(:, 1:end-6) - O.SDP(:, 1:end-6));
        TaskType(i_tsk).SDmean_SEM      =  nanstd(O.SD(:, 1:end-6) - O.SDP(:, 1:end-6))/ sqrt(length(TaskType(i_tsk).SDmean(:, 1:end-6) )) ;
        hold on
        if i_tsk == 1
            lineProps={'color','r','linewidth',3};
        else
            lineProps={'color','b','linewidth',3};
        end
        shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,TaskType(i_tsk).SDmean ,TaskType(i_tsk).SDmean_SEM ,lineProps,1);
        %text(BINS(10),max([Output.(target).Task.SD(u,:), Output.(target).Rest.SD(u,:)])*0.1, ['Task: trials = ' ,num2str(Output.(target).Task.NrTrials) ],'Color','red')
        %text(BINS(10),max([Output.(target).Task.SD(u,:), Output.(target).Rest.SD(u,:)]), ['Rest: trials = ' ,num2str(Output.(target).Rest.NrTrials) ],'Color','blue')

    end
end
filename= ['Average_ECG_triggered_PSTH_' (TargetBrainArea{i_BrArea})];

if savePlot; export_fig([basepath_to_save, filesep, filename], '-pdf','-transparent'); end % pdf by run
end