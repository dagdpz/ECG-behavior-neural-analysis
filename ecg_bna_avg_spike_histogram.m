function ecg_bna_avg_spike_histogram(SPK_PSTH,session_info)
% Here comes some sort of across population plot i assume?


ECG_event=-1;
keys.PSTH_WINDOWS={'ECG',ECG_event,-0.5,0.5};
keys.PSTH_binwidth=0.01;
keys.kernel_type='gaussian';
keys.gaussian_kernel=0.02;

basepath_to_save=[session_info(1).SPK_fldr filesep 'average'];

TaskTyp = {'Rest', 'Task'};
figure;
TargetBrainArea = fieldnames(SPK_PSTH);
for i_BrArea = 1: numel(TargetBrainArea)
    
    title(['Mean_', (TargetBrainArea{i_BrArea})],'interpreter','none');
    for i_tsk = 1: numel(TaskTyp)
        O = [SPK_PSTH.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
        TaskType(i_tsk).SDmean          =   mean(O.SD - O.SDP);
        TaskType(i_tsk).SDmean_SEM      =  std(O.SD - O.SDP)/ sqrt(length(TaskType(i_tsk).SDmean )) ;
        hold on
        if i_tsk == 1
            lineProps={'color','r','linewidth',3};
        else
            lineProps={'color','b','linewidth',3};
        end
        shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,TaskType(i_tsk).SDmean ,TaskType(i_tsk).SDmean_SEM ,lineProps,1);
    end
end
filename= ['Average_ECG_triggered_PSTH_' (TargetBrainArea{i_BrArea})];

if savePlot; export_fig([basepath_to_save, filesep, filename], '-pdf','-transparent'); end % pdf by run
end