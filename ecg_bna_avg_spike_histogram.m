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
    for i_brain = 1 : length(BrainArea)
        if ~ismember(TargetBrainArea, BrainArea(i_brain))
            TargetBrainArea = [TargetBrainArea, BrainArea(i_brain)];
        end
    end
end

%% Create function to concatenate the variables per TargetBrainArea

% initate empty structure
for i_BrArea = 1: numel(TargetBrainArea)
    for i_tsk = 1: numel(TaskTyp)
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).CmbDataSig_SD  = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD             = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD_SEM         = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDP            = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_period     = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR             = [];
        Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).quantSNR       = [];
        
    end
end


for i = 1: length(SPK_PSTH)
    BrainArea = fieldnames(SPK_PSTH{i});
    for i_BrArea = 1: numel(BrainArea)
        for i_tsk = 1: numel(TaskTyp)
            O = [SPK_PSTH{1,i}.(BrainArea{i_BrArea}).(TaskTyp{i_tsk})];
            
             O.SDP              =   O.SDP( ~isnan(O.SDP(:,end)),:) ; 
             O.SD_SEM           =   O.SD_SEM( ~isnan(O.SD_SEM(:,end)),:) ; 
             O.sig_FR_diff      =   O.sig_FR_diff(~isnan(O.SD(:,end))) ; 
             O.sig_sign         =   O.sig_sign(~isnan(O.SD(:,end))) ; 
             O.sig_n_bins       =   O.sig_n_bins(~isnan(O.SD(:,end))) ; 
             O.sig_time         =   O.sig_time(~isnan(O.SD(:,end))) ; 
             O.sig_period       =   O.sig(~isnan(O.SD(:,end)),:) ;  % significance in a given period
             O.sig_all          =   O.sig_all(~isnan(O.SD(:,end)),:) ; 
             O.FR               =   O.FR(~isnan(O.SD(:,end)),:) ; 
             O.quantSNR         =   O.quantSNR(~isnan(O.SD(:,end)),:) ; 

             O.SD               =   O.SD( ~isnan(O.SD(:,end)),:) ; 

            NrUnit = size(O.SD, 1) ;
            
            % remove NAN % zeros
            CmbDataSig_SD      = [O.SD,    repmat(NrUnit,  size(O.SD,1),1)+ (1:size(O.SD,1))',    repmat(i_tsk,  size(O.SD,1),1) ,    repmat(i_BrArea,  size(O.SD,1),1) , O.sig_FR_diff, O.sig_sign , O.sig_n_bins, O.sig_time ];

            Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).CmbDataSig_SD = [Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).CmbDataSig_SD; CmbDataSig_SD];
            Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD = [Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD; O.SD];
            Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDP = [Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDP; O.SDP];
            Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD_SEM = [Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).SD_SEM; O.SD_SEM];
            Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_period = [Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).sig_period; O.sig_period];
            Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).quantSNR = [Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).quantSNR; O.quantSNR];
            Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR = [Out.(BrainArea{i_BrArea}).(TaskTyp{i_tsk}).FR; O.FR];

        end
    end
end
%%
for i_BrArea = 1: length(Out)
    hf = figure('Name',sprintf(TargetBrainArea{i_BrArea}),'Position',[200 100 1400 1200],'PaperPositionMode', 'auto');
    
    
    for i_tsk = 1: numel(TaskTyp)
        out = [Out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk})];
        % mean(SD - SDP)
        out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean          =   nanmean(out.SD - out.SDP);
        % standard error of the SDmean
        out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SEM      =  nanstd(out.SD(:, 1:end) - out.SDP(:, 1:end))/ sqrt(length(out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean(:, 1:end) )) ;
        
        %% Calculate the amount of modulation in FR from the significant interval 
        Signal =    out.SD - out.SDP;   out.Modulation = NaN(1,size(Signal,1));
        for i = 1: size(Signal,1)    
            if any(out.sig_period(i,:))
               out.Modulation(i) = max(abs(Signal(logical(out.sig_period(i,:))))) +  abs(min(Signal(i,:)));
            end
        end        
        out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).PcSig_FR_diff = [sum(out.CmbDataSig_SD(:, end-2)== -1)/ length(out.CmbDataSig_SD(:, end)), sum(out.CmbDataSig_SD(:, end-2)== 1)/ length(out.CmbDataSig_SD(:, end)), sum(out.CmbDataSig_SD(:, end-2)== 0)/ length(out.CmbDataSig_SD(:, end))];
        out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).PcTtoPeak = [sum(out.CmbDataSig_SD(:, end)< 0)/ length(out.CmbDataSig_SD(:, end)), sum(out.CmbDataSig_SD(:, end)> 0)/ length(out.CmbDataSig_SD(:, end)), sum(out.CmbDataSig_SD(:, end-2)== 0)/ length(out.CmbDataSig_SD(:, end))];
        
        
        Sign = -1;
        out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SigDiffFR_Neg    =   (out.SD(out.CmbDataSig_SD(:, end-2)== Sign, 1:end) - out.SDP(out.CmbDataSig_SD(:, end-2)== Sign, 1:end));
        out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SEM_SigDiffFR_Neg      =  (out.SD(out.CmbDataSig_SD(:, end-2)== Sign, 1:end) - out.SDP(out.CmbDataSig_SD(:, end-2)== Sign, 1:end))/ sqrt(length(out.CmbDataSig_SD(out.CmbDataSig_SD(:, end-2)== Sign, end )));
        %!!! CHECK THE DEVISOR !!!
        Sign = 1;
        out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SigDiffFR_Pos    =   (out.SD(out.CmbDataSig_SD(:, end-2)== Sign, 1:end) - out.SDP(out.CmbDataSig_SD(:, end-2)== Sign, 1:end));
        out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SEM_SigDiffFR_Pos      =  (out.SD(out.CmbDataSig_SD(:, end-2)== Sign, 1:end) - out.SDP(out.CmbDataSig_SD(:, end-2)== Sign, 1:end))/ sqrt(length(out.CmbDataSig_SD(out.CmbDataSig_SD(:, end-2)== Sign, end ))) ;
        
        
        
        hold on
        ha1 = subplot(2,4,[1:2]); %
        if i_tsk == 1
            lineProps={'color','r','linewidth',3};
        else
            lineProps={'color','b','linewidth',3};
        end
        
        shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean ,out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SEM ,lineProps,1);
       title(['Population:  (all)' (TargetBrainArea{i_BrArea}) ' units'],'interpreter','none');

        
        
        % bar plot how many are significant & positive and negative FR?
        
        if i_tsk == 1
            text(-400,-1, ['Rest: units = ' ,num2str(max(out.CmbDataSig_SD(:, end-6))) ],'Color','red')
            
            for i = 1: size(out.SD(:, 1:end),1)
                line((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000 , (out.SD(i, 1:end) - out.SDP(i, 1:end)), 'color',[1 0 0 0.3]);
            end
            ha1 = subplot(2,4,3);% pie plot how many
        else
            text(-400,-1.5, ['Task: units = ' ,num2str(max(out.CmbDataSig_SD(:, end-6))) ],'Color','blue')
            
            for i = 1: size(out.SD(:, 1:end),1)
                line((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000 , out.SD(i, 1:end) - out.SDP(i, 1:end), 'color',[0 0 1 0.3] ,'LineWidth', 1);
            end
            ha1 = subplot(2,4,4);
            
        end
        %Bar graph
        
        barpairs =  [out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).PcSig_FR_diff; out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).PcTtoPeak];
        b = bar(barpairs,'stacked', 'Facecolor','flat' ); 
        title([(TaskTyp(i_tsk))],'interpreter','none');
        b(3).FaceColor = [0.5 .5 .5];
        legend({'sig','sig', 'non sign'},'location','Best');
        set(gca,'XTickLabel',{'diff FR (-/+) ' 'Before/After Rpeak'},'fontsize',10);


        % display only significant units showing a increase in FR
        ha1 = subplot(2,4,[5]); %
        if i_tsk == 1
            lineProps={'color','r','linewidth',3};
             text(-400,1, ['Rest: units = ' ,num2str(numel(out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SigDiffFR_Pos(:, end))) ],'Color','red')

        else
            lineProps={'color','b','linewidth',3};
            text(-400,1.5, ['Task: units = ' ,num2str(numel(out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SigDiffFR_Pos(:, end))) ],'Color','blue')

        end
        shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SigDiffFR_Pos) , nanstd(out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SEM_SigDiffFR_Pos) ,lineProps,1);
        set(gca,'ylim',[-2.5, 3]);
        title(['units showing a sig. INCREASE in FR'],'interpreter','none');

       % display only significant units showing a decrease in FR
       ha1 = subplot(2,4,[6]); %
        if i_tsk == 1
            lineProps={'color','r','linewidth',3};
            text(-400,1, ['Rest: units = ' ,num2str(numel(out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SigDiffFR_Neg(:, 1))) ],'Color','red')
        Color = [1 0 0]; 
        else
            lineProps={'color','b','linewidth',3};
            text(-400,1.5, ['Task: units = ' ,num2str(numel(out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SigDiffFR_Neg(:, 1))) ],'Color','blue')
        Color = [0 0 1]; 
        end
        
        shadedErrorBar((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,nanmean(out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SigDiffFR_Neg) , nanstd(out.(TargetBrainArea{i_BrArea}).(TaskTyp{i_tsk}).SDmean_SEM_SigDiffFR_Neg) ,lineProps,1);
        set(gca,'ylim',[-2.5, 3]);
        title(['units showing a sig. DECREASE in FR'],'interpreter','none');
        
        
     ha1 = subplot(2,4,[7]); %
     scatter(out.quantSNR , out.Modulation, 'filled', 'MarkerFaceColor',Color)
     ylabel('Amount of Modulation (spike/s)','fontsize',14,'fontweight','b' );
     xlabel('Signal-to-Noise','fontsize',14,'fontweight','b' );
     
     ha1 = subplot(2,4,[8]); %
     scatter(out.FR , out.Modulation, 'filled', 'MarkerFaceColor',Color)
     ylabel('Amount of Modulation (spike/s)','fontsize',14,'fontweight','b' );
     xlabel('mean firing rate','fontsize',14,'fontweight','b' );


    end
    filename= ['Average_ECG_triggered_spike_' (TargetBrainArea{i_BrArea})];
    
if savePlot; 
      print(gcf,['Y:\Projects\Pulv_distractor_spatial_choice\ECG_triggered_spikes\ver1\per_unit\Population' filesep, filename,'.pdf'],'-dpdf','-r0');
   % export_fig(['Y:\Projects\Pulv_distractor_spatial_choice\ECG_triggered_spikes\ver1\per_unit\Population' filesep ], '-pdf','-transparent'); 

   % export_fig([basepath_to_save, filesep,'Average', filesep, filename], '-pdf','-transparent'); end % pdf by run

end

end