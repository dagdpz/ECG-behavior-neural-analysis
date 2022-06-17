function ecg_histogram(session_info,ecg_bna_cfg, Col)

load(session_info.Input_ECG);
                  
Set.min_R2R                         = 0.25; % s % needs to be changed when different monkey analyzed
T_hist_R2R_valid                    = [];
for idx = 1: length(out)
    [hist_R2R_valid,bins]   = hist(out(idx).R2R_valid ,[Set.min_R2R:0.01:1]);
    T_hist_R2R_valid = hist_R2R_valid + hist_R2R_valid;
end
     
    
   subplot(1,2,1);hold on; % subplot per BrainRegion  
   plot(bins, ig_hist2per(T_hist_R2R_valid),'k','LineWidth', 3, 'color',Col);
   % title(sprintf('%d valid R2R',length(R2R), length(R2R_valid)));
   xlabel('R2R (s)','fontsize',12 );
   ylabel('Percentage of R-peaks','fontsize',12 );
   set(gca,'xlim',[0, 1]); hold on; box on; axis square

  
%      filename= ['ECG_HistogramForEachSession'];

%     if savePlot;
%         export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
%         close all;
%     end 
