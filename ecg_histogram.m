function ecg_histogram(session_info,ecg_bna_cfg, Col)

load(session_info.Input_ECG);
                  
Set.min_R2R                         = 0.25; % s % needs to be changed when different monkey analyzed
T_hist_R2R_valid                    = [];
for idx = 1: length(out)
    [hist_R2R_valid,bins]   = hist(out(idx).R2R_valid ,[Set.min_R2R:0.01:1]);
    T_hist_R2R_valid = hist_R2R_valid + hist_R2R_valid;
    
    if max(out(idx).R2R_valid) > 0.7
    Interval = []; 
    Interval =       (find(out(idx).R2R_valid == max(out(idx).R2R_valid))-5) :1: (find(out(idx).R2R_valid == max(out(idx).R2R_valid))+5  ); 

    disp([idx, { session_info.Input_ECG(end - 33: end) } , max(out(idx).R2R_valid) ,out(idx).R2R_t(find(out(idx).R2R_valid == max(out(idx).R2R_valid)))])
    disp(out(idx).R2R_valid(Interval))
    end
    
end
     

   subplot(1,2,2);hold on; % subplot per BrainRegion  
   plot(bins, ig_hist2per(T_hist_R2R_valid),'k','LineWidth', 3, 'color',Col);
   % title(sprintf('%d valid R2R',length(R2R), length(R2R_valid)));
   xlabel('R2R (s)','fontsize',12 );
   ylabel('Percentage of R-peaks','fontsize',12 );
   set(gca,'xlim',[0  , 1]); hold on; box on; axis square

   subplot(1,2,1);hold on; % subplot per BrainRegion  
   plot(bins, ig_hist2per(T_hist_R2R_valid),'k','LineWidth', 3, 'color',Col);
   % title(sprintf('%d valid R2R',length(R2R), length(R2R_valid)));
   xlabel('R2R (s)','fontsize',12 );
   ylabel('Percentage of R-peaks','fontsize',12 );
   set(gca,'xlim',[Set.min_R2R   , 1]); hold on; box on; axis square

   
   
  
%      filename= ['rightPul_ECG_HistogramForEachSession'];

%     if savePlot;
%       basepath_to_save=[session_info(1).SPK_fldr filesep 'Population'];

%         export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
%         close all;
%     end 

