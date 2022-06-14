function bsa_histogram(sessions_info(i))
for v = 1:length(versions)
    version = versions{v};

    ecg_bna_cfg = ecg_bna_define_settings(project,version);

    %% Get info about sessions to be analysed

    % Read the info about sessions to analyse
    sessions_info = ecg_bna_cfg.session_info;

  for i = 1:length(session_info)
      
              session_name = [session_info(i).Monkey '_' session_info(i).Date];
              load(session_info.Input_ECG);
              
                  
                  Set.min_R2R                         = 0.25; % s
                  T_hist_R2R_valid = []; 
                  for idx = 1: length(out)
                  [hist_R2R_valid,bins]   = hist(out(idx).R2R_valid ,[Set.min_R2R:0.01:1]);
                  T_hist_R2R_valid = hist_R2R_valid + hist_R2R_valid; 
                  end
                  
    subplot(1,2,1);hold on; 
    %plot(nanmean([out.mean_R2R_valid_bpm]),0,'mv','MarkerSize',6);hold on; 
    plot(bins,T_hist_R2R_valid,'k');
   % title(sprintf('%d valid R2R',length(R2R), length(R2R_valid)));
    xlabel('R2R (s)');
    ylabel('count of R-peaks');
  end
  
      filename= ['ECG_HistogramForEachSession'];

    if savePlot;
        export_fig([basepath_to_save,filesep ,filename], '-pdf'); %,'-transparent'
        close all;
    end 
end
