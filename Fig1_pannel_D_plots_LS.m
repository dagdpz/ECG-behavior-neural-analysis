clc
%%
folder_load_ECG='Y:\Projects\Pulv_bodysignal\ECG\ECG_Bacchus_TaskRest_generalTrig_ECG\Per_Site_smpl\';
folder_load_MUA='Y:\Projects\Pulv_bodysignal\MUA\ECG_Bacchus_Complete\Per_Site\';
folder_load_LFP='Y:\Projects\Pulv_bodysignal\LFP\ECG_Bacchus_Complete\Per_Site\';
sessions = {'20210720','20210730', '20220211'};
sites_per_session = {3,8,4}; 




h = figure('units','normalized','position',[0 0 1 1]);

for ses = 1:length(sessions)
    session = sessions{ses};
    Site_ID = ['Bac_',session,'_Site_01'];
    %%
    load([folder_load_ECG,Site_ID,'.mat'])
    
    ecg_task = triggered_site_data.condition(2).event(1).real.ecg  ;
    ecg_rest = triggered_site_data.condition(1).event(1).real.ecg  ;
    ecgtime = triggered_site_data.condition(1).event(1).time  ;
    %     PlotMode = {'real','normalized'};
    PlotMode = {'observed'};
    for Pm = 1:length(PlotMode)
        
        allSesData_mua = dir([folder_load_MUA,'*_',session,'*.mat']);
        cd(folder_load_MUA)
        s=    sites_per_session{ses};
            
        load(allSesData_mua(s).name)
        
        triggered_site_data=ecg_bna_add_sterr_posthoc(triggered_site_data);
        if isempty(triggered_site_data.condition(1).event) || isempty(triggered_site_data.condition(2).event)
            disp([sessions{ses} 'site' num2str(s)])
        end
        mua_rest = triggered_site_data.condition(1).event(1).(PlotMode{Pm}).evoked;
        mua_task = triggered_site_data.condition(2).event(1).(PlotMode{Pm}).evoked;
        
        target=triggered_site_data.target;
        
        allSesData_lfp = dir([folder_load_LFP,'*_',session,'*.mat']);
        cd(folder_load_LFP)
        
        load(allSesData_lfp(s).name)
        
        triggered_site_data=ecg_bna_add_sterr_posthoc(triggered_site_data);
        lfp_rest = triggered_site_data.condition(1).event(1).(PlotMode{Pm}).evoked;
        lfp_task = triggered_site_data.condition(2).event(1).(PlotMode{Pm}).evoked;
        time = triggered_site_data.condition(1).event(1).time  ;
        
        clear triggered_site_data
        
        % Smoothing Kernel here:
        cfg.lfp.smoothWin = 5;
        win = 1:cfg.lfp.smoothWin; win=win-(numel(win)+1)/2;
        half_win = ceil(size(win,2)/2)-1;
        gaussian_kernel=normpdf(win,0,numel(win)/6);
        gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);
        
        subplot(3,3,ses) % plotting ECG   
        title([target(1:end-2) '-' allSesData_lfp(s).name(5:end-4)], 'interpreter','none')
        linePropsR={'color',[0 0 1]};
        linePropsT={'color',[1 0 0]};
        % lineProps={'color',[0.1 0.1 0.1]};
        meanecg_T = squeeze(ecg_task.mean)'*1000;
        %stdecgT = squeeze(ecg_task.std)';
        sterrecgT = squeeze(ecg_task.sterr)'*1000;
        
        meanecg_R = squeeze(ecg_rest.mean)'*1000;
        %stdecgR = squeeze(ecg_rest.std)';
        sterrecgR = squeeze(ecg_rest.sterr)'*1000;
        
        % perc5 = squeeze(ecg.conf95(1,:,:))';
        % perc95 = squeeze(ecg.conf95(2,:,:))';
        shadedErrorBar(ecgtime,meanecg_R,sterrecgR,linePropsR,1);
        shadedErrorBar(ecgtime,meanecg_T,sterrecgT,linePropsT,1);
        
        xlabel('Time(s)'); ylabel(' ECG (mV)');
        ylim([-2 3])
        xlim([time(1) time(end)])
        axis square
        line([0 0], ylim, 'color', 'k');
        
        
        subplot(3,3,ses+3) % plotting LFP
        
        linePropsR={'color',[0 0 1]};
        linePropsT={'color',[1 0 0]};
        meanR = squeeze(lfp_rest.mean)';
        sterrR = squeeze(lfp_rest.sterr)';        
        meanT = squeeze(lfp_task.mean)';
        sterrT = squeeze(lfp_task.sterr)';
        shadedErrorBar(time,meanR,sterrR,linePropsR,1); hold on
        shadedErrorBar(time,meanT,sterrT,linePropsT,1);
        
        xlabel('Time(s)');
        switch PlotMode{Pm}
            case 'real'
                ylabel(' LFP (mV)');
            case 'normalized'
                ylabel(' normalized LFP evoked');
        end
        switch ses
            case 1
        ylim([-6 6])
            case 2
        ylim([-6 11])
            case 3
        ylim([-50 30])
        end
        xlim([time(1) time(end)])
        axis square
        line([0 0], ylim, 'color', 'k');
        
        
        subplot(3,3,ses+6) % plotting MUA        
        linePropsR={'color',[0 0 1]};
        linePropsT={'color',[1 0 0]};
        
        meanR=squeeze(mua_rest.mean);
        meanT=squeeze(mua_task.mean);
                
        sterrR = squeeze(mua_rest.sterr);
        sterrT = squeeze(mua_task.sterr);
        shadedErrorBar(time,meanR,sterrR,linePropsR,1); hold on
        shadedErrorBar(time,meanT,sterrT,linePropsT,1);
                
        xlabel('Time(s)');
        ylabel(' MUA (z-score)');
        ylim([-0.15 0.15])
        xlim([time(1) time(end)])
        axis square
        line([0 0], ylim, 'color', 'k');
    end
    
end
        
figname = fullfile(['Y:\Projects\Pulv_bodysignal\Figures',filesep,'ecg_lfp_mua_allsites','_',PlotMode{Pm}]);
export_fig(h,[figname,'.pdf']);