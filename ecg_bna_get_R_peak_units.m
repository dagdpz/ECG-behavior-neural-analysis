function unit_list = ecg_bna_get_R_peak_units(cfg)

% not bad
consider_tfr=0;
freq_range=[30 50];
time_range=[-0.01 0.01];
sig_prob_thr=1;


% % not bad
% freq_range=[30 50];
% time_range=[-0.02 0.005];
% sig_prob_thr=0.5;


% freq_range=[30 50];
% time_range=[-0.01 0.01];
% sig_prob_thr=0.5;

data_path = cfg.fldr.LFP_sites;
%     data_path = 'Y:\Projects\Pulv_bodysignal\LFP\ECG_Magnus_TaskRest_test_generalTrig_Reref_shamim\sample session check';
cd(data_path)
% read the files
all_lfp_data = dir('*.mat');
% storing all the Triggered parameters of the sites in the results folder:
unit_list={};
for f = 1:length(all_lfp_data)
    load(all_lfp_data(f).name)
    site_ID = triggered_site_data.site_ID;
    addsite=1;
    for c = 1:length(triggered_site_data.condition)
        con=triggered_site_data.condition(c);
        if isempty(con.event)
            continue;
        end
        
        % identify R-peak epoch
        e=find(ismember(cfg.analyse_states(:,2),{'Rpeak'}));
        event=con.event(e);
        % think about event present?
        
        nTriggers = event.observed.ntriggers;
        time      = event.time;
        tfr_time  = event.tfr_time;
        
        freq_idx=cfg.lfp.foi>=freq_range(1) & cfg.lfp.foi<=freq_range(2);
        time_idx=tfr_time>=time_range(1) & tfr_time<=time_range(2);
        
        
        sig=event.significance.itpc;        
        sig_around_tf=squeeze(sig(1,freq_idx,time_idx));
        sigbp=event.significance.itpcbp; 
        sigbp_around_tf=squeeze(sigbp(1,5,time_idx));     
        
        mean_bins_sig=sum(sum(sig_around_tf))/numel(sig_around_tf);
        mean_bins_sigbp=sum(sum(sigbp_around_tf))/numel(sigbp_around_tf);
        
        if consider_tfr
           tocheck= mean_bins_sig>=sig_prob_thr && mean_bins_sigbp>=sig_prob_thr;
        else
           tocheck= mean_bins_sigbp>=sig_prob_thr;
        end
        
        if tocheck
            addsite=0;
        end
    end
    if addsite==1
            unit_list=[unit_list; {site_ID}];
    end
    
end
unit_list=unique(unit_list);
end