% load('Y:\Projects\Pulv_bodysignal\LFP\ECG_Magnus_complete\Magnus_Rpeak_Triggered_target_wise_Grand_grand_avg_sessions_sitesall.mat')
% %load('Y:\Projects\Pulv_bodysignal\LFP\ECG_Magnus_complete\Magnus_Rpeak_Triggered_target_wise_Grand_grand_avg_sessions_targetsall.mat')
% 
% sites=sites(ismember({sites.target},{'VPL_R','dPul_R','MD_R','VPL_L','dPul_L','MD_L'}));
% sites=ecg_bna_add_siglabels_posthoc(sites);
% 


targets={'VPL_L','MD_L','dPul_L','VPL_R','MD_R','dPul_R'};
target_combined={'VPL','MD','dPul'};

monkeys={'Magnus','Bacchus'};
for m=1:numel(monkeys)
monkey=monkeys{m};

data_path_LFP = ['Y:\Projects\Pulv_bodysignal\LFP\ECG_' monkey '_complete\Per_Site\'];
data_path_MUA = ['Y:\Projects\Pulv_bodysignal\MUA\ECG_' monkey '_complete\Per_Site\'];
%     data_path = 'Y:\Projects\Pulv_bodysignal\LFP\ECG_Magnus_TaskRest_test_generalTrig_Reref_shamim\sample session check';
% read the files
all_lfp_data = dir([data_path_LFP '*.mat']);
% storing all the Triggered parameters of the sites in the results folder:
s=0;
clear sites
for f = 1:length(all_lfp_data)
    load([data_path_LFP all_lfp_data(f).name])    
    site_LFP=ecg_bna_add_siglabels_posthoc(triggered_site_data);
    load([data_path_MUA all_lfp_data(f).name])    
    site_MUA=ecg_bna_add_siglabels_posthoc(triggered_site_data);
    if site_LFP.all_conditions_present && ismember(site_LFP.target,targets)
        s=s+1;
        for c=1:numel(site_LFP.condition)
            for e=1:numel(site_LFP.condition(c).event)
                site_LFP.condition(c).event(e).siglabels.lfp=site_LFP.condition(c).event(e).siglabels.evoked;
                site_LFP.condition(c).event(e).siglabels.mua=site_MUA.condition(c).event(e).siglabels.evoked;
            end
        end
        sites(s)=site_LFP;
    end
end


params={'lfp','mua','pow','itpc'};
e=1;
for t = 1: length(target_combined)
    sites_for_this_target=arrayfun(@(x) any(strfind(x.target,target_combined{t})),sites);
    target_sites = sites(sites_for_this_target);
    
    
    tarcon=vertcat(target_sites.condition);
    tarcon1=tarcon(:,1);
    tarcon2=tarcon(:,2);
    tarev1=vertcat(tarcon1.event);tarev1=tarev1(:,e);
    tarev2=vertcat(tarcon2.event);tarev2=tarev2(:,e);
    tarev1=[tarev1.siglabels];
    tarev2=[tarev2.siglabels];
    
    for p=1:numel(params)
        PAR=params{p};
        par1=[tarev1.(PAR)];
        par2=[tarev2.(PAR)];
        sigsign1=[par1.label];
        sigsign2=[par2.label];
        
        incongruent=(sigsign1==1 & sigsign2==-1) | (sigsign1==-1 & sigsign2==1);
        anypos = (sigsign1==1 | sigsign2==1)   & ~incongruent;
        anyneg = (sigsign1==-1 | sigsign2==-1) & ~incongruent;
        bothpos = (sigsign1==1 & sigsign2==1) ;
        bothneg = (sigsign1==-1 & sigsign2==-1) ;
        
        parout.incongruent=sum(incongruent);
        parout.anypos=sum(anypos);
        parout.anyneg=sum(anyneg);
        parout.bothpos=sum(bothpos);
        parout.bothneg=sum(bothneg);
        parout.total=numel(incongruent);
        OUT.(monkey).(target_combined{t}).(PAR)=parout;
    end
end
end

%% make a table here
FN1=fieldnames(OUT);
for m=1:numel(FN1)
    monk=FN1{m};
    FN2=fieldnames(OUT.(monk));
    for t=1:numel(FN2)
        tar=FN2{t};
        FN3=fieldnames(OUT.(monk).(tar));
        for p=1:numel(FN3)
        par=FN3{p};
        FN4=fieldnames(OUT.(monk).(tar).(par));
            for s=1:numel(FN4)
                subp=FN4{s};
                col=(m-1)*numel(FN2)+t+1;
                row=(p-1)*numel(FN4)+s+1;
                
                fulltable{1,col}=[monk '_' tar];
                fulltable{row,1}=[par '_' subp];
                fulltable{row,col}=OUT.(monk).(tar).(par).(subp);
                
                
            end
        
        
        end
    end
end

aa=1;