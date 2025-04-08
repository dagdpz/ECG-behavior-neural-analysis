function T=ecg_bna_SPK_criteria(SPK_PEPH,cfg)
criteria={'ID','target','SNR','Stability','FR','Ntrials','Ntrigs','AMPFRcc','AMPFRsig','FRbelow0','AMPbelow0','MI','MIpc','Amp05B2','Amp10B2','Amp20B2','Amp05B5','Amp10B5','Amp20B5'};
areas={'VPL','MD','dPul'};
for n=1:numel(criteria)
    ix.(criteria{n})=n;
end

N_a=numel(areas);
N_c=numel(cfg.condition);

for a=1:N_a
    A=areas{a};
    for c=1:N_c
        L=cfg.condition(c).name;
        T.(A).(L)=criteria;
    end
end

for a=1:N_a
    
    A=areas{a};
    a_idx=cellfun(@(x) any(strfind(x,A)),{SPK_PEPH.target});
    units=SPK_PEPH(a_idx);
    
    rows(1:N_c)=1;
    for u=1:numel(units)
        U=units(u);
        
        for c=1:N_c
            L=cfg.condition(c).name;
            criteria_field=cfg.condition(c).saccadeTask;
            if U.(L).NrTrials>0
                rows(c)=rows(c)+1;
                T.(A).(L){rows(c),ix.ID}        =U.unitId;
                T.(A).(L){rows(c),ix.target}    =U.target;
                T.(A).(L){rows(c),ix.SNR}       =U.criteria.(['SNR_' criteria_field]);
                T.(A).(L){rows(c),ix.Stability} =U.criteria.(['stability_' criteria_field]);
                
                T.(A).(L){rows(c),ix.FR}=mean(U.(L).SD);
                
                T.(A).(L){rows(c),ix.Ntrials}=U.(L).NrTrials;
                T.(A).(L){rows(c),ix.Ntrigs}=U.(L).NrTrigs;
                
                T.(A).(L){rows(c),ix.AMPFRcc}=U.(L).FR_amp_cc;
                T.(A).(L){rows(c),ix.AMPFRsig}=U.(L).FR_amp_p;
                T.(A).(L){rows(c),ix.AMPbelow0}=any(U.(L).AMP_byBin)<=0;
                T.(A).(L){rows(c),ix.FRbelow0}=any(U.(L).SD)<=0;
                T.(A).(L){rows(c),ix.MI}=max(U.(L).SDsubtractedSDP)-min(U.(L).SDsubtractedSDP);
                T.(A).(L){rows(c),ix.MIpc}=(max(U.(L).SDsubtractedSDP)-min(U.(L).SDsubtractedSDP))/mean(U.(L).SD);
                                
                Amp05=find(cumsum(U.(L).AMP_voltageBinned)>sum(U.(L).AMP_voltageBinned)*0.05,1);
                Amp10=find(cumsum(U.(L).AMP_voltageBinned)>sum(U.(L).AMP_voltageBinned)*0.1,1);
                Amp20=find(cumsum(U.(L).AMP_voltageBinned)>sum(U.(L).AMP_voltageBinned)*0.2,1);
                
                thr=U.thresholds_microV(2);
                T.(A).(L){rows(c),ix.Amp05B2}=Amp05>thr+2;
                T.(A).(L){rows(c),ix.Amp10B2}=Amp10>thr+2;
                T.(A).(L){rows(c),ix.Amp20B2}=Amp20>thr+2;
                T.(A).(L){rows(c),ix.Amp05B5}=Amp05>thr+5;
                T.(A).(L){rows(c),ix.Amp10B5}=Amp10>thr+5;
                T.(A).(L){rows(c),ix.Amp20B5}=Amp20>thr+5;
            end
        end
    end
end

% keys.tt.avg_stability                   =[2.5, Inf];       % min and max value accepted
% keys.tt.avg_single_rating               =[0, Inf];         % min and max value accepted
% keys.tt.avg_SNR                         =[4, 26];     % min and max value accepted


%%here we plot
f1=figure;
f2=figure;
for a=1:N_a    
    A=areas{a};
    for c=1:N_c
        L=cfg.condition(c).name;
        col=cfg.condition(c).color;
        TT=T.(A).(L)(2:end,:);
        
        
        O.(A).(L).enough_FR=[TT{:,ix.FR}]>=1;
        O.(A).(L).enough_stability=[TT{:,ix.Stability}]>=2.5;
        O.(A).(L).amp_above0=[TT{:,ix.AMPbelow0}]==0;
        O.(A).(L).corr_nonsig=[TT{:,ix.AMPFRsig}]>0.05;
        
        O.(A).(L).units=TT(:,ix.ID);
        O.(A).(L).enough_triggers=[TT{:,ix.Ntrigs}]>=400;
        O.(A).(L).enough_SNR=[TT{:,ix.SNR}]>=4 & [TT{:,ix.SNR}]<=26;
        
        to_use=O.(A).(L).enough_triggers & O.(A).(L).enough_SNR; 
        sigfrcc=[TT{:,ix.AMPFRsig}]<=0.05 & to_use;
        bigamp=[TT{:,ix.Amp10B5}] & to_use;
        ampfrcc=[TT{:,ix.AMPFRcc}] ;
        MI=[TT{:,ix.MI}];
        MIpc=[TT{:,ix.MIpc}];
        
               
        O.(A).(L).AMPFRcc0=sum(ampfrcc<0);
        O.(A).(L).AMPFRcc1=sum(ampfrcc<0.1);
        O.(A).(L).AMPFRcc2=sum(ampfrcc<0.2);
        
        O.(A).(L).Amp05B2=[TT{:,ix.Amp05B2}];
        O.(A).(L).Amp10B2=[TT{:,ix.Amp10B2}];
        O.(A).(L).Amp20B2=[TT{:,ix.Amp20B2}];
        O.(A).(L).Amp05B5=[TT{:,ix.Amp05B5}];
        O.(A).(L).Amp10B5=[TT{:,ix.Amp10B5}];
        O.(A).(L).Amp20B5=[TT{:,ix.Amp20B5}];
        
        O.(A).(L).excluded=O.(A).(L).units(sigfrcc & ampfrcc>0 & ~bigamp);
        
        figure(f1);
        subplot(N_a,N_c,(a-1)*N_c+c)
        hold on
        scatter(ampfrcc(to_use),MI(to_use),30,col)
        scatter(ampfrcc(sigfrcc),MI(sigfrcc),30,'k','filled')
        scatter(ampfrcc(bigamp),MI(bigamp),30,col,'filled')
        xlabel('correlation coefficient')
        ylabel('modulation index')
        title([A ' ' L]);      
        
        figure(f2);
        subplot(N_a,N_c,(a-1)*N_c+c)
        hold on
        scatter(ampfrcc(to_use),MIpc(to_use),30,col)
        scatter(ampfrcc(sigfrcc),MIpc(sigfrcc),30,'k','filled')
        scatter(ampfrcc(bigamp),MIpc(bigamp),30,col,'filled')
        xlabel('correlation coefficient')
        ylabel('modulation index (% Change)')
        title([A ' ' L]);
    end
end

end
