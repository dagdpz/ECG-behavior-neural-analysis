function ecg_bna_temp_siglabels

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
                    site_LFP.condition(c).event(e).max_t_after.lfp=site_LFP.condition(c).event(e).max_t_after.evoked;
                    site_LFP.condition(c).event(e).max_t_after.mua=site_MUA.condition(c).event(e).max_t_after.evoked;
                    site_LFP.condition(c).event(e).max_T_after.lfp=site_LFP.condition(c).event(e).max_T_after.evoked;
                    site_LFP.condition(c).event(e).max_T_after.mua=site_MUA.condition(c).event(e).max_T_after.evoked;
                    nsigbins_lfp=site_LFP.condition(c).event(e).nsigbins;
                    nsigbins_mua=site_MUA.condition(c).event(e).nsigbins;%% find a way to remove pre and post from nsigbins so it doesnt get confusing
                    site_LFP.condition(c).event(e).nsigbins.lfp=nsigbins_lfp;
                    site_LFP.condition(c).event(e).nsigbins.mua=nsigbins_mua;
                end
            end
            sites(s)=site_LFP;
        end
    end
    
    
    
    
    %% make a table here
    
    fulltable=[{'site'};{sites.site_ID}'];
    fulltable(:,2)=[{'target'};{sites.target}'];
    params={'lfp','mua','pow','itpc'};
    conditionnames={'Rest','Task'};
    for p=1:numel(params)
        par=params{p};
        for c=1:numel(conditionnames)
            con=conditionnames{c};
            CON=vertcat(sites.condition);
            EVE=vertcat(CON(:,c).event);
            LAB=[EVE(:,1).siglabels];
            PAR=[LAB.(par)];
            sig={PAR.label}';
            Tsum={PAR.Tsum}';
            
            namesig=[par '_' con '_' 'signed_significance'];
            nameTsum=[par '_' con '_' 'Tsum'];
            ncol=size(fulltable,2)+1;
            fulltable(:,ncol)=[{namesig};sig];
            ncol=size(fulltable,2)+1;
            fulltable(:,ncol)=[{nameTsum};Tsum];
            
        end
    end
    writetable(cell2table(fulltable),['Y:\Projects\Pulv_bodysignal\Figures\significance_table.xls' ],'Sheet',[monkey '_' 'siglabels_per_site']);
    
    
    fulltable=[{'site'};{sites.site_ID}'];
    fulltable(:,2)=[{'target'};{sites.target}'];
    params={'lfp','mua'};
    conditionnames={'Rest','Task'};
    for p=1:numel(params)
        par=params{p};
        for c=1:numel(conditionnames)
            con=conditionnames{c};
            CON=vertcat(sites.condition);
            EVE=vertcat(CON(:,c).event);
            LAB=[EVE(:,1).nsigbins];
            PAR=[LAB.(par)];
            pre={PAR.pre}';
            post={PAR.post}';
            
            namepre=[par '_' con '_' 'pre'];
            namepost=[par '_' con '_' 'post'];
            ncol=size(fulltable,2)+1;
            fulltable(:,ncol)=[{namepre};pre];
            ncol=size(fulltable,2)+1;
            fulltable(:,ncol)=[{namepost};post];
            
        end
    end
    writetable(cell2table(fulltable),['Y:\Projects\Pulv_bodysignal\Figures\significance_table.xls' ],'Sheet',[monkey '_' 'nsigbins_per_site']);
    
    
    fulltable=[{'site'};{sites.site_ID}'];
    fulltable(:,2)=[{'target'};{sites.target}'];
    params={'lfp','mua'};
    conditionnames={'Rest','Task'};
    for p=1:numel(params)
        par=params{p};
        for c=1:numel(conditionnames)
            con=conditionnames{c};
            CON=vertcat(sites.condition);
            EVE=vertcat(CON(:,c).event);
            
            LAB=[EVE(:,1).max_t_after];
            PAR=[LAB.(par)]';            
            namet=[par '_' con '_' 't_max'];
            ncol=size(fulltable,2)+1;
            fulltable(:,ncol)=[{namet};num2cell(PAR)];
            
            
            
            LAB=[EVE(:,1).max_T_after];
            PAR=[LAB.(par)]';
            nameT=[par '_' con '_' 'T_max'];            
            ncol=size(fulltable,2)+1;
            fulltable(:,ncol)=[{nameT};num2cell(PAR)];            
        end
    end
    writetable(cell2table(fulltable),['Y:\Projects\Pulv_bodysignal\Figures\significance_table.xls' ],'Sheet',[monkey '_' 'max_t']);
    
    
    
    frequency_bands = {'1-3';'4-8';'8-14';'14-30';'30-50';'50-70';'70-120'};
    fulltable=[{'site'};{sites.site_ID}'];
    fulltable(:,2)=[{'target'};{sites.target}'];
    params={'powbp','itpcbp'};
    conditionnames={'Rest','Task'};
    for p=1:numel(params)
        par=params{p};
        for c=1:numel(conditionnames)
            con=conditionnames{c};
            CON=vertcat(sites.condition);
            EVE=vertcat(CON(:,c).event);
            LAB=[EVE(:,1).meannormalized];
            
            PAR=[LAB.(par)];
            pre=num2cell(vertcat(PAR.pre));
            post=num2cell(vertcat(PAR.post));
            
            for f=1:numel(frequency_bands)
                fb=frequency_bands{f};
            namepre=[par '_' con '_' fb '_pre'];
            namepost=[par '_' con '_' fb '_post'];
            ncol=size(fulltable,2)+1;
            fulltable(:,ncol)=[{namepre};pre(:,f)];
            ncol=size(fulltable,2)+1;
            fulltable(:,ncol)=[{namepost};post(:,f)];
            end
            
        end
    end
    writetable(cell2table(fulltable),['Y:\Projects\Pulv_bodysignal\Figures\significance_table.xls' ],'Sheet',[monkey '_' 'fb_normal']);
    
    
    
    frequency_bands = {'1-3';'4-8';'8-14';'14-30';'30-50';'50-70';'70-120'};
    fulltable=[{'site'};{sites.site_ID}'];
    fulltable(:,2)=[{'target'};{sites.target}'];
    params={'powbp','itpcbp'};
    conditionnames={'Rest','Task'};
    for p=1:numel(params)
        par=params{p};
        for c=1:numel(conditionnames)
            con=conditionnames{c};
            CON=vertcat(sites.condition);
            EVE=vertcat(CON(:,c).event);            
            for f=1:numel(frequency_bands)
                fb=frequency_bands{f};
                namet=[par '_' con '_' fb '_t_max'];
                nameT=[par '_' con '_' fb '_T_max'];
                
                LAB=[EVE(:,1).max_t_after];
                PAR_all=vertcat(LAB.(par));
                PAR=PAR_all(:,f);
                ncol=size(fulltable,2)+1;
                fulltable(:,ncol)=[{namet};num2cell(PAR)];
                
                LAB=[EVE(:,1).max_T_after];
                PAR_all=vertcat(LAB.(par));
                PAR=PAR_all(:,f);
                ncol=size(fulltable,2)+1;
                fulltable(:,ncol)=[{nameT};num2cell(PAR)];
                
            end
            
        end
    end
    writetable(cell2table(fulltable),['Y:\Projects\Pulv_bodysignal\Figures\significance_table.xls' ],'Sheet',[monkey '_' 'fb_timing']);
    
    
    e=1;
    for tc=1:2
        if tc==1
            tars=targets;
            subfield='perhemisphere';
        elseif tc==2
            tars=target_combined;
            subfield='hemispherescombined';
        end
        
        for t = 1: length(tars)
            sites_for_this_target=arrayfun(@(x) any(strfind(x.target,tars{t})),sites);
            target_sites = sites(sites_for_this_target);
            
            if isempty(target_sites)
                
                for p=1:numel(params)
                    PAR=params{p};
                    parout.incongruent=0;
                    parout.onepos=0;
                    parout.oneneg=0;
                    parout.bothpos=0;
                    parout.bothneg=0;
                    parout.nosig=0;
                    parout.total=0;
                    
                    parout_percent.incongruent=0;
                    parout_percent.onepos=0;
                    parout_percent.oneneg=0;
                    parout_percent.bothpos=0;
                    parout_percent.bothneg=0;
                    parout_percent.nosig=0;
                    parout_percent.total=0;
                    
                    OUT_N.(subfield).(monkey).(tars{t}).(PAR)=parout;
                    OUT_P.(subfield).(monkey).(tars{t}).(PAR)=parout_percent;
                end
                continue;
            end
            
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
                
                nonsig=(sigsign1==0 & sigsign2==0);
                incongruent=(sigsign1==1 & sigsign2==-1) | (sigsign1==-1 & sigsign2==1);
                bothpos = (sigsign1==1 & sigsign2==1) ;
                bothneg = (sigsign1==-1 & sigsign2==-1) ;
                onepos = (sigsign1==1 | sigsign2==1)   & ~incongruent & ~bothpos;
                oneneg = (sigsign1==-1 | sigsign2==-1) & ~incongruent & ~bothneg;
                
                parout.incongruent=sum(incongruent);
                parout.onepos=sum(onepos);
                parout.oneneg=sum(oneneg);
                parout.bothpos=sum(bothpos);
                parout.bothneg=sum(bothneg);
                parout.nosig=sum(nonsig);
                parout.total=numel(incongruent);
                
                parout_percent.incongruent=round(mean(incongruent)*1000)/10;
                parout_percent.onepos=round(mean(onepos)*1000)/10;
                parout_percent.oneneg=round(mean(oneneg)*1000)/10;
                parout_percent.bothpos=round(mean(bothpos)*1000)/10;
                parout_percent.bothneg=round(mean(bothneg)*1000)/10;
                parout_percent.nosig=round(mean(nonsig)*1000)/10;
                parout_percent.total=numel(incongruent);
                
                OUT_N.(subfield).(monkey).(tars{t}).(PAR)=parout;
                OUT_P.(subfield).(monkey).(tars{t}).(PAR)=parout_percent;
            end
        end
    end
end


clear fulltable

%% make a table here
for numbers_are = {'N_','percent_'}
    N=numbers_are{:};
    switch N
        case  'N_'
            OUT=OUT_N;
        case  'percent_'
            OUT=OUT_P;
    end
    
    
    for tc=1:2
        if tc==1
            subfield='perhemisphere';
        elseif tc==2
            subfield='hemispherescombined';
        end
        clear fulltable
        FN1=fieldnames(OUT.(subfield));
        for m=1:numel(FN1)
            monk=FN1{m};
            FN2=fieldnames(OUT.(subfield).(monk));
            for t=1:numel(FN2)
                tar=FN2{t};
                FN3=fieldnames(OUT.(subfield).(monk).(tar));
                for p=1:numel(FN3)
                    par=FN3{p};
                    FN4=fieldnames(OUT.(subfield).(monk).(tar).(par));
                    for s=1:numel(FN4)
                        subp=FN4{s};
                        col=(m-1)*numel(FN2)+t+1;
                        row=(p-1)*numel(FN4)+s+1;
                        
                        fulltable{1,col}=[monk '_' tar];
                        fulltable{row,1}=[par '_' subp];
                        fulltable{row,col}=OUT.(subfield).(monk).(tar).(par).(subp);
                    end
                end
            end
        end
        
        
        writetable(cell2table(fulltable),['Y:\Projects\Pulv_bodysignal\Figures\significance_table.xls' ],'Sheet',[N '_' subfield]);
        
        aa=1;
        
    end
end

aa=1;
end