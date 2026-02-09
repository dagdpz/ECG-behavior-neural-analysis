function site_out=ecg_bna_add_siglabels_posthoc(site)

site_out.all_conditions_present = 1;


%for s=1:numel(sites)
%    site=sites(s);
for c= 1:numel(site.condition)
    con=site.condition(c);
    
    if isempty(con.event) %|| con.event(1).observed.ntriggers<1000 %% very temporary condition
        site_out.all_conditions_present = 0;
        continue;
    end
    
    for e=1:numel(con.event)
        ev=con.event(e);
        t0=find(ev.time==0);
        
        FN=fieldnames(ev.significance);
        for f=1:numel(FN)
            fn=FN{f};
            obsmean=ev.observed.(fn).mean;
            surmean=ev.surrogate.(fn).mean;
            surstd=ev.surrogate.(fn).std;
            
            T_map=(obsmean-surmean)./surstd;
            sig=ev.significance.(fn);
            
            % only consider stuff after trigger
            T_map=T_map(:,:,t0:end);
            sig=sig(:,:,t0:end);
            
            sigpos=sig==1;
            signeg=sig==-1;
            
            pos_clusters = bwconncomp(squeeze(sigpos));
            neg_clusters = bwconncomp(squeeze(signeg));
            
            sumTpos = cellfun(@(x) sum(T_map(x)),pos_clusters.PixelIdxList);
            sumTneg = cellfun(@(x) sum(T_map(x)),neg_clusters.PixelIdxList);
            
            if ~any(sumTpos) && ~any(sumTneg)
                Tsum=0;
                siglabel=0;
            elseif any(sumTpos) && any(sumTneg)
                [Tsum, ix]=max([max(-1*sumTneg) max(sumTpos) ]);
                siglabel=round((ix-1.5)*2);
            elseif any(sumTpos)
                Tsum=max(sumTpos);
                siglabel=1;
            elseif any(sumTneg)
                Tsum=max(-1*sumTneg);
                siglabel=-1;
            end
            siglabels.label=siglabel;
            siglabels.Tsum=Tsum;
            site_out.condition(c).event(e).siglabels.(fn)=siglabels;
            site_out.site_ID =site.site_ID;
            site_out.target  =site.target;            
        end
        %site.condition(c).event(e).normalized.evoked.sterr= ev.observed.evoked.sterr./ev.surrogate.evoked.std; %/sqrt(N_surrogoates);
    end
end
%end


end