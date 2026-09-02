function site_out=ecg_bna_add_siglabels_posthoc(site)

site_out.all_conditions_present = 1;

bins1=16:26;
bins2=37:47;
for c= 1:numel(site.condition)
    con=site.condition(c);
    
    if isempty(con.event) %|| con.event(1).observed.ntriggers<1000 %% very temporary condition
        site_out.all_conditions_present = 0;
        continue;
    end
    
    for e=1:numel(con.event)
        ev=con.event(e);
        t0=find(ev.time==0);
        
        nsigbins.pre=sum(abs(ev.significance.evoked(bins1)));
        nsigbins.post=sum(abs(ev.significance.evoked(bins2)));
        site_out.condition(c).event(e).nsigbins=nsigbins;
        
        FN=fieldnames(ev.significance);
        for f=1:numel(FN)
            fn=FN{f};
            obsmean=ev.observed.(fn).mean;
            
            
            meannormalized.pre=mean(ev.normalized.(fn).mean(:,:,bins1),3);
            meannormalized.post=mean(ev.normalized.(fn).mean(:,:,bins2),3);
            
            surmean=ev.surrogate.(fn).mean;
            surstd=ev.surrogate.(fn).std;
            
            T_map=(obsmean-surmean)./surstd;
            sig=ev.significance.(fn);
            
            % only consider stuff after trigger
            T_map=T_map(:,:,t0:end);
            sig=sig(:,:,t0:end);
            
            sigpos=sig==1;
            signeg=sig==-1;
            sig=logical(sig);
            
            pos_clusters = bwconncomp(squeeze(sigpos));
            neg_clusters = bwconncomp(squeeze(signeg));
            sig_clusters = bwconncomp(squeeze(sig));
            
            sumTpos = cellfun(@(x) sum(T_map(x)),pos_clusters.PixelIdxList);
            sumTneg = cellfun(@(x) sum(T_map(x)),neg_clusters.PixelIdxList);
            sumT = cellfun(@(x) sum(abs(T_map(x))),sig_clusters.PixelIdxList);
            
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
            
            if ismember(fn,{'itpcbp','powbp'})
                max_t=[];
                max_T=[];
                for fb=1:size(sig,2)
                    
                    sig_clusters = bwconncomp(squeeze(sig(:,fb,:)));
                    T_map_fb=squeeze(T_map(:,fb,:));
                    sumT = cellfun(@(x) sum(abs(T_map_fb(x))),sig_clusters.PixelIdxList);
                    if any(sumT)
                        [Tsum_abs, ix]=max(sumT);
                        pix=sig_clusters.PixelIdxList{ix};
                        T_map_abs=abs(T_map_fb);
                        T_map_abs(setdiff(1:numel(T_map_fb),pix))=0;
                        [max_T_fb, max_t_bin]=max(T_map_abs);
                        max_t_fb=(max_t_bin-1)*10;
                        max_t=[max_t max_t_fb];
                        max_T=[max_T max_T_fb];
                    else
                        max_t=[max_t NaN];
                        max_T=[max_T NaN];
                    end
                end
                
                aaaa=1;
            else
                if any(sumTpos) || any(sumTneg)
                    % exclude points NOT from the largest clusters
                    [Tsum_abs, ix]=max(sumT);
                    pix=sig_clusters.PixelIdxList{ix};
                    T_map_abs=abs(T_map);
                    T_map_abs(setdiff(1:numel(T_map),pix))=0;
                    [max_T, max_t_idx]=max(T_map_abs);
                    [~,~,max_t_bin]=ind2sub(size(T_map_abs),max_t_idx);
                    max_t=(max_t_bin-1)*10;
                else
                    max_T=NaN;
                    max_t=NaN;
                end
            end
            
            siglabels.label=siglabel;
            siglabels.Tsum=Tsum;
            site_out.condition(c).event(e).siglabels.(fn)=siglabels;
            site_out.condition(c).event(e).meannormalized.(fn)=meannormalized;
            site_out.condition(c).event(e).max_t_after.(fn)=max_t;
            site_out.condition(c).event(e).max_T_after.(fn)=max_T;
            site_out.site_ID =site.site_ID;
            site_out.target  =site.target;            
        end
        %site.condition(c).event(e).normalized.evoked.sterr= ev.observed.evoked.sterr./ev.surrogate.evoked.std; %/sqrt(N_surrogoates);
    end
end
%end


end