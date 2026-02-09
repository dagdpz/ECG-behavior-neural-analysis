function site=ecg_bna_add_sterr_posthoc(site)
for c= 1:numel(site.condition)
    con=site.condition(c);
    for e=1:numel(con.event)
     ev=con.event(e);       
     site.condition(c).event(e).normalized.evoked.sterr= ev.observed.evoked.sterr./ev.surrogate.evoked.std; %/sqrt(N_surrogoates);
    end
end
end