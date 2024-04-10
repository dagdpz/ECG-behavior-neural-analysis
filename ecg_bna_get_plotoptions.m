function EPO=ecg_bna_get_plotoptions(cfg,spokorlf)
nevents=size(cfg.analyse_states,1);

ttickpos={};
tbinpos={};
ttick={};
%tbin={};
%histbin=[];
histbinpos={};
histtickpos={};
histtick={};
Nbins_to_add=5;
added=0;
added_hist=0;

for e=1:nevents
    eventbins=(cfg.analyse_states{e,4}:cfg.(spokorlf).PSTH_binwidth:cfg.analyse_states{e,5})*1000;
    %tbin=[tbin {eventbins}];
    
    bins=eventbins-eventbins(1)+added;
    eventalignment(e)=bins(eventbins==0);
    tbinpos=[tbinpos {bins}];
    added=bins(end)+(bins(2)-bins(1))*Nbins_to_add;
    
    ttickpos=[ttickpos {[bins(1) eventalignment(e) bins(end)]}];
    ttick   =[ttick    {[eventbins(1) 0 eventbins(end)]}];
    
    histbins=cfg.(spokorlf).histbins;
    %histbin=[histbin histbins];
    bins=histbins-histbins(1)+added_hist;
    histbinpos=[histbinpos {bins}];
    histtickpos=[histtickpos {[bins(1) bins(end)]}];
    histtick=[histtick {[histbins(1) histbins(end)]}];
    added_hist=bins(end)+(bins(2)-bins(1))*Nbins_to_add;
end
EPO.xticks=[ttickpos{:}];
EPO.xtickslabels=[ttick{:}];
EPO.eventalignment=eventalignment;
EPO.tbinpos=tbinpos;
EPO.histbinpos=histbinpos;
EPO.histticks=[histtickpos{:}];
EPO.histtickslabels=[histtick{:}];

% EPO.tbin=tbin;
% EPO.histbin=histbin;
% EPO.histtickpos=histtickpos;


% EPO.ttickpos=ttickpos;
% EPO.ttick=ttick;
% EPO.tbin=tbin;
% EPO.histbin=histbin;
% EPO.histtickpos=histtickpos;
end