function ecg_bna_population_comparisons(cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is updated on 23.07.25
% based on the configurations of the data that was generated for:
% "grand_average_wo combine_hemispheres_LFP real_New GrandAvg version "
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
reprocess=0;
reprocess2=0;
normversion = 'normalized'; %'real', 'subtraction'


%target_comparisons={'dPul','VPL';'MD','VPL';'MD','dPul'};
target_comparisons=cfg.lfp.compare_targets;
condition_comparisons=cfg.lfp.compare_conditions;


targets=lower(cfg.targets);
target_combined=targets;
if cfg.combine_hemispheres
    for t=1:numel(targets)
        xx=strfind(targets{t},'_');
        if any(xx)
            target_combined{t}=targets{t}(1:xx-1);
        end
    end
end
target_combined=unique(target_combined);

monkeys=cfg.monkeys; %{'Bacchus','Magnus'};
combine_monkeys=cfg.combine_monkeys;

switch combine_monkeys
    case 0
        monkey_label='per_monkey';
    case 1
        monkey_label='monkeys_combined';
        monkeys={'all'};
end


conditions={cfg.condition.name};



%% frequencies !
startfreq=3.5;
colsbp=cool(length(cfg.lfp.frequency_bands));
fbx=find(cfg.lfp.frequency_bands(:,1)>=startfreq,1,'first');
frequency_bands=cfg.lfp.frequency_bands(fbx:end,:);
%freqName = cfg.lfp.freqName(fbx:end);
colsbp=colsbp(fbx:end,:);
freqb=num2cell(strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz'),2);
fx=find(cfg.lfp.foi>=startfreq,1,'first');
freq = cfg.lfp.foi(fx:end);

if reprocess
    for m=1:numel(monkeys)
        monkey=monkeys{m};
        
        %data_path_LFP = ['Y:\Projects\Pulv_bodysignal\LFP\ECG_' monkey '_complete\Per_Site\'];
        data_path_LFP = cfg.fldr.LFP_sites;
        
        switch combine_monkeys
            case 0
                all_lfp_data = dir([data_path_LFP filesep '*' monkey(1:3) '*.mat']);
            case 1
                all_lfp_data = dir([data_path_LFP filesep '*.mat']);
        end
        
        % storing all the Triggered parameters of the sites in the results folder:
        s=0;
        clear sites
        for f = 1:length(all_lfp_data)
            load([data_path_LFP filesep all_lfp_data(f).name])
            site_LFP=triggered_site_data;
            
            site_LFP.all_conditions_present = 1;
            
            for c= 1:numel(site_LFP.condition)
                con=site_LFP.condition(c);
                
                ntrigs=vertcat(con.event);
                ntrigs=vertcat(ntrigs.observed);
                ntrigs=vertcat(ntrigs.ntriggers);
                
                if isempty(con.event) || any(ntrigs<cfg.lfp.min_triggers_per_condition)
                    site_LFP.all_conditions_present = 0;
                    continue;
                end
            end
            
            if site_LFP.all_conditions_present && ismember(lower(site_LFP.target),targets)
                s=s+1;
                D.(monkey).sites(s)=site_LFP;
            end
        end
    end
    
    for t=1:numel(target_combined)
        target=target_combined{t};
        for m=1:numel(monkeys)
            monkey=monkeys{m};
            sites=D.(monkey).sites;
            tars=lower({sites.target});
            tar=~cellfun(@isempty,(strfind(tars,target)));
            sites=sites(tar);
            cons=vertcat(sites.condition);
            for c=1:numel(conditions)
                condition=conditions{c};
                ev=vertcat(cons(:,c).event);
                for e=1:size(ev,2)
                    norm=vertcat(ev(:,e).(normversion));
                    itpc=vertcat(norm.itpc);
                    pow=vertcat(norm.pow);
                    DAT.(monkey).(target).(condition).event(e).itpc=vertcat(itpc.mean);
                    DAT.(monkey).(target).(condition).event(e).pow=vertcat(pow.mean);
                end
            end
        end
    end
    
    save([cfg.fldr.LFP_root filesep 'Population_for_comparisons_' monkey_label],'DAT');
else
    
    load([cfg.fldr.LFP_root filesep 'Population_for_comparisons_' monkey_label],'DAT');
end

if reprocess2
    params={'pow','itpc'};
    clear comps
    for m=1:numel(monkeys)
        monkey=monkeys{m};
        for p=1:numel(params)
            par=params{p};
            %% condition_comparisons
            for t=1:numel(target_combined)
                target=target_combined{t};
                for e=1:size(cfg.analyse_states,1)
                    for c = 1:numel(condition_comparisons)
                        
                        COMP=condition_comparisons{c};
                        comps1=COMP{1};
                        comps2=COMP{2};
                        
                        
                        c1=zeros(0,0,0,0);
                        n1='';
                        for cc=1:numel(comps1)
                            nc=conditions{comps1(cc)};
                            c1(cc,:,:,:)=DAT.(monkey).(target).(nc).event(e).(par); %rest
                            n1=[n1 '+' nc];
                        end
                        c1=squeeze(mean(c1,1));
                        n1=n1(2:end);
                        
                        n2='';
                        c2=zeros(0,0,0,0);
                        for cc=1:numel(comps2)
                            nc=conditions{comps2(cc)};
                            c2(cc,:,:,:)=DAT.(monkey).(target).(conditions{comps2(cc)}).event(e).(par); %rest
                            n2=[n2 '+' nc];
                        end
                        c2=squeeze(mean(c2,1));
                        n2=n2(2:end);
                        
                        compname=[target '_' n1 '-' n2];
                        
                        [mean_diff,sig,t_map]=cluster_test(c1,c2,'paired');
                        
                        comps.(monkey).(par).conditions(c).(target).event(e).mean_diff=mean_diff;
                        comps.(monkey).(par).conditions(c).(target).event(e).sig=sig;
                        comps.(monkey).(par).conditions(c).(target).event(e).t_map=t_map;
                        comps.(monkey).(par).conditions(c).(target).event(e).name=compname;
                        comps.(monkey).(par).conditions(c).(target).event(e).nsites=num2str(size(c1,1));
                    end
                end
            end
            
            %% target_comparisons
            count=0;
            for t=1:size(target_comparisons,1)
                tar1=target_comparisons{t,1};
                tar2=target_comparisons{t,2};
                
                for c=1:numel(conditions)
                    count=count+1;
                    condition=conditions{c};
                    
                    for e=1:size(cfg.analyse_states,1)
                        
                        c1=DAT.(monkey).(tar1).(condition).(par).event(e);
                        c2=DAT.(monkey).(tar2).(condition).(par).event(e);
                        [mean_diff,sig,t_map]=cluster_test(c1,c2,'unpaired');
                        
                        compname=[tar1 '-' tar2];
                        
                        comps.(monkey).(par).nuclei(count).(condition).event(e).mean_diff=mean_diff;
                        comps.(monkey).(par).nuclei(count).(condition).event(e).sig=sig;
                        comps.(monkey).(par).nuclei(count).(condition).event(e).t_map=t_map;
                        comps.(monkey).(par).nuclei(count).(condition).event(e).name=compname;
                        comps.(monkey).(par).nuclei(count).(condition).event(e).nsites=[num2str(size(c1,1)) '/' num2str(size(c2,1))];
                    end
                    
                end
            end
        end
    end
    
    save([cfg.fldr.LFP_root filesep 'Comparisons_' monkey_label],'comps');
else
    load([cfg.fldr.LFP_root filesep 'Comparisons_' monkey_label],'comps');
end



%% append events for concatenated plotting (all analyse_states on one subplot)
% Per-event windows are concatenated with NaN gaps; t=0 is marked per segment
% (same convention as ecg_bna_plots_per_site.m).

%% plot
comparisons={'nuclei','conditions'};
params={'pow','itpc'};

for p=1:numel(params);
    par=params{p};
    for c=1:numel(comparisons)
        comp=comparisons{c};
        comp_valid=1;
        
        clear sph cb
        h = figure('units','normalized','position',[0 0 1 1]);
        
        for m=1:numel(monkeys)
            monkey=monkeys{m};
            M=monkey(1);
            if ~ismember(comp,fieldnames(comps.(monkey).(par)))
                comp_valid=0;
                continue;
            end
            
            DM=comps.(monkey).(par).(comp);
            for d=1:numel(DM)
                DMC=DM(d);
                
                switch comp
                    case 'nuclei'
                        second_levels=conditions;
                    case 'conditions'
                        second_levels=target_combined;
                end
                
                for t=1:numel(second_levels)
                    DMT=DMC.(second_levels{t});
                    [mean_diff_cat, sig_cat, tfr_events] = concat_comparison_events(DMT.event, cfg);
                    DME = DMT.event(1);
                    
                    switch comp
                        case 'nuclei'
                            nc=numel(cfg.condition);
                            nr=numel(monkeys)*numel(target_comparisons);
                            spi=(m-1)*numel(DM)+(t-1)*numel(DM)*numel(monkeys)+d;
                        case 'conditions'
                            nc=numel(condition_comparisons);
                            nr=numel(monkeys)*numel(target_combined);
                            spi=(m-1)*numel(DM)+(t-1)*numel(DM)*numel(monkeys)+d;
                    end
                    [I,J] = ind2sub([nc,nr],spi);
                    sph(J,I)=subplot(nr,nc,spi);
                    ax=sph(J,I);
                    if m==1 && t==1 && c==2
                        title(DME.name,'interpreter','none');
                    end
                    if m==1 && d==1 && c==1
                        title(second_levels{t},'interpreter','none')
                    end
                    if d==1 && c==2
                        ylabel({['Monkey ' M],[second_levels{t} '-' DME.nsites ' sites'],'Frequency (Hz)'});
                    end
                    if t==1 && c==1
                        ylabel({['Monkey ' M],[DME.name  '-' DME.nsites ' sites'],'Frequency (Hz)'});
                    end
                    
                    if m==numel(monkeys) && t==numel(second_levels)
                        xlabel('Time relative to event (s)');
                    end
                    toplot=comparison_log_plot_scale(mean_diff_cat);
                    sigplot=sig_cat;
                    
                    axis square
                    hold on
                    
                    image(toplot,'CDataMapping','scaled');
                    set(gca,'YDir','normal');
                    significance = double(sigplot);
                    ecg_bna_draw_outlines(significance==1,'r')
                    ecg_bna_draw_outlines(significance==-1,'b')
                    
                    % horizontal lines to separate frequency bands
                    fbandstart = unique(frequency_bands(:))';
                    fbandstart_idx = zeros(size(fbandstart));
                    for f = fbandstart
                        f_idx = find(abs(freq - f) == min(abs(freq - f)), 1, 'first');
                        line([tfr_events.startsamples' tfr_events.endsamples']-1/2, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                        fbandstart_idx(fbandstart == f) = f_idx;
                    end
                    
                    set(ax,'TickDir','out')
                    set(ax, 'ytick', fbandstart_idx);
                    set(ax, 'yticklabel', fbandstart);
                    set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                    add_comparison_ticks_and_labels(tfr_events, [0.5, numel(freq) + 0.5], 8);
                    set(ax, 'xlim', [0 size(toplot,2)]);
                    box on
                end
            end
        end
        
        if ~comp_valid
            close(h)
            continue;
        end
        
        %% adjusting limits
        C_lim=get(sph,'clim');
        if isa(C_lim,'cell')
            c_lim=max(abs([C_lim{:}]));
        else
            c_lim=max(abs(C_lim));
        end
        
        figure(h)
        for sp2= 1:size(sph,1)
            for sp3= 1:size(sph,2)
                ax=sph(sp2,sp3);
                if isempty(ax) || ~ishandle(ax)
                    continue
                end
                subplot(ax);
                cm=customcolormap_preset('red-yellow-blue',256);
                colormap(ax,cm);
                cb = colorbar;
                [cb_ticks, cb_labels] = comparison_log_colorbar_ticks(c_lim);
                set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                set(get(cb,'title'),'string', ['n. ' par], 'fontsize',8);
                set(ax,'clim',[-c_lim c_lim]);
                cb.Ticks = cb_ticks;
                cb.TickLabels = cb_labels;
            end
        end
        
        event_suffix = strjoin(cfg.analyse_states(:,1)', '_');
        results_file=[cfg.fldr.LFP_population filesep monkey_label '_' par '_' comp '_' event_suffix];
        wanted_size=[50 30];
        set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size]);
        export_fig(h, results_file, '-pdf');
        close(h);
    end
end

end

function [mean_diff, significance, tfr_events] = concat_comparison_events(events, cfg)
% Concatenate per-event comparison maps with NaN gaps (per-site convention).
NaNseparator = round(100/25);
sr = 1017.2526041;
sl = 1/sr*10;
tick_step = 0.2; % x-axis tick spacing (s), as in pre-concat population plots

mean_diff = [];
significance = [];
tfr_time = [];
ticksamples_tfr = [];
ticksamples_grid = [];
ticks_grid = [];

for e = 1:numel(events)
    md = squeeze(events(e).mean_diff);
    sg = squeeze(events(e).sig);

    t_before = cfg.analyse_states{e,4};
    t_after = cfg.analyse_states{e,5};
    bins_before = floor(t_before/sl);
    bins_after = ceil(t_after/sl);
    tfr_time_e = (sl*bins_before):sl:(sl*bins_after);

    segment_offset = numel(tfr_time);
    ticklabelsx = round((t_before:tick_step:t_after) * 10) / 10;
    tlstepsize = numel(tfr_time_e) / (tfr_time_e(end) - tfr_time_e(1));
    for tl = 1:numel(ticklabelsx)
        tt = ticklabelsx(tl);
        local_x = tt * tlstepsize - tfr_time_e(1) * tlstepsize;
        ticksamples_grid(end+1) = segment_offset + local_x; %#ok<AGROW>
        ticks_grid(end+1) = tt; %#ok<AGROW>
    end

    onset_s = find(tfr_time_e <= 0, 1, 'last');
    start_s = 1;
    finish_s = numel(tfr_time_e);
    ticksamples_tfr = [ticksamples_tfr, segment_offset + [start_s, onset_s, finish_s]];

    mean_diff = cat(2, mean_diff, md, nan(size(md,1), NaNseparator));
    significance = cat(2, significance, sg, nan(size(sg,1), NaNseparator));
    tfr_time = [tfr_time, tfr_time_e, nan(1, NaNseparator)];
end

tfr_events.onset = find(tfr_time == 0);
tfr_events.name = cfg.analyse_states(:,1)';
tfr_events.startsamples = ticksamples_tfr(1:3:end);
tfr_events.endsamples = ticksamples_tfr(3:3:end);
tfr_events.ticksamples = ticksamples_grid;
tfr_events.ticks = ticks_grid;
end

function add_comparison_ticks_and_labels(events, ylm, stp)
% Event onset lines, event name labels, and regular time-axis ticks.
event_onsets = events.onset;
event_names = events.name;

for so = event_onsets
    line([so so], ylm, 'color', 'k');
    event_name = event_names{event_onsets == so};
    event_name = strrep(event_name, '_', ' ');
    text(so+1, ylm(1)+stp, event_name, 'fontsize', 8, 'fontweight', 'bold');
end

set(gca, 'xtick', events.ticksamples);
set(gca, 'xticklabel', events.ticks, 'fontsize', 6);
end

function [mean_diff,significance,t_map]=cluster_test(c1,c2,pairedorunpaired)


numPermutation = 1000;
pthreshold = 0.05;
clusterthreshold = 0.05;

n1=size(c1,1);
n2=size(c2,1);

switch pairedorunpaired
    case 'paired'
        c_all=c1-c2;
    case 'unpaired'
        c_all=[c1;c2];
end

T_pos=zeros(1,numPermutation);
T_neg=zeros(1,numPermutation);

for s = 1:numPermutation % loop through shuffles
    switch pairedorunpaired
        case 'paired'
            inverseindexes = repmat((randi(2,1,size(c_all,1))-1.5)'*2,[1,size(c_all,2),size(c_all,3)]);
            realRtmp=c_all.*inverseindexes;
            [~,P]=ttest(realRtmp);
            mean_diff = mean(realRtmp, 1);
            stderr_diff = std(c_all, 0, 1) / sqrt(size(c_all,1));
            
        case 'unpaired'
            N_total=size(c_all,1);
            [~,six]=sort(randn(1,N_total));
            ix1=six(1:n1);
            ix2=six(n1+1:end);
            [~,P]=ttest2(c_all(ix1,:,:),c_all(ix2,:,:));
            mean_diff = mean(c_all(ix1,:,:), 1) - mean(c_all(ix2,:,:), 1);
            stderr_diff = sqrt(var(c_all(ix1,:,:),0,1)/n1 + var(c_all(ix2,:,:),0,1)/n2);
    end
    
    
    
    P= squeeze(P);
    
    t_map = squeeze(mean_diff ./ (stderr_diff + eps));
    
    pos_thresh = (P < pthreshold) & t_map>0 ;
    neg_thresh = (P < pthreshold) & t_map<0 ;
    
    % Find clusters using connected components
    pos_clusters = bwconncomp(pos_thresh); % i would prefer 4 here
    neg_clusters = bwconncomp(neg_thresh); % i would prefer 4 here
    
    % max sum of the tmap
    sumTpos = cellfun(@(x) sum(t_map(x)),pos_clusters.PixelIdxList);
    sumTneg = cellfun(@(x) sum(t_map(x)),neg_clusters.PixelIdxList);
    T_pos(s) = max([sumTpos 0]);
    T_neg(s) = min([sumTneg 0]);
end


%% now find significant clusters in the real data


switch pairedorunpaired
    case 'paired'
        [~,P]=ttest(c_all);
        mean_diff=mean(c1-c2,1);
        stderr_diff = std(c1-c2, 0, 1) / sqrt(size(c1,1));
        
    case 'unpaired'
        [~,P]=ttest2(c1,c2);
        mean_diff=mean(c1,1)-mean(c2,1);
        stderr_diff = sqrt(var(c1,0,1)/n1 + var(c2,0,1)/n2);
end

P= squeeze(P);
t_map = squeeze(mean_diff ./ (stderr_diff + eps));

pos_thresh = (P < pthreshold) & t_map>0 ;
neg_thresh = (P < pthreshold) & t_map<0 ;


% Find clusters using connected components
pos_clusters = bwconncomp(pos_thresh); % i would prefer 4 here
neg_clusters = bwconncomp(neg_thresh); % i would prefer 4 here

% max sum of the tmap
sumTpos = cellfun(@(x) sum(t_map(x)),pos_clusters.PixelIdxList);
sumTneg = cellfun(@(x) sum(t_map(x)),neg_clusters.PixelIdxList);


% Compute p-values using vectorized operations
% For positive clusters
p_values_pos = zeros(1, numel(sumTpos));
for i = 1:numel(sumTpos)
    p_values_pos(i) = mean(T_pos >= sumTpos(i));
end
% For negative clusters
p_values_neg = zeros(1, numel(sumTneg));
for i = 1:numel(sumTneg)
    p_values_neg(i) = mean(T_neg <= sumTneg(i));
end

sigpixpos=vertcat(pos_clusters.PixelIdxList{p_values_pos<clusterthreshold});
sigpixneg=vertcat(neg_clusters.PixelIdxList{p_values_neg<clusterthreshold});

significance=zeros(size(t_map));
significance(sigpixpos)=1;
significance(sigpixneg)=-1;

end

function y = comparison_log_plot_scale(x)
% Piecewise log display scale for normalized comparison maps.
% |x| <= 1: unchanged; |x| > 1: sign(x) * (log2(|x|) + 1)
y = x;
pos = x > 1;
neg = x < -1;
y(pos) = log2(x(pos)) + 1;
y(neg) = -(log2(abs(x(neg))) + 1);
end

function [tick_pos, tick_labels] = comparison_log_colorbar_ticks(c_lim)
% Symmetric colorbar ticks (plot scale) with pre-transform integer labels.
orig_ticks = [-1, 0, 1];
if c_lim > 1
    max_power = floor(c_lim - 1);
    pos_powers = 2.^(1:max_power);
    orig_ticks = [-fliplr(pos_powers), orig_ticks, pos_powers];
end
tick_pos = comparison_log_plot_scale(orig_ticks);
valid = abs(tick_pos) <= c_lim + 1e-6;
tick_pos = tick_pos(valid);
orig_ticks = orig_ticks(valid);
[tick_pos, sort_ix] = sort(tick_pos);
tick_labels = arrayfun(@num2str, orig_ticks(sort_ix), 'UniformOutput', false);
end
