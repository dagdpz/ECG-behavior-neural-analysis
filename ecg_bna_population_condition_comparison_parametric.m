function ecg_bna_population_condition_comparison_parametric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is updated on 23.07.25
% based on the configurations of the data that was generated for:
% "grand_average_wo combine_hemispheres_LFP real_New GrandAvg version "
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
normversion = 'subtraction'; %'real', 'subtraction'

targets={'VPL_L','MD_L','dPul_L','VPL_R','MD_R','dPul_R'};
target_combined={'VPL','MD','dPul'};

monkeys={'Bacchus','Magnus'};

conditions={'Rest','Task'};



%% frequencies !

cfg.lfp.foi             = logspace(log10(4), log10(120), 60); %% add something to select frequency to start plotting...
cfg.lfp.frequency_bands = [1 3; 4 8; 8 14; 14 30; 30 50; 50 70; 70 120];

startfreq=3.5;
colsbp=cool(length(cfg.lfp.frequency_bands));
fbx=find(cfg.lfp.frequency_bands(:,1)>=startfreq,1,'first');
frequency_bands=cfg.lfp.frequency_bands(fbx:end,:);
%freqName = cfg.lfp.freqName(fbx:end);
colsbp=colsbp(fbx:end,:);
freqb=num2cell(strcat(num2str(round(frequency_bands(:,1))), '-',num2str(round(frequency_bands(:,2))), ' Hz'),2);
fx=find(cfg.lfp.foi>=startfreq,1,'first');
freq = cfg.lfp.foi(fx:end);



sr=1017.2526041; %hz
sl=1/sr*10;
lfp_time=[-sl*26:sl:sl*26];
ticklabelsx=[-0.2,0,0.2];
for tl=1:numel(ticklabelsx)
    tt=ticklabelsx(tl);
    tlstepsize=numel(lfp_time)/(lfp_time(end)-lfp_time(1));
    ticksx(tl)=tt*tlstepsize-lfp_time(1)*tlstepsize;
end

e=1;

for m=1:numel(monkeys)
    monkey=monkeys{m};
    
    data_path_LFP = ['Y:\Projects\Pulv_bodysignal\LFP\ECG_' monkey '_complete\Per_Site\'];
    all_lfp_data = dir([data_path_LFP '*.mat']);
    % storing all the Triggered parameters of the sites in the results folder:
    s=0;
    clear sites
    for f = 1:length(all_lfp_data)
        load([data_path_LFP all_lfp_data(f).name])
        site_LFP=triggered_site_data;
        
        site_LFP.all_conditions_present = 1;

        for c= 1:numel(site_LFP.condition)
            con=site_LFP.condition(c);
            if isempty(con.event) %|| con.event(1).observed.ntriggers<1000 %% very temporary condition
                site_LFP.all_conditions_present = 0;
                continue;
            end
        end
        
        if site_LFP.all_conditions_present && ismember(site_LFP.target,targets)
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
        tars={sites.target};
        tar=~cellfun(@isempty,(strfind(tars,target)));
        sites=sites(tar);
        cons=vertcat(sites.condition);
        for c=1:numel(conditions)
            condition=conditions{c};
            ev=vertcat(cons(:,c).event);
            norm=vertcat(ev(:,e).normalized);
            itpc=vertcat(norm.itpc);
            pow=vertcat(norm.pow);
            DAT.(monkey).(target).(condition).itpc=vertcat(itpc.mean);
            DAT.(monkey).(target).(condition).pow=vertcat(pow.mean);
        end
    end
end

save('Y:\Projects\Pulv_bodysignal\Figures\data_reordered','DAT');
target_comparisons={'dPul','VPL';'MD','VPL';'MD','dPul'};
params={'pow','itpc'};
clear comps
for m=1:numel(monkeys)
    monkey=monkeys{m};
    for p=1:numel(params)
        par=params{p};
        %% condition_comparisons
        for t=1:numel(target_combined)
            target=target_combined{t};
            
            
            c2=DAT.(monkey).(target).(conditions{1}).(par); %rest
            c1=DAT.(monkey).(target).(conditions{2}).(par); %task
            
            [mean_diff,sig,t_map]=cluster_test(c1,c2,'paired');
            
            comps.(monkey).(par).conditions(t).mean_diff=mean_diff;
            comps.(monkey).(par).conditions(t).sig=sig;
            comps.(monkey).(par).conditions(t).t_map=t_map;
            comps.(monkey).(par).conditions(t).name=target;
        end
        
        %% target_comparisons
        count=0;
        for t=1:size(target_comparisons,1)
            tar1=target_comparisons{t,1};
            tar2=target_comparisons{t,2};
            
            for c=1:numel(conditions)
                count=count+1;
                condition=conditions{c};
                
                c1=DAT.(monkey).(tar1).(condition).(par);
                c2=DAT.(monkey).(tar2).(condition).(par);
                [mean_diff,sig,t_map]=cluster_test(c1,c2,'unpaired');
                
                compname=[tar1 '-' tar2 ',' condition];
                
                comps.(monkey).(par).nuclei(count).mean_diff=mean_diff;
                comps.(monkey).(par).nuclei(count).sig=sig;
                comps.(monkey).(par).nuclei(count).t_map=t_map;
                comps.(monkey).(par).nuclei(count).name=compname;
                
            end
        end
    end
end

save('Y:\Projects\Pulv_bodysignal\Figures\nucleus_and_condition_comparisons','comps');


%% plot
comparisons={'nuclei','conditions'};
params={'pow','itpc'};

for p=1:numel(params);
    par=params{p};
    for c=1:numel(comparisons)
        comp=comparisons{c};
        h = figure('units','normalized','position',[0 0 1 1]);
        clear sph
        for m=1:numel(monkeys)
            monkey=monkeys{m};
            M=monkey(1);
            DM=comps.(monkey).(par).(comp);
            for d=1:numel(DM)
                
                switch comp
                    case 'nuclei'
                        spi=(ceil(d/2)-1)*2+2*(m-1)+d;
                        [I,J] = ind2sub([2,6],spi);
                        sph(J,I)=subplot(6,2,spi);
                        ylabel({['Monkey ' M],[DM(d).name],'Frequency (Hz)'});
                        ax=sph(J,I);
                    case 'conditions'
                        spi=(d-1)*2 + m;
                        sph(spi,1)=subplot(6,1,spi);
                        ylabel({['Monkey ' M],[DM(d).name],'Frequency (Hz)'});
                        ax=sph(spi,1);
                end
                if m==1 && d==1 && c==2
                    title(['Task-Rest']);
                end
                if m==1 && d==1 && c==1
                    title(conditions(1))
                end
                if m==1 && d==2 && c==1
                    title(conditions(2))
                end
                if m==numel(monkeys) && t==numel(targets)
                    xlabel('Time relative to R peak (s)');
                end
                toplot=squeeze(DM(d).mean_diff);
                sigplot=DM(d).sig;
                
                xlim([0 size(toplot,2)]);
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
                    line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                    fbandstart_idx(fbandstart == f) = f_idx;
                end
                line([ticksx(ticklabelsx==0) ticksx(ticklabelsx==0)], ylim, 'color', 'k');
                
                set(ax,'TickDir','out')
                set(ax, 'ytick', fbandstart_idx);
                set(ax, 'yticklabel', fbandstart);
                set(ax, 'xtick', ticksx);
                set(ax, 'xticklabel', ticklabelsx);
                
                
                set(ax, 'ylim', [0.5,numel(freq) + 0.5]);
                box on
                
            end
        end
        
        
        %% adjusting limits
        for sp1= 1:size(sph,1)
                C_lim=get(sph(sp1,:),'clim');
                if isa(C_lim,'cell')
                c_lim=max(abs([C_lim{:}]));
                else
                c_lim=max(abs(C_lim));
                end
            for sp2= 1:size(sph,2)
                % sig %
                ax=sph(sp1,sp2);
                %set(ax,'ylim',y_lim);
                subplot(ax);
                %line([0 0], y_lim, 'color', 'k');
                subplot(ax);
                cm=customcolormap_preset('red-yellow-blue',256);
                colormap(ax,cm);
                cb = colorbar;
                set(cb,'position',get(cb,'position')+[0.05 0 0 0]);
                set(get(cb,'title'),'string', ['n. ' par], 'fontsize',8);
                set(ax,'clim',[-c_lim c_lim]);
            end
        end
        results_file=['Y:\Projects\Pulv_bodysignal\Figures\' par '_' comp];
        wanted_size=[50 30];
        set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])    %
        export_fig(h, results_file, '-pdf');
        
    end
end
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
