function ecg_bna_plots_per_site( data, con,cfg, PlotMethod,varargin )


results_folder=[cfg.fldr.([cfg.to_trigger '_sites']) filesep];

% colorbar title
if strcmp(cfg.lfp.normalization, 'zscore')
    cbtitle = '(P - \mu) / std'; 
elseif strcmp(cfg.lfp.normalization, 'division')
    cbtitle = 'P / \mu';
elseif strcmp(cfg.lfp.normalization, 'subtraction')
    cbtitle = 'P - \mu';
elseif strcmp(cfg.lfp.normalization, 'relchange')
    cbtitle = '(P - \mu) / \mu';
elseif strcmp(cfg.lfp.normalization, 'none')
    cbtitle = 'not-Normalized';
end

% number of subplots required
plot_names={'POW','ITPC','Power_BP','ITPC_BP','LFP_Evoked','average phase'};
nsubplots=numel(plot_names);
ncolumns=2;
nrows=ceil(nsubplots/ncolumns);

%% Smoothing Kernel here:
win = 1:cfg.lfp.smoothWin; win=win-(numel(win)+1)/2;
half_win = ceil(size(win,2)/2)-1;
gaussian_kernel=normpdf(win,0,numel(win)/6);
gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);

%% loop through conditions
for cn= 1:numel(data.(con))
    if isempty(fieldnames(data.(con)(cn))) || isempty(data.(con)(cn).event) % &&~isempty([sites_data(i).(con)(cn).event.evoked])
        continue
    end
    % this here is eventually indicating number of triggers for each alignment
    con_data=data.(con)(cn).event;
    
    % create figure
    h = figure('units','normalized','position',[0 0 1 1]);
    
    %    states_valid=[];
    collim{1}=[];
    collim{2}=[];
    collim{6}=[];
    
    % concatenate tfs for different state windows for plotting
    concat.pow = [];
    concat.pha = [];
    concat.itpc = [];
    concat.itpcbp = [];
    concat.powbp = [];
    concat.evoked = [];
    concat.pow_sgnf = [];
    concat.itpc_sgnf = [];
    concat.itpcbp_sgnf = [];
    concat.powbp_sgnf = [];
    concat.evoked_sgnf = [];
    concat.tfr_time = [];
    concat.evoked_time = [];
    concat.evoked_shufmean = [];
    concat.evoked_shufstd = [];
    concat.evoked_shufconf = [];
    concat.freq  = cfg.lfp.foi; %% this is actually in the settings...
    
    %event_info = struct();
    ticksamples_tfr=[];
    ticksamples_evoked=[];
    for e = 1:size(cfg.analyse_states,1)%size(con_data, 2)
        shufmean=con_data(e).surrogate.evoked.mean;
        shufstd=con_data(e).surrogate.evoked.std;
        shufconf=con_data(e).surrogate.evoked.conf95;
        
        pow_sgnf        =con_data(e).significance.pow;
        itpc_sgnf       =con_data(e).significance.itpc;
        itpcbp_sgnf     =con_data(e).significance.itpcbp;
        powbp_sgnf      =con_data(e).significance.powbp;
        evoked_sgnf     =con_data(e).significance.evoked;
        
        pow_mean        =con_data(e).(PlotMethod).pow.mean;
        pha_mean        =con_data(e).(PlotMethod).pha.mean;
        itpc_mean       =con_data(e).(PlotMethod).itpc.mean;
        itpcbp_mean     =con_data(e).(PlotMethod).itpcbp.mean;
        powbp_mean      =con_data(e).(PlotMethod).powbp.mean;
        evoked_mean     =con_data(e).(PlotMethod).evoked.mean;
        
        tfr_time=con_data(e).tfr_time;
        evoked_time=con_data(e).time;
        
        onset_s = find(tfr_time <= 0, 1, 'last'); % state onset time
        start_s = 1; % state start time
        finish_s = length(tfr_time); % start end sample
        ticksamples_tfr=[ticksamples_tfr length(concat.tfr_time)+[start_s onset_s finish_s]];
        
        onset_s = find(evoked_time <= 0, 1, 'last'); % state onset time
        start_s = 1; % state start time
        finish_s = length(evoked_time); % start end sample
        ticksamples_evoked=[ticksamples_evoked length(concat.evoked_time)+[start_s onset_s finish_s]];
        
        %% smooth here !
        itpcbp_mean =smooth_bp(itpcbp_mean, gaussian_kernel,half_win);
        powbp_mean  =smooth_bp(powbp_mean,  gaussian_kernel,half_win);
        evoked_mean    =smooth_bp(evoked_mean,    gaussian_kernel,half_win);
        
        % concatenate across events with a NaN separation in between
        NaNseparator=100/25;
        concat.pow          = cat(3, concat.pow,        pow_mean,   nan(size(pow_mean, 1),      size(pow_mean, 2),   NaNseparator));
        concat.pha          = cat(3, concat.pha,        pha_mean,   nan(size(pha_mean, 1),      size(pha_mean, 2),   NaNseparator));
        concat.itpc         = cat(3, concat.itpc,       itpc_mean,  nan(size(itpc_mean, 1),     size(itpc_mean, 2),  NaNseparator));
        concat.itpcbp       = cat(3, concat.itpcbp,     itpcbp_mean,nan(size(itpcbp_mean, 1),   size(itpcbp_mean, 2),NaNseparator));
        concat.powbp        = cat(3, concat.powbp,      powbp_mean, nan(size(powbp_mean, 1),    size(powbp_mean, 2), NaNseparator));
        concat.evoked       = cat(3, concat.evoked,     evoked_mean,nan(size(evoked_mean, 1),   size(evoked_mean, 2),NaNseparator));
        concat.pow_sgnf     = cat(3, concat.pow_sgnf,   pow_sgnf,   nan(size(pow_sgnf, 1),      size(pow_sgnf, 2),   NaNseparator));
        concat.itpc_sgnf    = cat(3, concat.itpc_sgnf,  itpc_sgnf,  nan(size(itpc_sgnf, 1),   	size(itpc_sgnf, 2),  NaNseparator));
        concat.itpcbp_sgnf  = cat(3, concat.itpcbp_sgnf,itpcbp_sgnf,nan(size(itpcbp_sgnf, 1),   size(itpcbp_sgnf, 2),NaNseparator));
        concat.powbp_sgnf   = cat(3, concat.powbp_sgnf, powbp_sgnf, nan(size(powbp_sgnf, 1),    size(powbp_sgnf, 2), NaNseparator));
        concat.evoked_sgnf  = cat(3, concat.evoked_sgnf,evoked_sgnf,nan(size(evoked_sgnf, 1),   size(evoked_sgnf, 2),NaNseparator));
        
        
        concat.evoked_shufconf = cat(3, concat.evoked_shufconf,  shufconf,  nan(size(shufconf, 1),     size(shufconf, 2),  NaNseparator));
        concat.evoked_shufmean = cat(3, concat.evoked_shufmean,  shufmean,  nan(size(shufmean, 1),     size(shufmean, 2),  NaNseparator));
        concat.evoked_shufstd  = cat(3, concat.evoked_shufstd,   shufstd,   nan(size(shufstd, 1),      size(shufstd, 2),   NaNseparator));
        
        concat.tfr_time     = [concat.tfr_time, tfr_time, nan(1, NaNseparator)];
        concat.evoked_time     = [concat.evoked_time, evoked_time,     nan(1, NaNseparator)];
        
        %         % somehow needed for (not) labelling not existing alignments
        %         if ~all(isnan(tfr_time))
        %             states_valid=[states_valid e];
        %         end
    end
    
    %% plot
    tfr_events.onset        = find(concat.tfr_time == 0);
    %tfr_events.name         ={con_data(states_valid).event_name};
    tfr_events.name         ={con_data.event_name};
    tfr_events.ticksamples  = sort(ticksamples_tfr);
    tfr_events.startsamples = ticksamples_tfr(1:3:end);
    tfr_events.endsamples   = ticksamples_tfr(3:3:end);
    tfr_events.ticks        =round(concat.tfr_time(tfr_events.ticksamples)*100)/100;
    
    evoked_events.onset        = find(concat.evoked_time == 0);
    %evoked_events.name         ={con_data(states_valid).event_name};
    evoked_events.name         ={con_data.event_name};
    evoked_events.ticksamples  = sort(ticksamples_evoked);
    evoked_events.ticks        =round(concat.evoked_time(evoked_events.ticksamples)*10)/10;
    
    
    %% POW and ITPC
    toplot={real(concat.pow),real(concat.itpc),real(concat.powbp),real(concat.itpcbp),0,real(concat.pha)};
    sigplot = {concat.pow_sgnf, concat.itpc_sgnf,concat.powbp_sgnf,concat.itpcbp_sgnf};
    for sp=1:2 % frequency spectra
        sph(sp)=subplot(nrows, ncolumns, sp);
        image(1:size(toplot{sp},3), 1:numel(concat.freq), squeeze(toplot{sp}),'CDataMapping','scaled'); %% complex for task differences??
        set(gca,'YDir','normal');
        hold on;
        % calculate the significance here:
        significance = double(squeeze(sigplot{sp})); %+repmat(randi([0 1], 60,19),[1,11]); checking the plot results
        %significance(significance==0)=NaN;
        % plotting a contour around the significant parts:
        %[cout, hand] = contour(1:size(toplot{sp},3),1:numel(concat.freq),significance,1,'linecolor','k')
        ecg_bna_draw_outlines(significance==1,'r')
        ecg_bna_draw_outlines(significance==-1,'b')
        
        nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
        collim{sp}=[min([collim{sp}(:); nonnan(:)]) max([collim{sp}(:); nonnan(:)])];
%         if  strcmp(PlotMethod,'surrogate')
%             collim{1}=[min([collim{sp}(:); nonnan(:)]) max([collim{sp}(:); nonnan(:)])];
%             collim{2}=[0 0.55];
%         else
%             collim{sp}=[min([collim{sp}(:); nonnan(:)]) max([collim{sp}(:); nonnan(:)])];
%         end
%         
        % horizontal lines to separate frequency bands
        fbandstart = unique(cfg.lfp.frequency_bands(:))';
        fbandstart(fbandstart<min(cfg.lfp.foi))=[];
        fbandstart_idx = zeros(size(fbandstart));
        for f = fbandstart
            f_idx = find(abs(concat.freq - f) == min(abs(concat.freq - f)), 1, 'first');
            line([tfr_events.startsamples' tfr_events.endsamples']-1/2, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
            fbandstart_idx(fbandstart == f) = f_idx;
        end
        
        set(gca,'TickDir','out')
        set(gca, 'ytick', fbandstart_idx);
        set(gca, 'yticklabel', fbandstart);
        % add 0.5 at end since the time value is the center of the bin
        % add 0 at beginning to make x-axis visible
        set(gca, 'ylim', [0.5,numel(concat.freq) + 0.5]);
        add_ticks_and_labels(tfr_events,[0.5,numel(concat.freq) + 0.5],8)
        
        set(gca, 'xlim', [0 tfr_events.ticksamples(end)] + 0.5);
        ylabel('Frequency (Hz)');
        title(plot_names{sp},'Interpreter', 'none', 'fontsize',8);
    end
    
    %% Bandpassed POW and ITPC
    for sp=3:4
        sph(sp)=subplot(nrows, ncolumns, sp);             %% change color order to something nicer
        hold on;
        set(gca,'ColorOrder',cool(size(toplot{sp},2)));
        %plot(repmat(concat.evoked_time,size(concat.powbp,2),1)', squeeze(concat.powbp)')
        plot(squeeze(toplot{sp})')
        xlabel('Time(s)'); ylabel('Power (W)');
        
        % adding the signifiance horizontal lines:
        ylm = get(gca,'Ylim');
        stp = (ylm(2)-ylm(1))/20;
        ylm(1)=ylm(1)-size(toplot{sp},2)*stp;
        set(gca,'Ylim',ylm);
        
        significance = double(squeeze(sigplot{sp}));         % i needed to create concat.itpcbp_sgnf, it basically appends Nans for a (potential) separator with a second alignment
        significance(significance==0)=NaN;                          % replacing zeros with Nans means once we plot, lines will be discontinoous there
        multiplicator= ylm(1)+(1:size(significance,1))*stp;              % multiplicator basically defines position of significance line
        significance=significance.*repmat(multiplicator',1,size(significance,2));
        X=~all(isnan(diff(significance,1,2)),2);
        plot(significance(X,:)','linewidth',3);
        
        add_ticks_and_labels(evoked_events,ylm,stp)
        
        legend({strcat(num2str(round(cfg.lfp.frequency_bands(:,1))), '-',num2str(round(cfg.lfp.frequency_bands(:,2))), ' Hz')},'fontsize',3);
        title(plot_names{sp},'Interpreter', 'none', 'fontsize',8);
        set(gca, 'xlim', [0 tfr_events.ticksamples(end)] + 0.5); %%should be from evoked_events
    end
    
    %% Phase
    sp=6; % frequency spectra
    sph(sp)=subplot(nrows, ncolumns, sp);
    image(1:size(toplot{sp},3), 1:numel(concat.freq), squeeze(toplot{sp}),'CDataMapping','scaled');
    set(gca,'YDir','normal');
    hold on;
    
    nonnan=toplot{sp};nonnan(isnan(nonnan))=[];
    collim{sp}=[min([collim{sp}(:); nonnan(:)]) max([collim{sp}(:); nonnan(:)])];
    
    % horizontal lines to separate frequency bands
    fbandstart = unique(cfg.lfp.frequency_bands(:))';
    fbandstart(fbandstart<min(cfg.lfp.foi))=[];
    fbandstart_idx = zeros(size(fbandstart));
    for f = fbandstart
        f_idx = find(abs(concat.freq - f) == min(abs(concat.freq - f)), 1, 'first');
        line([tfr_events.startsamples' tfr_events.endsamples']-1/2, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
        fbandstart_idx(fbandstart == f) = f_idx;
    end
    
    set(gca,'TickDir','out')
    set(gca, 'ytick', fbandstart_idx);
    set(gca, 'yticklabel', fbandstart);
    % add 0.5 at end since the time value is the center of the bin
    % add 0 at beginning to make x-axis visible
    set(gca, 'ylim', [0.5,numel(concat.freq) + 0.5]);
    add_ticks_and_labels(tfr_events,[0.5,numel(concat.freq) + 0.5],8)
    
    set(gca, 'xlim', [0 tfr_events.ticksamples(end)] + 0.5);
    ylabel('Frequency (Hz)');
    title(plot_names{sp},'Interpreter', 'none', 'fontsize',8);
    
    
    %% Evoked LFP
    sp=5;
    sph(sp)=subplot(nrows, ncolumns, sp);
    hold on;
    plot(squeeze(concat.evoked)','linewidth',1.5)
    line([0 0], ylim, 'color', 'k');
    xlabel('Time(s)'); ylabel('Voltage (V)');
    title(plot_names{sp},'Interpreter', 'none', 'fontsize',8);
    
    if strcmp(PlotMethod,'observed')
        lineprops={};
        
        % plot 95% confidence?
        conmean=concat.evoked_shufmean;
        conf=concat.evoked_shufconf;
        conf(1,:,:)=-conmean+conf(1,:,:);
        conf(2,:,:)=conmean-conf(2,:,:);
        shadedErrorBar(1:size(concat.evoked_shufmean,3), squeeze(concat.evoked_shufmean),squeeze(conf),lineprops,1);
        %shadedErrorBar(1:size(concat.evoked_shufmean,3), squeeze(concat.evoked_shufmean),squeeze(concat.evoked_shufstd),lineprops,1);
    end
    ylm = get(gca,'Ylim');
    significance = double(squeeze(concat.evoked_sgnf));
    signeg=significance;
    sigpos=significance;
    signeg(signeg~=-1)=NaN;
    sigpos(sigpos~=1)=NaN;
    signeg=abs(signeg).*ylm(1);
    sigpos=abs(sigpos).*ylm(1)+diff(ylm)/40;    
%     if any(~isnan(signeg) & [false; diff(signeg)==0])
%         plot(1:numel(signeg),signeg','linewidth',3);
%     end
%     if any(~isnan(sigpos) & [false; diff(sigpos)==0])
%         plot(1:numel(sigpos),sigpos','linewidth',3);
%     end
    
        plot(1:numel(signeg),signeg','b','linewidth',3);
        plot(1:numel(sigpos),sigpos','r','linewidth',3);
    add_ticks_and_labels(evoked_events,ylm,diff(ylim)/10)
    set(gca, 'xlim', [0 tfr_events.ticksamples(end)] + 0.5); %%should be from evoked_events
    % end
    %% format spectra colors
    collim{3}=collim{1};
    collim{4}=collim{1};
    collim{5}=collim{1};
    for sp=1:6
        subplot(sph(sp));
        set(gca,'CLim',real(collim{sp}))
        cm = colormap('jet');
        if nargin > 4
            cm = colormap(varargin{1});
        end
        cb = colorbar('EastOutside');%('North');
        if strcmp(PlotMethod,'normalized')
            set(get(cb,'title'),'string', cbtitle, 'fontsize',8);
        else
            if sp==1
                set(get(cb,'title'),'string', 'power(W)', 'fontsize',8);
            elseif sp==2
                set(get(cb,'title'),'string', 'ITPC', 'fontsize',8);
            end
        end
        switch sp
            case {1,2,6}
                colormap(cm);
            case {3,4,5}
                set(cb,'Visible','off');
        end
        
    end
    
    
    %% plot title...
    R=[con_data(:).observed];R=[R(:).ntriggers];
    S=[con_data(:).surrogate];S=[S(:).ntriggers];
    
    observedtriggers=[num2str(R') repmat('/',size(R'))]';observedtriggers=observedtriggers(:)';observedtriggers=strrep(observedtriggers,' ','');
    surrogatetriggers=[num2str(S') repmat('/',size(S'))]';surrogatetriggers=surrogatetriggers(:)';surrogatetriggers=strrep(surrogatetriggers,' ','');
    plottitle = [data.site_ID ' - ' data.target ' ' data.(con)(cn).label ' ' PlotMethod ', ' num2str(cfg.lfp.n_permutations) ' surrogates, ntriggers observed: ' observedtriggers ',surrogate: ' surrogatetriggers];
    %ntriggers = [ observedtriggers ' observed, ' surrogatetriggers ' surrogate'];
    
    results_file = fullfile(results_folder, [data.site_ID '_' data.(con)(cn).label ' ' PlotMethod]);
    if strcmp(PlotMethod,'normalized')
        mtit([plottitle   ' (' cfg.lfp.normalization ')'],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
    else
        mtit(plottitle,'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
    end
    
    
    wanted_size=[50 30];
    set(h, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
    
    export_fig(h, results_file, '-pdf'); %% how come this does not export most plots ??
end
end


function add_ticks_and_labels(events,ylm,stp)


%state_ticks=round(concat.tfr_time(event_samples)*10)/10;

event_onsets    =events.onset;
event_names     =events.name;
event_samples   =events.ticksamples;
state_ticks     =events.ticks;

for so = event_onsets
    line([so so], ylm, 'color', 'k');
    %if isfield(con_data(event_onsets == so), 'event_name') && ~isempty(event_names(event_onsets == so))
    event_name = event_names{event_onsets == so};
    event_name=strrep(event_name,'_',' ');
    text(so+1, ylm(1)+stp, event_name, 'fontsize', 8, 'fontweight', 'bold');
    %end
end

% mark state onsets
set(gca,'xtick',event_samples(~isnan(state_ticks)))
set(gca,'xticklabels', state_ticks(~isnan(state_ticks)), 'fontsize', 6)
% add 0.5 since the time value is the center of the bin
% add 0 at the beginning to make the y-axis visib
end

function smoothed=smooth_bp(input,gaussian_kernel,half_win)

jnk = [];
concat_input = cat(3,(input(:,:,half_win:-1:1)),(input(:,:,:)));
concat_input = cat(3,concat_input, (input(:,:,end:-1:end-half_win+1)));
for m = 1: size(input,1)
    for k=1:size(input,2)
        jnk(m,k,:)= conv(squeeze(concat_input(m,k,:)),gaussian_kernel,'same');%./max(conv(ones(100,1), gausswin(win)));
    end
end
clear concat_powbp
smoothed = jnk(:,:,half_win+1:end-half_win);
end


