function ecg_bna_plots_per_mua_site( data, con, cfg, PlotMethod,varargin )


% number of subplots required
plot_names={'MUA'};
results_folder=[cfg.fldr.([cfg.to_trigger '_sites']) filesep];

%% Smoothing Kernel here:
%win = 1:cfg.mua.smoothWin; win=win-(numel(win)+1)/2;
win = 1:cfg.lfp.smoothWin; win=win-(numel(win)+1)/2;


half_win = ceil(size(win,2)/2)-1;
gaussian_kernel=normpdf(win,0,numel(win)/6);
gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);

%% loop through conditions
for cn= 1:numel(data.(con))
    if isempty(fieldnames(data.(con)(cn))) || isempty(data.(con)(cn).event) % &&~isempty([sites_data(i).(con)(cn).event.mua])
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
    concat.mua = [];
    concat.mua_sgnf = [];
    concat.mua_time = [];    
    concat.mua_shufmean = [];
    concat.mua_shufstd = [];
    concat.mua_shufconf = [];
    
    ticksamples_mua=[];
    for e = 1:size(cfg.analyse_states,1)% size(con_data, 2)
        shufmean=con_data(e).surrogate.evoked.mean;
        shufstd=con_data(e).surrogate.evoked.std;
        shufconf=con_data(e).surrogate.evoked.conf95;
        
        mua_sgnf        =con_data(e).significance.evoked;
        mua_mean        =con_data(e).(PlotMethod).evoked.mean;
        mua_time        =con_data(e).time;
        
        onset_s = find(mua_time <= 0, 1, 'last'); % state onset time
        start_s = 1; % state start time
        finish_s = length(mua_time); % start end sample        
        ticksamples_mua=[ticksamples_mua length(concat.mua_time)+[start_s onset_s finish_s]];
        
        %% smooth here !
        mua_mean    =smooth_bp(mua_mean,    gaussian_kernel,half_win);
        
        % concatenate across events with a NaN separation in between
        NaNseparator=100/25;
        concat.mua          = cat(3, concat.mua,        mua_mean,   nan(size(mua_mean, 1),      size(mua_mean, 2),   NaNseparator));
        concat.mua_sgnf     = cat(3, concat.mua_sgnf,   mua_sgnf,   nan(size(mua_sgnf, 1),      size(mua_sgnf, 2),   NaNseparator));
                
        concat.mua_shufmean = cat(3, concat.mua_shufmean,  shufmean,  nan(size(shufmean, 1),     size(shufmean, 2),  NaNseparator));
        concat.mua_shufstd  = cat(3, concat.mua_shufstd,   shufstd,   nan(size(shufstd, 1),      size(shufstd, 2),   NaNseparator));
        concat.mua_shufconf  = cat(3, concat.mua_shufconf,   shufconf,   nan(size(shufconf, 1),      size(shufconf, 2),   NaNseparator));
        
        concat.mua_time     = [concat.mua_time, mua_time,     nan(1, NaNseparator)];
        
    end
    
    %% plot    
    mua_events.onset        = find(concat.mua_time == 0);
    %mua_events.name         ={con_data(states_valid).event_name};    
    mua_events.name         ={con_data.event_name};   
    mua_events.ticksamples  = sort(ticksamples_mua);
    mua_events.ticks        =round(concat.mua_time(mua_events.ticksamples)*100)/100;
        
    
    %% Evoked mua
    hold on;
    xlabel('Time(s)'); ylabel('Voltage (V)');
    title(plot_names{1},'Interpreter', 'none', 'fontsize',8);
    
    if strcmp(PlotMethod,'observed')
        lineprops={};
        % plot 95% confidence?
        conmean=concat.mua_shufmean;
        conf=concat.mua_shufconf;
        conf(1,:,:)=-conmean+conf(1,:,:);
        conf(2,:,:)=conmean-conf(2,:,:);
        shadedErrorBar(1:size(concat.mua_shufmean,3), squeeze(concat.mua_shufmean),squeeze(conf),lineprops,1);
    end
    plot(squeeze(concat.mua)','linewidth',1.5)
    line([0 0], ylim, 'color', 'k');
    ylm = get(gca,'Ylim');
    significance = double(squeeze(concat.mua_sgnf));
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
    
    add_ticks_and_labels(mua_events,ylm,diff(ylim)/10)
    set(gca, 'xlim', [0 mua_events.ticksamples(end)] + 0.5); %%should be from mua_events 
    axis square
    %% plot title...
    R=[con_data(:).observed];R=[R(:).ntriggers];
    S=[con_data(:).surrogate];S=[S(:).ntriggers];
    
    observedtriggers=[num2str(R') repmat('/',size(R'))]';observedtriggers=observedtriggers(:)';observedtriggers=strrep(observedtriggers,' ','');
    surrogatetriggers=[num2str(S') repmat('/',size(S'))]';surrogatetriggers=surrogatetriggers(:)';surrogatetriggers=strrep(surrogatetriggers,' ','');
    %plottitle = [data.site_ID ' - ' data.target ', ' data.(con)(cn).label ', ' num2str(cfg.mua.n_permutations) ' shuffles, ntriggers:' observedtriggers ' observed, ' surrogatetriggers ' surrogate'];
    plottitle = [data.site_ID ' - ' data.target ', ' data.(con)(cn).label ', ' num2str(cfg.lfp.n_permutations) ' shuffles, ntriggers:' observedtriggers ' observed, ' surrogatetriggers ' surrogate'];
    
    results_file = fullfile(results_folder, [data.site_ID '_' data.(con)(cn).label ' ' PlotMethod]);
    if strcmp(PlotMethod,'normalized')
        mtit([plottitle ' ' PlotMethod ' (' cfg.lfp.normalization ')'],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
        %mtit([plottitle ' ' PlotMethod ' (' cfg.mua.normalization ')'],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
    else
        mtit([plottitle ' ' PlotMethod],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
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


