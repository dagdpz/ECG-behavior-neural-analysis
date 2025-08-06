function ecg_bna_plots_per_ecg_site( data, cfg, PlotMethod,varargin )


% number of subplots required
plot_names={'ECG'};
results_folder=[cfg.sites_ecg_fldr filesep];

%% Smoothing Kernel here:
%win = 1:cfg.ecg.smoothWin; win=win-(numel(win)+1)/2;
win = 1:cfg.lfp.smoothWin; win=win-(numel(win)+1)/2;


half_win = ceil(size(win,2)/2)-1;
gaussian_kernel=normpdf(win,0,numel(win)/6);
gaussian_kernel=gaussian_kernel/sum(gaussian_kernel);

%% loop through conditions
for cn= 1:numel(data.condition)
    if isempty(fieldnames(data.condition(cn))) || isempty(data.condition(cn).event) % &&~isempty([sites_data(i).condition(cn).event.ecg])
        continue
    end    
    % this here is eventually indicating number of triggers for each alignment
    con_data=data.condition(cn).event;
    
    % create figure
    h = figure('units','normalized','position',[0 0 1 1]);
    
%    states_valid=[];
    collim{1}=[];
    collim{2}=[];
    collim{6}=[];
    
    % concatenate tfs for different state windows for plotting
    concat.ecg = [];
    concat.ecg_sgnf = [];
    concat.ecg_time = [];    
    concat.ecg_shufmean = [];
    concat.ecg_shufstd = [];
    
    ticksamples_ecg=[];
    for e = 1:size(con_data, 2)
        shufmean=con_data(e).shuffled.ecg.mean;
        shufstd=con_data(e).shuffled.ecg.std;
        
        ecg_sgnf        =con_data(e).significance.ecg;
        ecg_mean        =con_data(e).(PlotMethod).ecg.mean;
        ecg_time        =con_data(e).time;
        
        onset_s = find(ecg_time <= 0, 1, 'last'); % state onset time
        start_s = 1; % state start time
        finish_s = length(ecg_time); % start end sample        
        ticksamples_ecg=[ticksamples_ecg length(concat.ecg_time)+[start_s onset_s finish_s]];
        
        %% smooth here !
        ecg_mean    =smooth_bp(ecg_mean,    gaussian_kernel,half_win);
        
        % concatenate across events with a NaN separation in between
        NaNseparator=100/25;
        concat.ecg          = cat(3, concat.ecg,        ecg_mean,   nan(size(ecg_mean, 1),      size(ecg_mean, 2),   NaNseparator));
        concat.ecg_sgnf     = cat(3, concat.ecg_sgnf,   ecg_sgnf,   nan(size(ecg_sgnf, 1),      size(ecg_sgnf, 2),   NaNseparator));
                
        concat.ecg_shufmean = cat(3, concat.ecg_shufmean,  shufmean,  nan(size(shufmean, 1),     size(shufmean, 2),  NaNseparator));
        concat.ecg_shufstd  = cat(3, concat.ecg_shufstd,   shufstd,   nan(size(shufstd, 1),      size(shufstd, 2),   NaNseparator));
        
        concat.ecg_time     = [concat.ecg_time, ecg_time,     nan(1, NaNseparator)];
        
    end
    
    %% plot    
    ecg_events.onset        = find(concat.ecg_time == 0);
    %ecg_events.name         ={con_data(states_valid).event_name};    
    ecg_events.name         ={con_data.event_name};   
    ecg_events.ticksamples  = sort(ticksamples_ecg);
    ecg_events.ticks        =round(concat.ecg_time(ecg_events.ticksamples)*10)/10;
        
    
    %% Evoked ecg
    hold on;
    xlabel('Time(s)'); ylabel('Voltage (V)');
    title(plot_names{1},'Interpreter', 'none', 'fontsize',8);
    
    if strcmp(PlotMethod,'real')
        lineprops={};
        shadedErrorBar(1:size(concat.ecg_shufmean,3), squeeze(concat.ecg_shufmean),squeeze(concat.ecg_shufstd),lineprops,1);
    end
    plot(squeeze(concat.ecg)','linewidth',1.5)
    line([0 0], ylim, 'color', 'k');
    ylm = get(gca,'Ylim');
    significance = double(squeeze(concat.ecg_sgnf));
    significance(significance==0)=NaN;
    significance=significance.*ylm(1);
    if any(~isnan(significance) & [false; diff(significance)==0])
        plot(1:numel(significance),significance','linewidth',3);
    end
    add_ticks_and_labels(ecg_events,ylm,diff(ylim)/10)
    set(gca, 'xlim', [0 ecg_events.ticksamples(end)] + 0.5); %%should be from ecg_events 
    
    %% plot title...
    R=[con_data(:).real];R=[R(:).ntriggers];
    S=[con_data(:).shuffled];S=[S(:).ntriggers];
    
    realtriggers=[num2str(R') repmat('/',size(R'))]';realtriggers=realtriggers(:)';realtriggers=strrep(realtriggers,' ','');
    shuffledtriggers=[num2str(S') repmat('/',size(S'))]';shuffledtriggers=shuffledtriggers(:)';shuffledtriggers=strrep(shuffledtriggers,' ','');
    %plottitle = [data.site_ID ' - ' data.target ', ' data.condition(cn).label ', ' num2str(cfg.ecg.n_permutations) ' shuffles, ntriggers:' realtriggers ' real, ' shuffledtriggers ' shuffled'];
    plottitle = [data.site_ID ' - ' data.target ', ' data.condition(cn).label ', ' num2str(cfg.lfp.n_permutations) ' shuffles, ntriggers:' realtriggers ' real, ' shuffledtriggers ' shuffled'];
    
    results_file = fullfile(results_folder, [data.site_ID '_' data.condition(cn).label ' ' PlotMethod]);
    if strcmp(PlotMethod,'normalized')
        mtit([plottitle ' ' PlotMethod ' (' cfg.lfp.normalization ')'],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
        %mtit([plottitle ' ' PlotMethod ' (' cfg.ecg.normalization ')'],'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none')
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


