function allSitesData = ecg_bna_remove_rawLFP_outliers(sitesdir,sitefiles,cfg,ts_original,sr,Triggers,blocks,blockstart)
cfg.ica_sampling_rate = 256;  
ts=round(sr/cfg.ica_sampling_rate); %% cfg.lfp.timestep is binning used in process_lfp, at the moment it's 10 ms

%% but lfp 
session = strsplit(sitefiles(1).name,'_');
session = session{3};

reprocess = 0; % to generate the "allSitesData" and removing the noisy channels,
% and calculating the ICA weights. Otherweise directly load the
% "sampling_table" and "allSitesData" for each session and go further.

if reprocess==1
    allSitesData = [];
    
    temp_allSites = struct;
    nSites = length(sitefiles);
    for s = 1: nSites
        load([sitesdir,filesep,sitefiles(s).name],'sites');
        
        temp_allSites(s).site = sites;
        temp_allSites(s).site.LFP = sites.LFP;
        temp_allSites(s).name = sitefiles(s).name;
        temp_allSites(s).inconsistant_ch = 0;
        temp_allSites(s).noisy_site = 0;
    end
    
    targets = {'VPL_R', 'VPL_L', 'dPul_R', 'dPul_L','MD_L','MD_R'};
    allSites_targets = arrayfun(@(x) (x.site.target), temp_allSites, 'UniformOutput', false);
    inval_targets = ~ismember(allSites_targets,targets);
    temp_allSites(inval_targets) = [];
    
    trials_block = arrayfun(@(x) unique(x.site.block), temp_allSites, 'UniformOutput', false);
    blocks_with_LFP = unique([trials_block{:}]);
    
    dummy = 1;
    for s = 1: length(temp_allSites)
        samples_past=0;
        % you are not resampling here :P --> Now I am :D
        samples_past_resampled=0;
        
        for b = 1:length(blocks_with_LFP)
            B = blocks_with_LFP(b);
            if sum(ismember(unique(temp_allSites(s).site.block),B))>0
                sampling_table.site(dummy)  = s;
                sampling_table.block(dummy) = B;
                
                b_samples=[temp_allSites(s).site.block ==B];
                bs=samples_past_resampled+1;
                be=floor(sum(temp_allSites(s).site.LFP_samples(b_samples))/ts)+samples_past_resampled;
                bs_original=samples_past+1;
                be_original=sum(temp_allSites(s).site.LFP_samples(b_samples))+samples_past;
                
                sampling_table.start_sample(dummy)  = bs_original;
                sampling_table.end_sample(dummy)    = be_original;
                
                samples_past=be_original;
                %samples_past_resampled=be;
                %%
                % Here I generate the number of sample per blocks to feed for Triggers:
                sampling_table.tfs(s).n_ica_samples_per_block(:,b)=[B,be-bs+1];
                sampling_table.tfs(s).n_samples_per_block(:,b)=[B,be_original-bs_original+1];
                
                dummy = dummy+1;
            end
        end
        % Here I generate the triggers to check the ICA epochs later on:
        n_LFP_samples_per_block = sampling_table.tfs(s).n_samples_per_block; 
        sampling_table.site_triggers(s)= ecg_bna_resample_triggers2(Triggers,[blocks;blockstart(blocks)],n_LFP_samples_per_block,sr);
        n_LFP_samples_per_block = sampling_table.tfs(s).n_ica_samples_per_block; 
        sampling_table.site_triggers_ica(s)= ecg_bna_resample_triggers2(Triggers,[blocks;blockstart(blocks)],n_LFP_samples_per_block,cfg.ica_sampling_rate);
    end
    
    
    lfp_data = arrayfun(@(x) x.site.LFP, temp_allSites, 'UniformOutput', false);
    lfp_sizes = cellfun(@(x) size(x, 2), arrayfun(@(s) s.site.LFP, temp_allSites, 'UniformOutput', false));
    
    % idea was not to check across all sites
    all_LFP_concat = horzcat(lfp_data{:});
    all_sites_25perc = prctile(abs(all_LFP_concat),25);
    all_sites_75perc = prctile(abs(all_LFP_concat),75);
    
    noisy_site = [];
    for s = 1:length(temp_allSites)
        
        % but rather across all sites other than the one you are currently checking
        site_idx=1:numel(lfp_data);
        all_LFP_concat = horzcat(lfp_data{site_idx~=s});
        all_sites_25perc = prctile(abs(all_LFP_concat),25);
        all_sites_75perc = prctile(abs(all_LFP_concat),75);
        
        tmp_site_25perc = prctile(abs(lfp_data{s}),25);
        tmp_site_75perc = prctile(abs(lfp_data{s}),75);
        if (tmp_site_75perc < all_sites_25perc)||(tmp_site_25perc > all_sites_75perc)
            noisy_site = [noisy_site, s];
        end
    end
    
    if ~isempty(noisy_site)
        [temp_allSites(noisy_site).noisy_site] = deal(1);
    end
    sites_to_remove=[temp_allSites.noisy_site]==1;
    allSitesData = temp_allSites(~sites_to_remove);
    allSitesData = rmfield(allSitesData,{'inconsistant_ch','noisy_site'});
    
    
    %% sites that are removed dont need to be fixed, correct?
    %% if so, remove them from sampling table immediately
    %% alternatively, you can do noise detection BEFORE creating the sampling table, so that you do not run into this problem in the first place
    to_remove_from_sampling_table=ismember(sampling_table.site,find(sites_to_remove)); % please check this line, I'm not sure it does what it supposed to
    sampling_table.site(to_remove_from_sampling_table)  = [];
    sampling_table.block(to_remove_from_sampling_table) = [];
    sampling_table.start_sample(to_remove_from_sampling_table)  = [];
    sampling_table.end_sample(to_remove_from_sampling_table)    = [];
    
    % clear temp_allSites, because it's huge and should not be needed any more ?
    clear temp_allSites
    
    % storing the all sites sampling_table and allSitesData
    smplTblfilename = fullfile([cfg.analyse_lfp_folder,filesep,...
        cfg.session_info(1).Monkey(1:3),'_', session , '_sampleTble_allSitesdata']);
    save(smplTblfilename,'sampling_table','allSitesData','-v7.3');
    
else
    smplTblfilename = fullfile([cfg.analyse_lfp_folder,filesep,...
        cfg.session_info(1).Monkey(1:3),'_', session , '_sampleTble_allSitesdata']);
    load([smplTblfilename,'.mat']);
end
%%

lfp_tagets = arrayfun(@(x) x.site.target, allSitesData, 'UniformOutput', false);
site_names = arrayfun(@(x) x.name, allSitesData, 'UniformOutput', false); %% are these site names?
target_list = unique(lfp_tagets);
%%
lfp_data = arrayfun(@(x) x.site.LFP, allSitesData, 'UniformOutput', false);
%%
% ICA processing parameters:
reprocess = 0;
doTimecourse = 0;
icaw_path = [cfg.analyse_lfp_folder,filesep,'Per_Site_ICAweights'];
allICAweights = dir([icaw_path,filesep,'*.mat']);
%
if (isfield(cfg.lfp, 'runICA') && cfg.lfp.runICA==1) && (reprocess==1)
    sites_hemisphere=contains(lfp_tagets,'_L')+2*contains(lfp_tagets,'_R'); %look for the _ as well, because f.e. VPL contains L
    hemisphere_list = unique(sites_hemisphere);

    for h = 1: length(hemisphere_list)
        
        hemisphere_sites = sites_hemisphere==hemisphere_list(h);
        lfp_data_concat = vertcat(lfp_data{hemisphere_sites});
        hemisphere_site_names = site_names(hemisphere_sites);
        
        data = [];
        data.label = arrayfun(@(x) sprintf('LFP%02d', x), 1:size(lfp_data_concat, 1), 'UniformOutput', false);
        data.fsample = sr;  % Sampling frequency
        data.time = {(0:size(lfp_data_concat,2)-1) / sr}; % Cell array containing time vector
        data.trial = {lfp_data_concat}; % Cell array containing LFP data matrix
        data.sampleinfo = [1 size(lfp_data_concat,2)]; % Start and end sample indices
        
        ft_datatype_raw(data) % To Check if the structure is correctly formatted
        data_orig = data;
        
        cfg_ica              = [];
        cfg_ica.resamplefs   = cfg.ica_sampling_rate;
        cfg_ica.detrend      = 'no';
        data = ft_resampledata(cfg_ica, data_orig);
        
        cfg_ica              = [];
        cfg_ica.method       = 'runica';
        cfg_ica.numcomponent = round(size(lfp_data_concat,1)/2); % half the size of nSites
        cfg_ica.channel      = {'LFP'};
        data_comp = ft_componentanalysis(cfg_ica, data);
        data_to_plot = data_comp;
        
%         data_to_plot            = rmfield(data_comp,{'topo','topolabel'});
        
        % storing the sites ICA weights
        icawfile = [cfg.session_info(1).Monkey(1:3),'_', session , '_',target_list{h},'_ICAweigths_w ',num2str(length(hemisphere_site_names)),'cmp'];
        ICAfilename = fullfile([icaw_path,filesep,icawfile]);
        save(ICAfilename,'data_to_plot','hemisphere_site_names','-v7.3');
        
        fprintf("\n ============================================ \n investigating %s\n",session)
        triggers = sampling_table.site_triggers_ica(1).R_real;
        [ica_epochs, ica_avg] = ecg_bna_epochICAaroundTriggers([cfg.analyse_lfp_folder,filesep,'Per_Site_ICAweights'],session,data_to_plot, triggers, 0.25, 0.25, true,true);
        
        site_depth_idx = ismember({allSitesData.name},hemisphere_site_names');
        site_depth = arrayfun(@(x) x.site.electrode_depth, allSitesData(site_depth_idx), 'UniformOutput', false);
        site_depth = horzcat(site_depth{:});
        depths = site_depth;
        n_channels = length(data_to_plot.label);
        % ICA weight per Site plots
        fig1 = figure;
        for cmp = 1:n_channels
            
            comp_index = cmp;  
            weights = data_to_plot.unmixing(comp_index,:)';
            % Normalize for color mapping
            weights_norm = (weights - min(weights)) / (max(weights) - min(weights)); 
            subplot(6,6,cmp)
            scatter(weights, -depths, 100, weights_norm, 'filled');  % invert Y for depth
            colormap(jet);
            colorbar;
            xlabel('ICA weight');
            ylabel('Electrode Depth');
            title(sprintf('ICA Component %d Spatial Weights (LFP Channels)', comp_index),'FontSize',7);
            set(gca, 'YDir', 'normal');
            grid on;
        end
        fig1.WindowState = 'maximized';
        sgtitle(icawfile(1:end-12),'interpreter','none')
        figname2save = fullfile([icaw_path,filesep,icawfile(1:end-4)]);
        saveas(fig1,figname2save,'fig');
        saveas(fig1,figname2save,'png');
        
        % ICA Inv weight per Site plots
        weights_inv = pinv(data_to_plot.unmixing);
        fig2 = figure;
        for cmp = 1:n_channels
            
            comp_index = cmp;          
            weights_inv_2plot = weights_inv(:,comp_index)';
            weights_inv_norm = (weights_inv_2plot - min(weights_inv_2plot)) / (max(weights_inv_2plot) - min(weights_inv_2plot));
            
            subplot(6,6,cmp)
            scatter(weights_inv_2plot, -depths, 100, weights_inv_norm, 'filled');  % invert Y for depth 
            colormap(jet);
            colorbar;
            xlabel('ICA Inverse weight');
            ylabel('Electrode Depth');
%             ylabel('ICA cmp');
%             title(sprintf('Inverse ICA Component in Site: %d ', comp_index),'FontSize',7);
            title(sprintf('Inv ICA Component %d ', comp_index),'FontSize',7);
            set(gca, 'YDir', 'normal');
            grid on;
        end
        fig2.WindowState = 'maximized';
        sgtitle([icawfile(1:end-12),'_Inv'],'interpreter','none')
        figname2save = fullfile([icaw_path,filesep,icawfile(1:end-4),'_Inv2']);
        saveas(fig2,figname2save,'fig');
        saveas(fig2,figname2save,'png');
        
        
        if doTimecourse ==1
            cfg_ica                 = [];
            cfg_ica.continuous      = 'yes';
            cfg_ica.blocksize       = 30; % show 30 seconds at the time
            cfg_ica.plotevents      = 'no';
            cfg_ica.preproc.demean  = 'no';
            ica_com = ft_databrowser(cfg_ica, data_to_plot);
        else
            ica_com = data_to_plot;
        end
        
        BadICAcomps = inputdlg('Enter space-separated Selected ICA Components to remove:','Sample', [1 64]);
        BadICAcomps = str2num(BadICAcomps{1});
        close all
        
        cfg_ica = [];
        cfg_ica.component = BadICAcomps; % to be removed
        data_orig_clean = ft_rejectcomponent(cfg_ica, ica_com, data_orig);
        
        for sitei = 1:length(find(hemisphere_sites))
            site_idx = find(hemisphere_sites);
            allSitesData(site_idx(sitei)).site.LFP = data_orig_clean.trial{1, 1}(site_idx(sitei),:);
        end
        
        % % % Another version of ICA, directly from fastICA function:
        % % addpath('Y:\Projects\_Shamim\fieldtrip-20231130\external\fastica'); % feed the location of fastica from somewhere local!
        % % % --->> keep the same fieldtrip version for consistency
        % % % [icasig, A, W] = fastica (1000*lfp_data_concat); % outputs the estimated separating matrix W and the corresponding mixing matrix A.
        % % [icasig2, A2, W2] = fastica (1000*lfp_data_concat,  'numOfIC',20); % outputs the estimated separating matrix W and the corresponding mixing matrix A.
        % % % icaplot('complot',icasig,[1:20]);
        % % % sr = 1.0173e+03;
        % % icaplot('complot',icasig2);
        % % set(gca,'xlim',[0 10*sr])
        % ecg_bna_plot_PSD_LFP(icasig, sr)
        
    end
    
elseif (isfield(cfg.lfp, 'runICA') && cfg.lfp.runICA==1) && (reprocess==0)
    
    sites_hemisphere=contains(lfp_tagets,'_L')+2*contains(lfp_tagets,'_R');
    hemisphere_list = unique(sites_hemisphere);
    
    for h = 1: length(hemisphere_list)
        
        icaw_fileidx = contains({allICAweights.name},[session,'_',target_list{h}]);
        hemisphere_sites = sites_hemisphere==hemisphere_list(h);
        hemisphere_site_names = site_names(hemisphere_sites);
%         icawfile = [cfg.session_info(1).Monkey(1:3),'_', session , '_',target_list{h},'_ICAweigths_w ',num2str(length(hemisphere_site_names)),'cmp'];
        icawfile = allICAweights(icaw_fileidx).name;
        ICAfilename = fullfile([icaw_path,filesep,icawfile]);
        load(ICAfilename);
        data_to_plot.topo = data_to_plot.unmixing';
        data_to_plot.topolabel = data_to_plot.cfg.channel';
        lfp_data_concat = vertcat(lfp_data{hemisphere_sites});
        
        data = [];
        data.label = arrayfun(@(x) sprintf('LFP%02d', x), 1:size(lfp_data_concat, 1), 'UniformOutput', false);
        data.fsample = sr;  % Sampling frequency
        data.time = {(0:size(lfp_data_concat,2)-1) / sr}; % Cell array containing time vector
        data.trial = {lfp_data_concat}; % Cell array containing LFP data matrix
        data.sampleinfo = [1 size(lfp_data_concat,2)]; % Start and end sample indices
        
        ft_datatype_raw(data) % To Check if the structure is correctly formatted
        data_orig = data;
        
%         cfg_ica              = [];
%         cfg_ica.resamplefs   = 256;
%         cfg_ica.detrend      = 'no';
%         data = ft_resampledata(cfg_ica, data_orig);
        fprintf("\n ============================================ \n investigating %s\n",session)
        triggers = sampling_table.site_triggers_ica(1).R_real;
        [ica_epochs, ica_avg] = ecg_bna_epochICAaroundTriggers(icaw_path,icawfile(1:end-12),data_to_plot, triggers, 0.25, 0.25, true,true);
       
        site_depth_idx = ismember({allSitesData.name},hemisphere_site_names');
        site_depth = arrayfun(@(x) x.site.electrode_depth, allSitesData(site_depth_idx), 'UniformOutput', false);
        site_depth = horzcat(site_depth{:});
        depths = site_depth;
        n_channels = length(data_to_plot.label);
        % ICA weight per Site plots
        fig1 = figure;
        for cmp = 1:n_channels
            
            comp_index = cmp;  
            weights = data_to_plot.unmixing(comp_index,:)';
            % Normalize for color mapping
            weights_norm = (weights - min(weights)) / (max(weights) - min(weights)); 
            subplot(6,6,cmp)
            scatter(weights, -depths, 100, weights_norm, 'filled');  % invert Y for depth
            colormap(jet);
            colorbar;
            xlabel('ICA weight');
            ylabel('Electrode Depth');
            title(sprintf('ICA Component %d Spatial Weights (LFP Channels)', comp_index),'FontSize',7);
            set(gca, 'YDir', 'normal');
            grid on;
        end
        fig1.WindowState = 'maximized';
        sgtitle(icawfile(1:end-12),'interpreter','none')
        figname2save = fullfile([icaw_path,filesep,icawfile(1:end-4)]);
        saveas(fig1,figname2save,'fig');
        saveas(fig1,figname2save,'png');
        
        % ICA Inv weight per Site plots
        weights_inv = pinv(data_to_plot.unmixing);
        fig2 = figure;
        for cmp = 1:n_channels
            
            comp_index = cmp;          
            weights_inv_2plot = weights_inv(:,comp_index)';
            weights_inv_norm = (weights_inv_2plot - min(weights_inv_2plot)) / (max(weights_inv_2plot) - min(weights_inv_2plot));
            
            subplot(6,6,cmp)
            scatter(weights_inv_2plot, -depths, 100, weights_inv_norm, 'filled');  % invert Y for depth 
            colormap(jet);
            colorbar;
            xlabel('ICA Inverse weight');
            ylabel('Electrode Depth');
%             ylabel('ICA cmp');
%             title(sprintf('Inverse ICA Component in Site: %d ', comp_index),'FontSize',7);
            title(sprintf('Inv ICA Component %d ', comp_index),'FontSize',7);
            set(gca, 'YDir', 'normal');
            grid on;
        end
        fig2.WindowState = 'maximized';
        sgtitle([icawfile(1:end-12),'_Inv'],'interpreter','none')
        figname2save = fullfile([icaw_path,filesep,icawfile(1:end-4),'_Inv2']);
        saveas(fig2,figname2save,'fig');
        saveas(fig2,figname2save,'png');
        
        
        if doTimecourse ==1
            cfg_ica                 = [];
            cfg_ica.continuous      = 'yes';
            cfg_ica.blocksize       = 30; % show 30 seconds at the time
            cfg_ica.plotevents      = 'no';
            cfg_ica.preproc.demean  = 'no';
            ica_com = ft_databrowser(cfg_ica, data_to_plot);
        else
            ica_com = data_to_plot;
            ica_com.topolabel = data_to_plot.cfg.channel';
            ica_com.topo = pinv(data_to_plot.unmixing);
        end
        
        BadICAcomps = inputdlg('Enter space-separated Selected ICA Components to remove:','Sample', [1 64]);
        BadICAcomps = str2num(BadICAcomps{1});
        close all
        
        cfg_ica = [];
        cfg_ica.component = BadICAcomps; % to be removed
        data_orig_clean = ft_rejectcomponent(cfg_ica, ica_com, data_orig);
        
        for sitei = 1:length(find(hemisphere_sites))
            site_idx = find(hemisphere_sites);
            allSitesData(site_idx(sitei)).site.LFP = data_orig_clean.trial{1, 1}(site_idx(sitei),:);
        end
        
    end
    
end

%% okay here the point was to average across all sites in one hemisphere, not in one target only..

if isfield(cfg.lfp, 'Reref') && cfg.lfp.Reref==1
    % i have never used contains, bu ti guess i would do something like this:
    % sites_hemisphere is supposed to be 1 for left hemispheres and 2 for
    % right hemispheres
    sites_hemisphere=contains(lfp_tagets,'_L')+2*contains(lfp_tagets,'_R'); %look for the _ as well, because f.e. VPL contains L
    hemisphere_list = unique(sites_hemisphere);
    
    for h = 1: length(hemisphere_list)
        %target = target_list{tr};
        hemisphere_sites = find(sites_hemisphere==hemisphere_list(h));
        
        valid_blocks = [];
        invalid_blocks = [];
        for b=1:numel(blocks_with_LFP)
            B=blocks_with_LFP(b);
            %% thing is i dont trust this as it is im afraid (sampling table and temp_allSites have very different indexing i believe?)
            %% -> hence the removing noise sites from sampling table beforehand
            
            %% i tend to create the indexes i will need inside this loop at the start so it doesnt get confusing
            block_idx=ismember(sampling_table.block,B) & ismember(sampling_table.site,hemisphere_sites);
            site_with_blockB_idx = sampling_table.site(block_idx);
            chan_with_BlockB_idx = cell2mat(arrayfun(@(x) x.site.channel, allSitesData(site_with_blockB_idx), 'UniformOutput', false));
            %site_with_blockB_idx = site_with_blockB_idx([temp_allSites(site_with_blockB_idx).noisy_site]==0);
            %             if (strcmp(cfg.session_info.Monkey  ,'Bacchus') && length(site_with_blockB_idx)>2 && (length(site_with_blockB_idx)<11)&& length(unique(chan_with_BlockB_idx))>=3)...
            %                     || (strcmp(cfg.session_info.Monkey  ,'Magnus') && length(site_with_blockB_idx)>2 && (length(site_with_blockB_idx)<33) && length(unique(chan_with_BlockB_idx))>=3)
            if (length(chan_with_BlockB_idx) > length(unique(chan_with_BlockB_idx)))
                warning('This block has 2 sites from the same channel : Block = %d', B );
            elseif length(unique(chan_with_BlockB_idx))>=3
                valid_blocks = [valid_blocks,B];
                block_lfp_concat = [];
                for s = 1: length(site_with_blockB_idx)
                    allSitesData(site_with_blockB_idx(s)).valid_blocks = valid_blocks;
                    allSitesData(site_with_blockB_idx(s)).block_with_not_enough_ref_elec = invalid_blocks;
                    %% i believe its easier to track like this (just a matter of taste i guess)
                    idx=block_idx & sampling_table.site==site_with_blockB_idx(s);
                    start_sample = sampling_table.start_sample(idx);
                    end_sample   = sampling_table.end_sample(idx);
                    block_lfp_concat = [block_lfp_concat; allSitesData(site_with_blockB_idx(s)).site.LFP(start_sample:end_sample)];
                end
                block_lfp_concat_mean = squeeze(mean(block_lfp_concat,1));
                
                for s = 1: length(site_with_blockB_idx)
                    %% i believe its easier to track like this (just a matter of taste i guess)
                    idx=block_idx & sampling_table.site==site_with_blockB_idx(s);
                    start_sample = sampling_table.start_sample(idx);
                    end_sample   = sampling_table.end_sample(idx);
                    allSitesData(site_with_blockB_idx(s)).site.LFP(start_sample:end_sample) = allSitesData(site_with_blockB_idx(s)).site.LFP(start_sample:end_sample) - block_lfp_concat_mean;
                end
            elseif length(unique(chan_with_BlockB_idx))<3
                invalid_blocks = [invalid_blocks ,B];
                for s = 1: length(site_with_blockB_idx)
                    allSitesData(site_with_blockB_idx(s)).block_with_not_enough_ref_elec = invalid_blocks;
                end
                continue;
            end
        end
    end
end

