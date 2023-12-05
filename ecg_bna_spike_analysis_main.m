function ecg_bna_spike_analysis_main(cfg)
% This is a function that runs all the spike analysis related functions
% together.

%% setting section - hopefully temporary
compute_unit_subsets      = 0;
move_files                = 0;

compute_spike_histograms  = 1;
plot_spike_histograms     = 1;
compute_spike_phase       = 1;
plot_spike_phase          = 1;

population_analysis       = 1;


if move_files
    %% Temporary - copy selected units separately
    ecg_bna_copy_selected_units(cfg)
end
% cfg.unit_IDS=keys.tuning_table(2:end,1);

%% Get info about sessions to be analysed
% Read the info about sessions to analyse
sessions_info = cfg.session_info;

%% per session processing..
for i = 1:length(sessions_info)
    java.lang.System.gc() % added by Luba to control the number of graphical device interface handles (hope to handle the problem with freezing plots while creating figures)
    
    session_name = [sessions_info(i).Monkey '_' sessions_info(i).Date];
    
    % First make seed and ecg shuffles, then use those shuffles for all subfunctions
    seed_filename=[cfg.ECG_root_results_fldr filesep 'seed.mat']; %% not quite sure yet where to put seed
    if exist(seed_filename,'file')
        load(seed_filename);
        rng(seed);
    else
        seed=rng;
        save(seed_filename,'seed');
    end
    
    if compute_spike_histograms
        Rpeaks=ecg_bna_compute_session_shuffled_Rpeaks(sessions_info(i),cfg);
        ecg_bna_compute_session_spike_histogram(sessions_info(i),Rpeaks,cfg);
    end
    
    if plot_spike_histograms
        ecg_bna_plot_session_spike_histogram(sessions_info(i),cfg);
    end
    
    if compute_spike_phase
        Rpeaks=ecg_bna_compute_session_shuffled_Rpeaks(sessions_info(i),cfg);
        ecg_bna_compute_session_ECG_related_spikePhase(sessions_info(i),Rpeaks,cfg)
    end
    
    if plot_spike_phase
        ecg_bna_plot_session_ECG_related_spikePhase(sessions_info(i),cfg)
    end
end

%% average across sessions
if population_analysis
    SPK_PSTH=load_spikes(sessions_info,'SPK_fldr','','','per_unit','Output');
    ecg_bna_avg_spike_histogram(SPK_PSTH,sessions_info, cfg);
end

end

function Out=load_spikes(sessions_info,subfield,subfolder,namepart,per,varname)
condition_labels = {'Rest', 'Task'};
Out = cell(length(sessions_info), 1);
for i = 1:length(sessions_info)
    disp(['Session ' num2str(i) ' out of ' num2str(length(sessions_info))])
    if isempty(subfolder) %% unfortunate inconsistent naming
        monkey=sessions_info(i).Monkey(1:3);
    else
        monkey=['_' sessions_info(i).Monkey(1:3)];
    end
    warning('There is a hardcoded filepath')
%     results_folder = fullfile(sessions_info(i).(subfield),per,subfolder);
    results_folder = 'Y:\Projects\Pulv_bodysignal\ECG_triggered_spikes\ECG_TaskRest_Magnus_merged\per_unit_selected_stableTaskAndRest_600';
    to_load=dir(fullfile(results_folder, [namepart monkey '_' sessions_info(i).Date '*.mat']));
    if ~isempty(to_load)
        tmp = arrayfun(@(x) load([x.folder filesep x.name], varname), to_load, 'UniformOutput', false);
        
        % put common parameters in the structure
        Out{i}.unit_ID = cellfun(@(x) x.Output.unit_ID, tmp, 'UniformOutput', false);
        Out{i}.target = cellfun(@(x) x.Output.target, tmp, 'UniformOutput', false);
        Out{i}.quantSNR = cellfun(@(x) x.Output.quantSNR, tmp);
        Out{i}.Single_rating = cellfun(@(x) x.Output.Single_rating, tmp);
        Out{i}.stability_rating = cellfun(@(x) x.Output.stability_rating, tmp);
        
        field_names = fieldnames(tmp{1}.Output.Rest);
        for conditionNum = 1:2
            for fieldNum = 1:length(field_names)
                if strcmp(field_names{fieldNum}, 'Rts') || ...
                        strcmp(field_names{fieldNum}, 'Rds') || ...
                        strcmp(field_names{fieldNum}, 'Rds_perm') || ...
                        strcmp(field_names{fieldNum}, 'raster') || ...
                        ischar(tmp{1}.Output.(condition_labels{conditionNum}).(field_names{fieldNum}))
                    Out{i}.(condition_labels{conditionNum}).(field_names{fieldNum}) = ...
                        cellfun(@(x) cat(1,x.Output.(condition_labels{conditionNum}).(field_names{fieldNum})), tmp, 'UniformOutput', false);
                else
                    tmp_out = cellfun(@(x) x.Output.(condition_labels{conditionNum}).(field_names{fieldNum}), tmp, 'UniformOutput', false);
                    tmp_out = cat(1, tmp_out{:});
                    Out{i}.(condition_labels{conditionNum}).(field_names{fieldNum}) = tmp_out;
                    clear tmp_out
                end
            end
        end
    end
end
Out = Out(~cellfun('isempty',Out));
end