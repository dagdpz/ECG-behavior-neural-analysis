function [SD,RAST,PSTH_time,SD_1ms,RAST_1ms,PSTH_time_1ms]=ecg_bna_spike_density(AT,trial_onsets,trial_ends,cfg)
%% make SD across all trials appended (no average)!
switch cfg.kernel_type
    case 'gaussian'
        Kernel     = normpdf(-5*cfg.gaussian_kernel:cfg.PSTH_binwidth:5*cfg.gaussian_kernel,0,cfg.gaussian_kernel);
    case 'box'
        n_bins=round(2*cfg.gaussian_kernel/cfg.PSTH_binwidth);
        Kernel=ones(1,n_bins)/n_bins/cfg.PSTH_binwidth; %%*1000 cause a one full spike in one 1ms bin means 1000sp/s locally
end
PSTH_time     = trial_onsets(1):cfg.PSTH_binwidth:trial_ends(end);
PSTH_time_1ms = trial_onsets(1):0.001:trial_ends(end);
RAST          = hist(AT,PSTH_time);
RAST_1ms      = hist(AT,PSTH_time_1ms);
SD            = conv(RAST,Kernel,'same');
SD_1ms        = conv(RAST_1ms,Kernel,'same');
end
