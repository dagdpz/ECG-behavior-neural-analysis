function [normalized] = ecg_bna_compute_shufflePredictor_normalization_general(real,shuffled,cfg)
% ecg_bna_compute_shufflePredictor_normalization_general - normalizing the real tfs
% and evoked data based on the shuffle predictor results
%
% USAGE:
%	[normalized] = ecg_bna_compute_shufflePredictor_normalization(variable,ecg_bna_cfg)
%
% INPUTS:
%		real            - 1xN struct containing the real data
%		shuffled  	    - 1xN struct containing the shuffled data
%       ecg_bna_cfg     - struct containing configuration settings
%
% OUTPUTS:
%       normalized        - 1xN struct containing the normalized results of
%       pow, itpc, lfp, or itpcbp data
% ======================================================================= %


method = cfg.lfp.normalization;

switch cfg.to_trigger
    case 'LFP';
        parameters={'pow','itpc','lfp','itpcbp','powbp','pha'};
    case 'MUA';
        parameters ={'mua'};
    case 'ECG';
        parameters ={'ecg'};
end

for p=1:numel(parameters)
    parameter=parameters{p};
    if ismember(parameter,{'lfp','ecg','mua'})
        realmean = real.(parameter).complete;
        realstd = real.(parameter).std;
        shuffledmean=shuffled.(parameter).mean;
        shuffledstd=shuffled.(parameter).std;
        if strcmp(method , 'subtraction')
            normalized.(parameter).mean    = mean((realmean-shuffledmean),1);
            normalized.(parameter).std     = std((realmean-shuffledmean),0,1);
            normalized.(parameter).sterr   = sterr((realmean-shuffledmean),1);
        elseif strcmp(method , 'division')
            normalized.(parameter).mean    = mean((realmean./shuffledmean),1);
            normalized.(parameter).std     = std((realmean./shuffledmean),0,1);% realstd;%./shuffledmean;  %%??
            normalized.(parameter).sterr   = sterr((realmean./shuffledmean),1);
        elseif strcmp(method , 'zscore')
            normalized.(parameter).mean    = mean(((realmean-shuffledmean)./shuffledstd),1);
            normalized.(parameter).std     = std(((realmean-shuffledmean)./shuffledstd),0,1);% realstd;%./shuffledstd; %% ??
            normalized.(parameter).sterr   = sterr(((realmean-shuffledmean)./shuffledstd),1);
        elseif strcmp(method , 'not normalized')
            normalized.(parameter).mean    = realmean;
            normalized.(parameter).std    = realstd;
        end
    else
        realmean=real.(parameter).mean;
        realstd=real.(parameter).std;
        shuffledmean=shuffled.(parameter).mean;
        shuffledstd=shuffled.(parameter).std;
        if strcmp(method , 'subtraction')
            normalized.(parameter).mean    = realmean-shuffledmean;
            normalized.(parameter).std    = shuffledstd;  %%??
        elseif strcmp(method , 'division')
            normalized.(parameter).mean    = realmean./shuffledmean;
            normalized.(parameter).std    = realstd;%./shuffledmean;  %%??
        elseif strcmp(method , 'zscore')
            normalized.(parameter).mean    = (realmean-shuffledmean)./shuffledstd;
            normalized.(parameter).std    = realstd;%./shuffledstd; %% ??
        elseif strcmp(method , 'not normalized')
            normalized.(parameter).mean    = realmean;
            normalized.(parameter).std    = realstd;
        end
    end
end

