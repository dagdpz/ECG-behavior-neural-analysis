function [significance] = ecg_bna_compute_significance(real,shuffled,ecg_bna_cfg)
% ecg_bna_compute_significance - checking statistical significance of the real tfs
% and evoked data based on the shuffle predictor results
%
% USAGE:
%	[significance] = ecg_bna_compute_significance(real,shuffled,ecg_bna_cfg)
%
% INPUTS:
%		real            - 1xN struct containing the real data
%		shuffled  	    - 1xN struct containing the shuffled data
%       ecg_bna_cfg     - struct containing configuration settings
%		
% OUTPUTS: 
%       significance        - 1xN struct containing the statistical analysis 
%       results of pow, itpc, lfp, or itpcbp data
% ======================================================================= %

method = ecg_bna_cfg.significance_method;
parameters={'pow','itpc','lfp','itpcbp'};
for p=1:numel(parameters)
    parameter=parameters{p};
    realmean=real.(parameter).mean;
%     shuffledmean=shuffled.(parameter).mean;
%     shuffledstd=shuffled.(parameter).std;
    if strcmp(method , '95Conf_intrvl')
        suffledConf = shuffled.(parameter).conf95;
%         significance.(parameter).upper = find(realmean > suffledConf(1,:,:));
%         significance.(parameter).lower = find(realmean < suffledConf(2,:,:));
        significance.(parameter) = (suffledConf(1,:,:) < realmean < suffledConf(2,:,:));

    end
    
end  