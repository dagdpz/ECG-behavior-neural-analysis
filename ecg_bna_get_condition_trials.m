function cond_trials = ecg_bna_get_condition_trials(trials, condition)
%lfp_tfa_get_condition_trials - Function to get the indices of trials which
%satisfy a given trial condition
%
% USAGE:
%	cond_trials = lfp_tfa_get_condition_trials(site_lfp, condition)
%
% INPUTS:
%       site_lfp    - struct containing the condition
%       information about all trials for one site
%       condition   - struct containing the condition to be
%       analysed
%           Required fields:
%               type                - (integer) trial type
%               effector            - (integer) trial effector
%               choice              - choice (1) or instructed (1) trial
%               perturbation        - pre-injection (0) or post-injection (1)
%               perturbation_groups - blocks to be analysed for the given
%               perturbation, can be integer, integer array, 'all',
%               'allbutone'
%
% OUTPUTS:
%		cond_trials     - array of length (1 x N, N = number of
%		trials), with ones at indices of trials belonging to the given
%		condition
%
% REQUIRES:
%
% See also settings/lfp_tfa_settings, lfp_tfa_compare_conditions,
% lfp_tfa_site_average_tfr, lfp_tfa_site_evoked_LFP,
% lfp_tfa_site_powspctrum
%
% Author(s):	S.Nair, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2019-02-15:	Created function (Sarath Nair)
% 2019-03-05:	First Revision
% ...
% $Revision: 1.0 $  $Date: 2019-03-05 17:18:00 $

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%

cond_trials = ones(1, length(trials));
FN=fieldnames(condition);

for f=1:numel(FN)
    fn=FN{f};
    cond_trials = cond_trials & ismember([trials.(fn)],condition.(fn));
end