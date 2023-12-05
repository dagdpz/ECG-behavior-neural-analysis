function cond_trials = ecg_bna_get_condition_trials(trials, condition)
% INPUTS: condition from ecg_bna_cfg.condition
%         trials from trials structer
% OUTPUT: logical of size(trials) with given conditions;

cond_trials = ones(1, length(trials));
FN=fieldnames(condition);

for f=1:numel(FN)
    fn=FN{f};
    if isfield(trials,fn)
        cond_trials = cond_trials & ismember([trials.(fn)],condition.(fn));
    end
end