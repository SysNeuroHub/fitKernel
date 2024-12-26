function [spkOkTrials, spkOk_th] = getSpkOKtrials(spk_all_cat, t_r, trIdx_r, param)
%[spkOkTrials, spkOk_th] = getSpkOKtrials(spk_all_cat, t_r, trIdx_r, param)

% from fitPSTH_cv
detrend = 1; %22/7/22

dt_r = median(diff(t_r));

PSTH_r = getPSTH(spk_all_cat, t_r);
PSTH_f = filtPSTH(PSTH_r, dt_r, param.psth_sigma, 2, detrend);

spkOk_th = param.spkOk_th*mad(PSTH_f) + median(PSTH_f);

ng_idx = find(PSTH_f > spkOk_th);


%% only retrain trials which firing rate keeps below the threshold
spkOkTrials = [];
for itr = 1:numel(trIdx_r)
    if isempty(intersect(trIdx_r{itr}, ng_idx))
        spkOkTrials = [spkOkTrials itr];
    end
end

