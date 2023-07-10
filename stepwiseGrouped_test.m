detrend = 1; %22/7/22
t_r = predictorInfo.t_r;
sigma = param.psth_sigma;

dt_r = median(diff(t_r));
PSTH_r = getPSTH(spk_all_cat, t_r);
PSTH_f = filtPSTH(PSTH_r, dt_r, sigma, 2, detrend);

y = PSTH_f;
X = predictorInfo.predictors_r';
groups = [];
groups{1} = 1:predictorInfo.npredVars(1);
for ii=2:predictorInfo.nPredictors
    groups{ii} = groups{ii-1}(end)+(1:predictorInfo.npredVars(ii));
end

% super slow
[model, selectedGroups] = stepwiseGrouped(X, y, groups);


[ynew,ynewci] = predict(model,X);
