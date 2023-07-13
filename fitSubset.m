function [mdl, Rsqadjusted,rr,r0] = fitSubset(PSTH_f, predictorInfo, tgtGroups, param)

groups = [];
groups{1} = 1:predictorInfo.npredVars(1);
for iii=2:predictorInfo.nPredictors
    groups{iii} = groups{iii-1}(end)+(1:predictorInfo.npredVars(iii));
end

lagRanges = []; %lag ranges for all variables
for igroup = 1:numel(tgtGroups)
    lagRanges = cat(1, lagRanges, repmat(param.lagRange(tgtGroups(igroup),:),[numel(groups{tgtGroups(igroup)}) 1]));
end

tgtVars = cat(2, groups{tgtGroups});

y = PSTH_f;
[X, tlags] = getPredictorsDelayed(predictorInfo.t_r,...
    predictorInfo.predictors_r(tgtVars,:), lagRanges, predictorInfo.npredVars(tgtGroups),...
    param.predictorNames(tgtGroups));

%% fitting. use fitlm as it directly produces Rsq adjusted
mdl = fitlm(X,y);

beta = mdl.Coefficients.Estimate(2:end);

nVar = numel(tgtVars);
nLags = numel(tlags);

rr = reshape(beta,nLags,nVar);
r0 = mdl.Coefficients.Estimate(1);

Rsqadjusted = mdl.Rsquared.Adjusted;

