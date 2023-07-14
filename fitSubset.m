function [Rsqadjusted,rr,r0] = fitSubset(PSTH_f, predictorInfo, tgtGroups, param)
%[Rsqadjusted,rr,r0] = fitSubset(PSTH_f, predictorInfo, tgtGroups, param)

fitMethod = 'ridge';
useGPU = 0;

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

        nVar = numel(tgtVars);
        nLags = numel(tlags);

switch fitMethod
    case 'fitlm'


        %fitting. use fitlm as it directly produces Rsq adjusted
        mdl = fitlm(X,y);

        beta = mdl.Coefficients.Estimate(2:end);

        rr = reshape(beta,nLags,nVar);
        r0 = mdl.Coefficients.Estimate(1);

        Rsqadjusted = mdl.Rsquared.Adjusted;

    case 'ridge'
            b = rReg(X, y, param.ridgeParams, useGPU); %from ridgeXs.m

            r0 = b(1);
            rr = reshape(b(2:end),nLags,nVar);
            Xc = reshape(X, size(X,1), nLags, nVar);
            yHat = r0;
            for ibhv = 1:nVar
                yHat = yHat + Xc(:,:,ibhv)*rr(:,ibhv);
            end

         % [rr, r0, yHat] = ridgeXs(tpredictorInfo.t_r,...
         %    predictorInfo.predictors_r(tgtVars,:),  y, ...
         % param.lagRange(tgtVars,:), param.ridgeParam);
          
         nPredictors = numel(rr)+1;
         [Rsqadjusted] = getRsqadj(y, yHat, nPredictors);
end