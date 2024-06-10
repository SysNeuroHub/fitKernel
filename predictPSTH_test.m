function [predicted_all, predicted] = predictPSTH_test(predictorInfo, param, kernelInfo)

% from fitPSTH.m only applicable to option==5
predicted = zeros(numel(predictorInfo.npredVars), length(predictorInfo.t_r));
dt_r = median(diff(t_r));
for ivar = 1:numel(predictorInfo.npredVars)
    if ivar==1
        theseVarIdx = 1:predictorInfo.npredVars(1);
    else
        theseVarIdx = sum(predictorInfo.npredVars(1:ivar-1))+1:sum(predictorInfo.npredVars(1:ivar));
    end
    if size(param.lagRange,1)== 1
        thisLagRange = param.lagRange;
    else
        % thisLagRange = [min(param.lagRange(:,1)) max(param.lagRange(:,2))];
        thisLagRange = param.lagRange(ivar,:);
    end

    predicted(ivar, :) = 1/dt_r * predictXs(predictorInfo.t_r, predictorInfo.predictors_r(theseVarIdx,:), ...
        kernelInfo.intercept, kernelInfo.kernel{ivar}, thisLagRange);
end
predicted_all = sum(predicted-1/dt_r * kernelInfo.intercept,1) + 1/dt_r * kernelInfo.intercept;