function  kernelInfo_norm = getKernelInfo_norm(kernelInfo, predictorInfo)

kernelInfo_norm = kernelInfo;
npredVars = predictorInfo.npredVars;
predVars = 0;
for ivar = 1:numel(npredVars)
    predVars = predVars(end)+1:predVars(end)+predictorInfo.npredVars(ivar);
    kernelInfo_norm.kernel{ivar} = kernelInfo.kernel{ivar} ./ nanmean(predictorInfo.predictors_r(predVars,:),2)';
end
