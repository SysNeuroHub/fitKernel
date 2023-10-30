function [predictorInfo_split, param_split] = splitPredictorByCue(predictorInfo, dd, onsets_cat, param)

[rewardTimes, punishTimes] = getRewardTimes(dd);
trialEndTimes = sort([rewardTimes punishTimes]);

%from getCueDirMtx.m
eventTimes = [];
for itr = 1:dd.numTrials
    %only register trials with cue and feedback (success or failure)
    cueOnset = onsets_cat.cueOnset(itr);
    cOnset = onsets_cat.cOnset(itr);
    if isinf(cueOnset) || isnan(cOnset) || isnan(trialEndTimes(itr))
        continue;
    end

    eventTimes = [eventTimes; cueOnset trialEndTimes(itr)];
end

wCueTrace = event2Trace(predictorInfo.t_r, eventTimes)';
woCueTrace = 1-wCueTrace;

%% split predictor into trials with cue and trials without cue
predictorInfo_split = predictorInfo;
predictorInfo_split.npredVars = [predictorInfo.npredVars; predictorInfo.npredVars];
predictorInfo_split.nPredictors = predictorInfo.nPredictors * 2;
nRows_ori = size(predictorInfo.predictors_r,1);
predictorInfo_split.predictors_r(1:nRows_ori,:) = predictorInfo.predictors_r .* woCueTrace;
predictorInfo_split.predictors_r(1+nRows_ori:2*nRows_ori,:) = predictorInfo.predictors_r .* wCueTrace;

%% modify param
param_split = param;
predictorNames_ori = param.predictorNames;
param_split.predictorNames = [cellfun(@(x)([x '_wocue']),predictorNames_ori,'UniformOutput',false),...
    cellfun(@(x)([x '_wcue']),predictorNames_ori,'UniformOutput',false)];
param_split.lagRange = [param.lagRange; param.lagRange];


