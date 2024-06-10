function [predictorInfo_limit, param_limit] = limitPredictor(predictorInfo, dd, onsets_cat, param, option)
% created from splitPredictorByCue
% option 1: only use period with successful trials
% option 2: use period EXCEPT successful trials

[rewardTimes,punishTimes,~, trialOutcome] = getRewardTimes(dd);
trialEndTimes = sort([rewardTimes punishTimes]);

%from getCueDirMtx.m
eventTimes = [];

for itr = 1:dd.numTrials
    %only register trials with reward 
    fOnset = onsets_cat.fOnset(itr);
    cOnset = onsets_cat.cOnset(itr);
    if isnan(cOnset) || isnan(fOnset) || trialOutcome(itr) == -1
        continue;
    end

    eventTimes = [eventTimes; fOnset trialEndTimes(itr)];
end

tgtPeriod = event2Trace(predictorInfo.t_r, eventTimes)';

if option == 2
    tgtPeriod = 1 - tgtPeriod;
end

%% split predictor into trials with cue and trials without cue
predictorInfo_limit = predictorInfo;
%predictorInfo_limit.npredVars = [predictorInfo.npredVars; predictorInfo.npredVars];
%predictorInfo_limit.nPredictors = predictorInfo.nPredictors * 2;
%nRows_ori = size(predictorInfo.predictors_r,1);
predictorInfo_limit.predictors_r = predictorInfo.predictors_r .* tgtPeriod;
% predictorInfo_split.predictors_r(1+nRows_ori:2*nRows_ori,:) = predictorInfo.predictors_r .* tgtPeriod;

%% modify param
param_limit = param;
predictorNames_ori = param.predictorNames;
param_limit.predictorNames = [cellfun(@(x)([x '_success']),predictorNames_ori,'UniformOutput',false)];
%param_limit.lagRange = [param.lagRange; param.lagRange];


