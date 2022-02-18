function predictorInfo = preparePredictors(dd, eyeData_rmotl_cat, t_r, param, catEvTimes)
%param:
%     cardinalDir
%     dt_r
%     cutoffFreq
%
% predictorInfo:
%     predictors_r
%     t_r
%     npredVars
%     nPredictors



nPredictors = length(param.predictorNames);
predictors_r = [];
npredVars = zeros(nPredictors,1);
for ivar = 1:nPredictors
    disp(['preparePredictors: ' param.predictorNames{ivar}]);
    switch param.predictorNames{ivar}
        case 'vision'
            thisPredictor = getTgtDirMtx(dd, t_r, catEvTimes, param.cardinalDir);
        case 'cue' %NOT YET IMPLEMENTED
            thisPredictor = getCueDirMtx(dd, t_r, catEvTimes, param.cardinalDir);
        case 'eyeposition'
            thisPredictor = getEyeDirMtx(dd, t_r, eyeData_rmotl_cat, param.cardinalDir);
        case 'pdiam'
            [~, pdiam] = getPupilDiameter(eyeData_rmotl_cat);
            pdiam_r = interp1(eyeData_rmotl_cat.t, pdiam, t_r)';
            thisPredictor = hpFilt(pdiam_r', 1/param.dt_r, param.cutoffFreq)';
        case 'pdiam_prctile'
            pdiam_prctile = getPupilDiameter(eyeData_rmotl_cat);
            thisPredictor = interp1(eyeData_rmotl_cat.t, pdiam_prctile, t_r)';
        case 'parea'
            parea_r = interp1(eyeData_rmotl_cat.t, eyeData_rmotl_cat.parea, t_r)';
            thisPredictor = hpFilt(parea_r', 1/param.dt_r, param.cutoffFreq)';
        case 'blink'
            blinks = event2Trace(eyeData_rmotl_cat.t, [catEvTimes.blinkStartTimes catEvTimes.blinkEndTimes]);
            thisPredictor = interp1(eyeData_rmotl_cat.t, single(blinks), t_r)';
        case 'eyespeed' %resampling rate needs to be high (>100Hz)
            excludeTimes = [catEvTimes.blinkStartTimes catEvTimes.blinkEndTimes; ...
                catEvTimes.outlierStartTimes catEvTimes.outlierEndTimes];
            thisPredictor = getEyeSpdDirMtx(dd, t_r, eyeData_rmotl_cat, ...
                param.cardinalDir, excludeTimes);
        case 'saccade'
            % exculde blink and outlier periods
            excludePeriod = (event2Trace(eyeData_rmotl_cat.t, [catEvTimes.blinkStartTimes catEvTimes.blinkEndTimes; ...
                catEvTimes.outlierStartTimes catEvTimes.outlierEndTimes])>0);
            
            [startSacc, endSacc] = selectSaccades(catEvTimes.saccadeStartTimes, ...
                catEvTimes.saccadeEndTimes, eyeData_rmotl_cat.t, excludePeriod);
            thisPredictor = getSaccMtx(t_r, startSacc, endSacc, eyeData_rmotl_cat, param.cardinalDir);
        case 'reward'
            [rewardTimes, punishTimes] = getRewardTimes(dd);
            
            thisPredictor = cat(1,event2Trace(t_r, rewardTimes)',event2Trace(t_r, punishTimes)');
    end
    predictors_r = cat(1, predictors_r, thisPredictor);
    npredVars(ivar) = size(thisPredictor,1);
end
%predictors_r = norm_std_mean(predictors_r')';23/1/22

predictorInfo.predictors_r = predictors_r;
predictorInfo.npredVars = npredVars;
predictorInfo.t_r = t_r;
predictorInfo.nPredictors = nPredictors;
