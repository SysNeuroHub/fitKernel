function predictorInfo = preparePredictors(dd, eyeData_rmotl_cat, t_r, param, catEvTimes)
% predictorInfo = preparePredictors(dd, eyeData_rmotl_cat, t_r, param, catEvTimes)
% returns matrix (time x target direction) about predictors specified as param.predictorNames
%
%param:
%     cardinalDir
%     dt_r
%     cutoffFreq
%     predictorNames:
%     {'vision','cue','eyeposition','pdiam','pdiam_prctile','parea','blink','saccade','reward'}
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
            thisPredictor = getTgtDirMtx(dd, t_r, catEvTimes, param.cardinalDir);%,...
                %eyeData_rmotl_cat); %< currently ignoring trials with NAN
        case 'cue' %NOT YET IMPLEMENTED
            thisPredictor = getCueDirMtx(dd, t_r, catEvTimes, param.cardinalDir);
        case 'eyeposition'
            thisPredictor = getEyeDirMtx(dd, t_r, eyeData_rmotl_cat, param.cardinalDir);
        case 'pdiam'
            [~, pdiam] = getPupilDiameter(eyeData_rmotl_cat);
            pdiam_r = interp1(eyeData_rmotl_cat.t, pdiam, t_r)';
            thisPredictor = hpFilt(pdiam_r', 1/param.dt_r, param.cutoffFreq)';
            thisPredictor = lpFilt(thisPredictor', 1/param.dt_r, 2)'; %test 11/6/24
        case 'pdiam_prctile'
            pdiam_prctile = getPupilDiameter(eyeData_rmotl_cat);
            thisPredictor = interp1(eyeData_rmotl_cat.t, pdiam_prctile, t_r)';
        case 'pdiamspd'
            [~, pdiam] = getPupilDiameter(eyeData_rmotl_cat);
            pdiam_r = interp1(eyeData_rmotl_cat.t, pdiam, t_r)';
            pdiam_r = hpFilt(pdiam_r', 1/param.dt_r, param.cutoffFreq)';
            thisPredictor = [diff(pdiam_r) 0]/param.dt_r; %too noisy
        case 'parea'
            parea_r = interp1(eyeData_rmotl_cat.t, eyeData_rmotl_cat.parea, t_r)';
            thisPredictor = hpFilt(parea_r', 1/param.dt_r, param.cutoffFreq)';
        case 'blink'
            blinks = event2Trace(eyeData_rmotl_cat.t, [catEvTimes.blinkStartTimes catEvTimes.blinkEndTimes]);
            thisPredictor = interp1(eyeData_rmotl_cat.t, single(blinks), t_r)';
        case 'eyespeed' %resampling rate needs to be high (>100Hz)
            ex1 = event2Trace(eyeData_rmotl_cat.t,[catEvTimes.blinkStartTimes catEvTimes.blinkEndTimes]);
            ex2 = event2Trace(eyeData_rmotl_cat.t,[catEvTimes.outlierStartTimes catEvTimes.outlierEndTimes]);
            excludePeriod = single(ex1+ex2>0);
            excludePeriod_r = logical(interp1(eyeData_rmotl_cat.t, excludePeriod, t_r, 'nearest'));

            thisPredictor = getEyeSpdDirMtx(dd, t_r, eyeData_rmotl_cat, ...
                param.cardinalDir, excludePeriod_r);
        case 'saccade'
            % exculde blink and outlier periods
            %             excludePeriod = (event2Trace(eyeData_rmotl_cat.t, [catEvTimes.blinkStartTimes catEvTimes.blinkEndTimes; ...
            %                 catEvTimes.outlierStartTimes
            %                 catEvTimes.outlierEndTimes])>0); %NG!!
            ex1 = event2Trace(eyeData_rmotl_cat.t,[catEvTimes.blinkStartTimes catEvTimes.blinkEndTimes]);
            ex2 = event2Trace(eyeData_rmotl_cat.t,[catEvTimes.outlierStartTimes catEvTimes.outlierEndTimes]);
            excludePeriod = (ex1+ex2>0);
            
            [startSacc, endSacc] = selectSaccades(catEvTimes.saccadeStartTimes, ...
                catEvTimes.saccadeEndTimes, eyeData_rmotl_cat.t, excludePeriod);
            thisPredictor = getSaccMtx(t_r, startSacc, endSacc, eyeData_rmotl_cat, param.cardinalDir);
        case 'reward'
            [rewardTimes, punishTimes] = getRewardTimes(dd);
            
            thisPredictor = cat(1,event2Trace(t_r, rewardTimes)',event2Trace(t_r, punishTimes)');
        case 'fixationbhv'
            valid = ~isnan(catEvTimes.fOnset);
            thisPredictor = event2Trace(t_r, catEvTimes.fOnset(valid))';
    end
    predictors_r = cat(1, predictors_r, thisPredictor);
    npredVars(ivar) = size(thisPredictor,1);
end
%predictors_r = norm_std_mean(predictors_r')';23/1/22

predictorInfo.predictors_r = predictors_r;
predictorInfo.npredVars = npredVars;
predictorInfo.t_r = t_r;
predictorInfo.nPredictors = nPredictors;
