function [includeTrace, includeTrials] = getIncludeTrace(t_cat, t_r, t_tr, onsets_cat, spkOkTrials)
% obtain predictor only for trials without cue
%
% 5/6/24: NEED FIX woCue condition includes outside trials
%6/12/24 created from splitPredictorByCue

%trialEndTimes = getRewardTimes(dd);

[trIdx_r] = retrieveTrIdx_r(t_cat, t_r, t_tr);

%from getCueDirMtx.m
eventTimes = []; excludeTrials = [];
for itr = 1:numel(t_tr)
    cueOnset = onsets_cat.cueOnset(itr);
    % cOnset = onsets_cat.cOnset(itr);
    % if  isnan(cOnset) || isnan(trialEndTimes(itr))
    %     continue;
    % end

    % if ~isinf(cueOnset) 
    %     eventTimes = [eventTimes; cueOnset t_r(trIdx_r{itr}(end))];
    %     excludeTrials = [excludeTrials itr];
    % end
    if  sum(ismember(spkOkTrials, itr)) == 0 || ~isinf(cueOnset) 
        eventTimes = [eventTimes; t_r( trIdx_r{itr}(1)) t_r(min(trIdx_r{itr}(end)+1, numel(t_r)))];
        excludeTrials = [excludeTrials itr];
    end
end

excludeTrace = event2Trace(t_r, eventTimes);
includeTrace = 1-excludeTrace;
includeTrials = setdiff(1:numel(t_tr), excludeTrials);

