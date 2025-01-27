function prefDir = getPrefDir_wrapper(y_r, t_r, dd, catEvTimes, param, includeTrials)

%preferred direction determined from includeTrials
%originally implemented in showTonsetResp.m

[choiceOutcome] = getChoiceOutcome(dd);

figTWin = param.figTWin; %figure temporal window
compTWin = [figTWin(1) - 0.2 figTWin(2)]; %triggered response is computed

%% triggered by tOnsets
successEvents = intersect(find((choiceOutcome==1)), includeTrials);
successEvents = intersect(successEvents, find(catEvTimes.tOnset + compTWin(2) < max(t_r)));
%only use trials when the choices were registered.
%this is a temporary fix as my current algorithm assumes stimuli were NOT
%presented, causing no visual response in the model

onsetTimes = catEvTimes.tOnset(successEvents);
tgtDir = getTgtDir(dd.targetloc(successEvents), param.cardinalDir);

%[~,dirIdx] = intersect(param.cardinalDir, unique(tgtDir));

% [~, winSamps, singleOnsetResp, ...
%     sortedOnsetLabels, uniqueOnsetLabels] ...
%     = eventLockedAvg(y_r', t_r, onsetTimes, tgtDir, compTWin);
% tonsetRespAmp = characteriseResp(singleOnsetResp, ...
%     winSamps, param.tOnRespWin, param.baseWin, 'mean');

 prefDir = getPrefDir(y_r(:,1), t_r, onsetTimes, tgtDir, param);