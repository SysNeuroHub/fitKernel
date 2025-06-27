function [corr_tgt_avg, corr_tgt_avg_rel] = getCorr_tgt_avg(t_r, y_r, catEvTimes, dd, param, includeTrials)
% [corr_tgt_avg, corr_tgt_avg_rel] = getCorr_tgt_avg(t_r, y_r, catEvTimes, dd, param, includeTrials)


compTWin = param.tOnRespWin; %triggered response is computed

[choiceOutcome] = getChoiceOutcome(dd);


%% triggered by tOnsets
successEvents = intersect(find((choiceOutcome==1)), includeTrials);
successEvents = intersect(successEvents, find(catEvTimes.tOnset + compTWin(2) < max(t_r)));
%only use trials when the choices were registered.
%this is a temporary fix as my current algorithm assumes stimuli were NOT
%presented, causing no visual response in the model

onsetTimes = catEvTimes.tOnset(successEvents);
tgtDir = getTgtDir(dd.targetloc(successEvents), param.cardinalDir);

[~,dirIdx] = intersect(param.cardinalDir, unique(tgtDir));

[~, winSamps, singleOnsetResp, ...
    sortedOnsetLabels, uniqueOnsetLabels] ...
    = eventLockedAvg(y_r', t_r, onsetTimes, tgtDir, compTWin);
tonsetRespAmp = characteriseResp(singleOnsetResp, ...
    winSamps, param.tOnRespWin, param.baseWin, 'mean');

% prefDir = getPrefDir(y_r(:,1), t_r, onsetTimes, tgtDir, param);
% prefDir = getPrefDir_wrapper(y_r, t_r, dd, catEvTimes, param, includeTrials);
% 
% tonsetRespAmp_p = characteriseResp(singleOnsetResp, ...
%     winSamps, param.tOnRespWin, [], 'mean');
% tonsetRespAmp_b = characteriseResp(singleOnsetResp, ...
%     winSamps, param.baseWin, [], 'mean');


%% select trials
% if allTr
%     theseTrials = 1:numel(tgtDir);
% else
% %    theseTrials = prefDirTrials;%trials with cell's preferred direction
%     theseTrials = find(tgtDir == prefDir);
% end

for iDir = 1: numel(param.cardinalDir)
    theseTrials = find(tgtDir == param.cardinalDir(iDir));
    mtOnsetResp(:,:,iDir) = squeeze(mean(singleOnsetResp(theseTrials,:,:)));%avg response to preferred direction
end

noNanDirs = find(~isnan(squeeze(sum(sum(mtOnsetResp,1),2))));
mtOnsetResp_c = reshape(mtOnsetResp(:,:,noNanDirs), size(mtOnsetResp,1), size(mtOnsetResp,2)*numel(noNanDirs));
nVars = size(y_r, 2);
corr_tgt_avg = zeros(nVars,1);
for ivar = 1:nVars-1
 
    corr_tgt_avg(ivar,1) = corr(mtOnsetResp_c(1,:)', mtOnsetResp_c(ivar+1,:)');
end
corr_tgt_avg_rel = corr_tgt_avg./corr_tgt_avg(1,1);
