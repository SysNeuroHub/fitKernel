function [expval_tgt, corr_tgt] = getExpVal_avgtgt(PSTH_f, predicted, ...
    catEvTimes, t_r, tWin, cardinalDir, dd)
%  [expval_tgt, corr_tgt] = getExpVal_avgtgt(PSTH_f, predicted, ...
%     catEvTimes, t_r, tWin, cardinalDir, dd)

%% explained variance for target response
onsetTimes = catEvTimes.tOnset;
validEvents = find(~isnan(onsetTimes));
onsetTimes = onsetTimes(validEvents);

tgtDir = getTgtDir(dd.targetloc(validEvents), cardinalDir);

[~,dirIdx]=intersect(cardinalDir, unique(tgtDir));

y_r = cat(2,PSTH_f,predicted);

[avgTonset(dirIdx,:,:), winSamps_tonset] ...
    = eventLockedAvg(y_r', t_r, onsetTimes, tgtDir, tWin);
%avgTonsetByCue: direction x variable x time

expval_tgt = zeros(size(predicted,2),1);
mse_tgt = zeros(size(predicted,2),1);
corr_tgt = zeros(size(predicted,2),1);
for ivar = 1:size(predicted,2)
    
    observedAvg = avgTonset(:,1,:);
    observedAvg = observedAvg(:);
    predictedAvg = avgTonset(:,ivar+1,:);
    predictedAvg = predictedAvg(:);
    
    [expval_tgt(ivar,1), mse_tgt(ivar,1), corr_tgt(ivar,1)] = ...
        getExpVal(observedAvg, predictedAvg);
end