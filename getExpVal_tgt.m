function [expval_tgt, corr_tgt] = getExpVal_tgt(PSTH_f, predicted, catEvTimes, t_r, tWin)
% expval_tgt = getExpVal_tgt(PSTH_f, predicted, catEvTimes, t_r, tWin)

%% explained variance for target response
idxTgtOnsets = [];
onsetTimes = catEvTimes.tOnset;
onsetTimes = onsetTimes(~isnan(onsetTimes));
for ievent = 1:numel(onsetTimes)
    idxTgtOnsets = cat(1, idxTgtOnsets, ...
        find((t_r>=onsetTimes(ievent)+tWin(1)) .* (t_r<onsetTimes(ievent)+tWin(2))));
end

expval_tgt = zeros(size(predicted,2),1);
corr_tgt = zeros(size(predicted,2),1);
for ivar = 1:size(predicted,2)
    expval_tgt(ivar,1) = getExpVal(PSTH_f(idxTgtOnsets)-mean(PSTH_f(idxTgtOnsets)), ...
        predicted(idxTgtOnsets,ivar)-mean(predicted(idxTgtOnsets,ivar)));
    corr_tgt(ivar,1) = corr(PSTH_f(idxTgtOnsets),predicted(idxTgtOnsets,ivar));
end