function [Rsqadjusted_tgt] = getRsqadj_tgt(PSTH_f, predicted, catEvTimes, t_r, tWin, nPredictors)
% created from getExpVal_tgt(PSTH_f, predicted, catEvTimes, t_r, tWin)
% NOTE this function use pre-computed predicted traces

%% explained variance for target response
idxTgtOnsets = periTgtTidx(t_r, catEvTimes, tWin);

Rsqadjusted_tgt = zeros(size(predicted,2),1);
for ivar = 1:size(predicted,2)
    [Rsqadjusted_tgt(ivar,1)] = getRsqadj(PSTH_f(idxTgtOnsets), ...
        predicted(idxTgtOnsets,1), nPredictors);
end