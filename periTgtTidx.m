function idxTgtOnsets = periTgtTidx(t_r, catEvTimes, tWin)
% idxTgtOnsets = periTgtTidx(t_r, catEvTimes, tWin)
% created from getExpVal_tgt.m

idxTgtOnsets = [];
onsetTimes = catEvTimes.tOnset;
onsetTimes = onsetTimes(~isnan(onsetTimes));
for ievent = 1:numel(onsetTimes)
    idxTgtOnsets = cat(1, idxTgtOnsets, ...
        find((t_r>=onsetTimes(ievent)+tWin(1)) .* (t_r<onsetTimes(ievent)+tWin(2))));
end