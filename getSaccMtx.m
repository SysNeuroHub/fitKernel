function saccade = getSaccMtx(t_r, startSacc, endSacc, eyeData_cat, cardinalDir)
%returns [direction x time]

[~, dirIndex] = getSaccDir(startSacc, endSacc, eyeData_cat, cardinalDir);

[~, startSaccTidx] = arrayfun(@(x)(min(abs(t_r - x))), startSacc);
[~, endSaccTidx] = arrayfun(@(x)(min(abs(t_r - x))), endSacc);

%if start time = end time, delay the end time one time bin
[sameTime] = find(endSaccTidx == startSaccTidx);
endSaccTidx(sameTime) = endSaccTidx(sameTime)+1;
endSaccTidx(endSaccTidx>length(t_r)) = length(t_r);


saccade = zeros(length(cardinalDir),length(t_r));
% for isacc = 1:length(dirIndex)
%     saccade(dirIndex(isacc),startSaccTidx(isacc):endSaccTidx(isacc)) = 1;
% end

for idir = 1:length(cardinalDir)
    theseSaccs = find(dirIndex == idir);
    onTrace = zeros(length(t_r),1);
    onTrace(startSaccTidx(theseSaccs))=1;
    offTrace = zeros(length(t_r),1);
    offTrace(endSaccTidx(theseSaccs))=1;
    thisTrace = cumsum(onTrace) - cumsum(offTrace);
    saccade(idir,:) = thisTrace;
end

