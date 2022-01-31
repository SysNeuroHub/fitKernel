function [startSacc_after, endSacc_after] = selectSaccades(startSacc,endSacc, ...
    t_cat, excTimes, minSaccInterval)
% [startSacc_after, endSacc_after] = selectSaccades(startSacc,endSacc, ...
%     t_cat, excTimes, minSaccInterval)

% TODO: should exlude saccades which includes but  not end within excTimes


if nargin < 5
    minSaccInterval = [];
end

assert(length(startSacc)==length(endSacc));

%% omit successive saccades (omit initial leave last)
%minSaccInterval = 0.5;%[s]
if ~isempty(minSaccInterval)
    okSacc2 = find(diff(startSacc) > minSaccInterval);
    startSacc = startSacc(okSacc2);
    endSacc = endSacc(okSacc2);
end

%% adjust time axis of saccade times
[~, startSaccTidx] = arrayfun(@(x)(min(abs(t_cat - x))), startSacc);
[~, endSaccTidx] = arrayfun(@(x)(min(abs(t_cat - x))), endSacc);

%% omit saccades which starts within excTimes
[~,inclSaccIdx] = setdiff(startSaccTidx, find(excTimes));
startSaccTidx = startSaccTidx(inclSaccIdx);
endSaccTidx = endSaccTidx(inclSaccIdx);

%% omit saccades which ends within excTimes
%FIXME: should exlude saccades which includes but  not end within excTimes
[~,inclSaccIdx] = setdiff(endSaccTidx, find(excTimes));
startSaccTidx = startSaccTidx(inclSaccIdx);
endSaccTidx = endSaccTidx(inclSaccIdx);

startSacc_after = t_cat(startSaccTidx);
endSacc_after = t_cat(endSaccTidx);

