function [startSacc_after, endSacc_after] = selectSaccades(startSacc,endSacc, ...
    t_cat, excTrace, minSaccInterval, marginDur)
% [startSacc_after, endSacc_after] = selectSaccades(startSacc,endSacc, ...
%     t_cat, excTimes, minSaccInterval)

% TODO: should exlude saccades which includes but  not end within excTimes

if nargin < 6
    marginDur = 0;
end

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
%startSaccTidx=interp1(t_cat, 1:length(t_cat), startSacc, 'nearest');
%endSaccTidx=interp1(t_cat, 1:length(t_cat), endSacc, 'nearest');


%% omit saccades which starts within excTimes
%[~,inclSaccIdx] = setdiff(startSaccTidx, find(excTimes));
%startSaccTidx = startSaccTidx(inclSaccIdx);
%endSaccTidx = endSaccTidx(inclSaccIdx);

%% omit saccades which ends within excTimes
%FIXME: should exlude saccades which includes but  not end within excTimes
%[~,inclSaccIdx] = setdiff(endSaccTidx, find(excTimes));
%startSaccTidx = startSaccTidx(inclSaccIdx);
%endSaccTidx = endSaccTidx(inclSaccIdx);

eventTimes = [startSacc-0.5*marginDur endSacc+0.5*marginDur];
eventTimes(eventTimes<min(t_cat)) = min(t_cat);
eventTimes(eventTimes>max(t_cat)) = max(t_cat);

excEventTimes = trace2Event(excTrace, t_cat);
evOverlappedIdx = detectOverlapEvents(t_cat, eventTimes, excEventTimes);
OkIdx = setdiff(1:length(startSacc), evOverlappedIdx);

startSacc_after = startSacc(OkIdx);
endSacc_after = endSacc(OkIdx);

