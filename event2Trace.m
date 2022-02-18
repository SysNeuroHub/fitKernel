function trace = event2Trace(taxis, eventTimes, marginDur)
% trace = event2Trace(taxis, eventTimes)
%converts times of event into a single trace of 0/1
%
% trace = event2Trace(taxis, eventTimes, marginDur)
% imposes margin before and after each event 
%
% eventTimes: [events x 1] or [events x 2]
% if latter, first and 2nd columns represent start and end of an event

% 14/2/22 made faster

if nargin < 3
    marginDur = 0;
end

assert(size(eventTimes,1) >= size(eventTimes,2));

if size(eventTimes,2)==2
    assert(isempty(find(eventTimes(:,2)-eventTimes(:,1)<=0)));
end
%% implementation 1: SLOW
% trace = zeros(length(taxis),1);
% for itr = 1:length(eventTimes)
%     evStartT = max(taxis(1), eventTimes(itr,1) - 0.5*marginDur);
%     
%     evEndT = min(taxis(end), eventTimes(itr,end) + 0.5*marginDur);
% 
%     theseTimes = intersect(find(taxis>=evStartT), find(taxis<=evEndT));
%     
%     trace(theseTimes) = 1;
% end

%% implementation 2
%[~, evStartTidx] = arrayfun(@(x)(min(abs(taxis - x))), eventTimes(:,1) - 0.5*marginDur);
%evStartTidx(evStartTidx<1) = 1;
evStartTidx=interp1(taxis, 1:length(taxis), eventTimes(:,1) - 0.5*marginDur, 'nearest');%this is it!

% [~, evEndTidx] = arrayfun(@(x)(min(abs(taxis - x))), eventTimes(:,end) + 0.5*marginDur);
% evEndTidx(evEndTidx>length(taxis)) = length(taxis);
evEndTidx=interp1(taxis, 1:length(taxis), eventTimes(:,end) + 0.5*marginDur, 'nearest');%this is it!

[sameTime] = find(evStartTidx == evEndTidx);
evEndTidx(sameTime) = evEndTidx(sameTime)+1;
evEndTidx(evEndTidx>length(taxis)) = length(taxis);

onTrace = zeros(length(taxis),1);
onTrace(evStartTidx)=1;
offTrace = zeros(length(taxis),1);
offTrace(evEndTidx)=1;
trace = cumsum(onTrace) - cumsum(offTrace);
trace = (trace>0);

