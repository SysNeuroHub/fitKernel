function trace = event2Trace(taxis, eventTimes, marginDur)
% trace = event2Trace(taxis, eventTimes)
%converts times of event into a single trace of 0/1
%
if nargin < 3
    marginDur = 0;
end

trace = zeros(length(taxis),1);
for itr = 1:length(eventTimes)
    excStartT = max(taxis(1), eventTimes(itr) - 0.5*marginDur);
    
    excEndT = min(taxis(end), eventTimes(itr) + 0.5*marginDur);
    theseTimes = intersect(find(taxis>=excStartT), find(taxis<=excEndT));
    
    trace(theseTimes) = 1;
end