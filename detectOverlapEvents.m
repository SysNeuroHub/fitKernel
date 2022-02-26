function evOverlappedIdx = detectOverlapEvents(taxis, eventTimes, excEventTimes)
%evOverlappedIdx = detectOverlapEvents(taxis, eventTimes, excEventTimes)
%eventTimes, excEventTimes = [events x 2]

excTrace = event2Trace(taxis, excEventTimes);

evOverlapped = zeros(length(eventTimes),1);
for iev = 1:length(eventTimes)
    thisEvTrace = event2Trace(taxis, [eventTimes(iev,1) eventTimes(iev,end)]);
    evOverlapped(iev) = (thisEvTrace'*excTrace > 0);
end
evOverlappedIdx = find(evOverlapped);


%
% test = arrayfun(@(x)(event2Trace(t,x)), [drStartTimes_th drEndTimes_th], ...
%     'UniformOutput', false);
%cumsum(test{:,2}) - cumsum(test{:,1})
%drEndTraces=cellfun(@(x)cumsum(x), test(:,2), 'UniformOutput',
%false);outofmemory

%test=num2cell( [drStartTimes_th drEndTimes_th], 2);
%excMat = cellfun(@(x)(event2Trace(t,x)), test, 'UniformOutput', false);
%excludeTraces = cell2mat(excMat);%too large
