function [evTimes, evStart, evEnd] = trace2Event(trace, taxis)
% evTimes = trace2Event(trace)
% rerutns index of on and off times in trace
%
% evTimes = trace2Event(trace, taxis)
% rerutns times of on and off 
%
% see also: event2Trace

trace = logical(trace);

if size(trace,1)<size(trace,2)
    trace = trace';
end

diffTrace = [0; diff(trace)];

evStartIdx = find(diffTrace>0);
if trace(1) == 1
    evStartIdx = [1; evStartIdx];
end

evEndIdx = find(diffTrace<0);
if trace(end)==1
    evEndIdx = [evEndIdx; length(trace)];
end

assert(length(evStartIdx) == length(evEndIdx));

if nargin==2
    evStart = taxis(evStartIdx);
    evEnd = taxis(evEndIdx);
else
    evStart = evStartIdx;
    evEnd = evEndIdx;
end

evTimes = [evStart evEnd];

%% sanity check
% [evStartIdx, evEndIdx] = trace2Event(trace);
% trace_restored = event2Trace(1:length(trace),[evStartIdx evEndIdx]);
% plot(trace);hold on;
% plot(trace_restored)
% plot(evStartIdx,.5,'ro',evEndIdx,.5,'gx');


