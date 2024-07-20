function latency_bhv = getTgtBhvLatency(onsets_cat, dd, validEvents, option)
%latency_bhv = getTgtBhvLatency(onsets_cat, dd, validEvents. option)
%
% option = 0: use "srt" stored in dd
% option = 1: cOnset - tOnset

if nargin < 4
    option = 0;
end
if nargin < 3
    validEvents = 1:numel(onsets_cat.tOnset);
end

switch option
    case 0
        latency_bhv = onsets_cat.cOnset(validEvents) - onsets_cat.tOnset(validEvents);
    case 1
        latency_bhv = dd.srt(validEvents);
end

% latency_bhv=getTgtBhvLatency(onsets_cat, validEvents);
% srt = dd.srt(validEvents);
% plot(latency_bhv, srt,'.'); hold on;
% 
% discrepancy = find(~isnan(latency_bhv)&isnan(srt));
% plot(latency_bhv(discrepancy), .5,'r*');
% 
% discrepancy2 = find(isnan(latency_bhv)&~isnan(srt));
% plot(.5,srt(discrepancy2),'r*');
% 
% xlabel('cOnset - tOnset');
% ylabel('srt = sacccade.onset - tOnset');
% 
% squareplot;
% axis equal padded

