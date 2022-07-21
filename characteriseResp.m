function [amplitude] = characteriseResp(singleTrialResp, winSamps, respWin, baseWin, type)
%[amp] = characteriseResp(singleTrialResp, winSamps, respWin, baseWin)
% returns response amplitude at each trial basis
%
%INPUTS:
% singleTrialResp: events x kinds x times
% winSamps: 1xtimes
%
%OUTPUTS:
% amp: events x kinds
% latency: not yet implemented

if nargin < 5
    type = 'abspeak';
end
assert(size(singleTrialResp,3) == length(winSamps));

respIdx = intersect(find(respWin(1)<=winSamps), find(respWin(2)>=winSamps));

respSignal = (singleTrialResp(:,:,respIdx));
if isempty(baseWin)
    baseSignal = 0;
else
    baseIdx = intersect(find(baseWin(1)<=winSamps), find(baseWin(2)>=winSamps));
    baseSignal = (mean(singleTrialResp(:,:,baseIdx),3));
end

switch type
    case 'abspeak'
        amplitude = max(abs(respSignal - baseSignal),[],3);
    case 'mean'
        amplitude = mean(respSignal,3) - baseSignal;
end