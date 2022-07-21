function [fitPars, fitErr, fittedResp] = fitResponse(deg, singleResp, cardinalDir, initPars)
%[fitPars, fitErr, fittedResp] = fitResponse(singleResp, sortedLabels, cardinalDir, initPars)
%
% INPUTS:
% deg: stimulus angle in degree
%
% OUTPUTS:
% [Dp(peak angle), Rp(peak amp), Rn(2nd peak amp = 0), Ro(baseline), sigma(width)]
%
%
%created from resp_cueConditions
%fitoriWrapped from Matteobox

if nargin < 4
    initPars = [];
end
avgResp = [];
for idir = 1:numel(cardinalDir)
    avgResp(idir) = mean(singleResp(deg == cardinalDir(idir)));
end
minResp = min(avgResp);
[fitPars, fitErr] ...
    = fitoriWrapped(deg, singleResp,...
    [], [nan nan 0 minResp nan],'',20, initPars);

if nargout > 2
    fittedResp ...
        = orituneWrapped(fitPars, cardinalDir);
    
    %plot(cardinalDir,  avgResp);hold on
    %plot(cardinalDir, fittedResp)
end

