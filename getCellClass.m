function unitClass = getCellClass(PtonsetResp, PtonsetResp_paired, mtOnsetRespAmp, PsaccResp, Pth)
% unitClass = getCellClass(PtonsetResp, PtonsetResp_paired, mtOnsetRespAmp)
% returns a functional cell class, defined as follows:
%1: significant response to tOnset, driven by vision
%2: significant response to tOnset, driven by eye velocity
%3: significant response to tOnset, driven by both vision & eye velocity
%4: no response to tOnset, but significant response to saccade outside of the task
%5: no response to tOnset or saccades outside of the task
%
% INPUTS:
% mtOnsetRespAmp(1): response of vision kernel
% mtOnsetRespAmp(2): response of eye velocity kernel

if nargin < 5
    Pth = 0.05;
end

visionIdx = 1;
eyeVelIdx = 2;

if PtonsetResp < Pth
    if PtonsetResp_paired > Pth
        unitClass = 3;
    elseif mtOnsetRespAmp(visionIdx) < mtOnsetRespAmp(eyeVelIdx)
        unitClass = 2;
    elseif mtOnsetRespAmp(visionIdx) > mtOnsetRespAmp(eyeVelIdx)
        unitClass = 1;
    end
elseif PsaccResp < Pth
    unitClass = 4;
else
    unitClass = 5;
end
end
