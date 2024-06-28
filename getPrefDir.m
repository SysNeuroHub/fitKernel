function [prefDir_quantized, prefDirTrials, tonsetRespAmp] = getPrefDir(y_r, t_r, onsetTimes, tgtDir, param)
% [prefDir, tonsetRespAmp] = getPrefDir(singleOnsetResp, tgtDir, winSamps, tOnRespWin, baseWin, cardinalDir)
% returns preferred direction in degree
%
% created from showTonsetResp.m

[avgOnsetResp, winSamps, singleOnsetResp, ...
    tgtDir_valid, uniqueOnsetLabels] ...
    = eventLockedAvg(y_r', t_r, onsetTimes, tgtDir, param.figTWin);

tonsetRespAmp = characteriseResp(singleOnsetResp, ...
    winSamps, param.tOnRespWin, param.baseWin, 'mean');
%trials x kinds


%prefDir = 180/pi*circ_mean(tgtDir'/180*pi,tonsetRespAmp_p(:,psthIdx)); %NG
%psthIdx = find(strcmp(psthNames, 'psth'));
[fitPars, fitErr, r] = fitResponse(tgtDir_valid, tonsetRespAmp, param.cardinalDir);
prefDir = fitPars(1); %[deg]
[~, prefDirIdx] = min(abs(circ_dist(pi/180*param.cardinalDir, pi/180*prefDir)));
prefDir_quantized = param.cardinalDir(prefDirIdx);
prefDirTrials =  find(tgtDir_valid == prefDir_quantized);%trials with cell's preferred direction

