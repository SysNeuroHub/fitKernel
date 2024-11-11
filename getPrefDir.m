function [prefDir_quantized, prefDirTrials, tonsetRespAmp] = getPrefDir(y_r, t_r, onsetTimes, tgtDir, param)
% [prefDir, tonsetRespAmp] = getPrefDir(singleOnsetResp, tgtDir, winSamps, tOnRespWin, baseWin, cardinalDir)
% returns preferred direction in degree
%
% created from showTonsetResp.m
% TODO: deal with inhibition

[avgOnsetResp, winSamps, singleOnsetResp, ...
    tgtDir_valid, uniqueOnsetLabels] ...
    = eventLockedAvg(y_r', t_r, onsetTimes, tgtDir, param.figTWin);

tonsetRespAmp = characteriseResp(singleOnsetResp, ...
    winSamps, param.tOnRespWin, param.baseWin, 'mean');
%trials x kinds

% %% decide response polarity
% bWin =  intersect(find(winSamps_t >= param.baseWin(1)), find(winSamps_t <= param.baseWin(2)));
% rWin =  intersect(find(winSamps_t >= param.respWin(1)), find(winSamps_t <= param.respWin(2)));
% bSignal = mean(singleResp_t(tgtDir == 0, bWin),2);
% rSignal = mean(singleResp_t(tgtDir == 0, rWin),2);
% %p_visResp = signrank(bSignal, rSignal); %is there a modulation after target stimulation?
% respPolarity = sign(mean(rSignal) - mean(bSignal));


%prefDir = 180/pi*circ_mean(tgtDir'/180*pi,tonsetRespAmp_p(:,psthIdx)); %NG
%psthIdx = find(strcmp(psthNames, 'psth'));
[fitPars, fitErr, r] = fitResponse(tgtDir_valid, tonsetRespAmp, param.cardinalDir);
prefDir = fitPars(1); %[deg]
[~, prefDirIdx] = min(abs(circ_dist(pi/180*param.cardinalDir, pi/180*prefDir)));
prefDir_quantized = param.cardinalDir(prefDirIdx);
prefDirTrials =  find(tgtDir_valid == prefDir_quantized);%trials with cell's preferred direction

