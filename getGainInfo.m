function gainInfo = getGainInfo(t_r, y_r, cardinalDir, catEvTimes, dd, figTWin, onlySuccess, respWin)

cueDir = 0; %deg
[gain_dir, gain_distCue, winSamps_tonsetByCue, avgTonsetByCue] = getGainsByDist(t_r, ...
    y_r(:,1:2), cardinalDir, catEvTimes, dd, figTWin, onlySuccess, cueDir);

%avgTonsetByCue: direction x observed/mdl x time x w/wo cue
%prefDir = getPrefDir();
[~,respWinIdx(1)] = min(abs(winSamps_tonsetByCue - respWin(1)));
[~,respWinIdx(2)] = min(abs(winSamps_tonsetByCue - respWin(2)));

%% obtain quantized preferred direction from time-avg response of wo cue 
% resp = mean(avgTonsetByCue(:,1,respWinIdx(1):respWinIdx(2),2),3); 
resp = mean(avgTonsetByCue(:,1,respWinIdx(1):respWinIdx(2),1),3); 
[~,prefDirIdx] = max(resp);
prefDir = cardinalDir(prefDirIdx);

[~, gain_distPref] = getGainsByDist(t_r, ...
    y_r(:,1:2), cardinalDir, catEvTimes, dd, figTWin, onlySuccess, prefDir);

gainInfo.avgTonsetByCue = avgTonsetByCue;
gainInfo.prefDir = prefDir;
gainInfo.respWin = respWin;
gainInfo.winSamps = winSamps_tonsetByCue;
gainInfo.gain_dir = gain_dir;
gainInfo.gain_distCue = gain_distCue;
gainInfo.gain_distPref = gain_distPref;
gainInfo.cardinalDir = cardinalDir;
gainInfo.cardinalDist = [0 45 90 135 180];
end
