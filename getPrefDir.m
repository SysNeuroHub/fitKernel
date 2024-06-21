function [prefDir_quantized, prefDirTrials, tonsetRespAmp] = getPrefDir(y_r, t_r, onsetTimes, tgtDir, param)
% [prefDir, tonsetRespAmp] = getPrefDir(singleOnsetResp, tgtDir, winSamps, tOnRespWin, baseWin, cardinalDir)
% returns preferred direction in degree
%
% created from showTonsetResp.m

[avgOnsetResp, winSamps, singleOnsetResp, ...
    sortedOnsetLabels, uniqueOnsetLabels] ...
    = eventLockedAvg(y_r', t_r, onsetTimes, tgtDir, param.figTWin);

tonsetRespAmp = characteriseResp(singleOnsetResp, ...
    winSamps, param.tOnRespWin, param.baseWin, 'mean');
%trials x kinds


%prefDir = 180/pi*circ_mean(tgtDir'/180*pi,tonsetRespAmp_p(:,psthIdx)); %NG
%psthIdx = find(strcmp(psthNames, 'psth'));
[fitPars, fitErr, r] = fitResponse(tgtDir, tonsetRespAmp, param.cardinalDir);
prefDir = fitPars(1); %[deg]
[~, prefDirIdx] = min(abs(circ_dist(pi/180*param.cardinalDir, pi/180*prefDir)));
prefDir_quantized = param.cardinalDir(prefDirIdx);
prefDirTrials =  find(tgtDir == prefDir_quantized);%trials with cell's preferred direction



% % created from dirTuinig_pupil.m
% % UNDER CONSTRUCTION
% 
% %% direction tuning
% dirs = unique(dd.targetloc);
% mresp = zeros(length(dirs),3);
% mparea = [];
% for idir = 1:length(dirs)
% 
% 
%     theseTr = intersect(find(dd.successTrials), find(dd.targetloc==dirs(idir)));
% 
%     resp_c = [];
%     parea_c = [];
%     for itr = 1:length(theseTr)
%         thisTr = theseTr(itr);
%         theseTimes  = intersect(find(eyeData(thisTr).t - dd.tOnset(thisTr) > 1e-3*respWin(1)), ...
%             find(eyeData(thisTr).t - dd.tOnset(thisTr) < 1e-3*respWin(2)));
%         resp_c(itr) = sum(psth_tr{thisTr}(theseTimes))/diff(respWin)/1e-3; %[spikes/s]
% 
%         theseTimes_pre  = intersect(find(eyeData(thisTr).t - dd.tOnset(thisTr) > 1e-3*preWin(1)), ...
%             find(eyeData(thisTr).t - dd.tOnset(thisTr) < 1e-3*preWin(2)));
% 
%        parea_c(itr) = mean(eyeData_rmblk_tr(thisTr).parea(theseTimes_pre));
%     end
%     mresp(idir,1) = mean(resp_c);
%     mresp(idir,2) = nanmean(resp_c(parea_c<parea_avg));
%     mresp(idir,3) = nanmean(resp_c(parea_c>parea_avg));
%     mparea(idir) = mean(parea_c);
%     length(find(parea_c>parea_avg))
% end
% 
% plot(dirs, mresp);
% legend('all trials','parea<mean','parea>mean');
% xlabel('tgt direction [deg]');
% %polarplot(dirs/180*pi, mresp);
