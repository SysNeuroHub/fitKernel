function [f, cellclassInfo] = showTonsetResp(t_r, y_r, catEvTimes, dd, psthNames, ...
    startSaccNoTask, saccDirNoTask, param, figTWin)
%[f, cellclassInfo] = showTonsetResp(t_r, y_r, catEvTimes, dd, psthNames, ...
%    startSaccNoTask, saccDirNoTask, param, figTWin)
%
%created from classifyUnits_circ & resp_cueConditions.m


allTr = 0;

psthIdx = find(strcmp(psthNames, 'psth'));
allMdlIdx = find(strcmp(psthNames, 'predicted_all'));
visionIdx = find(strcmp(psthNames, 'vision'));
eyevelIdx = find(strcmp(psthNames, 'eyespeed'));

%% triggered by tOnsets
%validEvents = intersect(find(~isnan(onset)), find(dd.cueOn==icue-1));
validEvents = intersect(find(~isnan(catEvTimes.tOnset)), find(dd.successTrials));
%only use trials when the choices were registered.
%this is a temporary fix as my current algorithm assumes stimuli were NOT
%presented, causing no visual response in the model

onsetTimes = catEvTimes.tOnset(validEvents);
tgtDir = getTgtDir(dd.targetloc(validEvents), param.cardinalDir);

[~,dirIdx] = intersect(param.cardinalDir, unique(tgtDir));

[avgOnsetResp, winSamps, singleOnsetResp, ...
    sortedOnsetLabels, uniqueOnsetLabels] ...
    = eventLockedAvg(y_r', t_r, onsetTimes, tgtDir, param.figTWin);


% respTidx = intersect(find(param.tOnRespWin(1)<=winSamps), ...
%     find(param.tOnRespWin(2)>=winSamps));


tonsetRespAmp = characteriseResp(singleOnsetResp, ...
    winSamps, param.tOnRespWin, param.baseWin, 'mean');
%trials x kinds
tonsetRespAmp_p = characteriseResp(singleOnsetResp, ...
    winSamps, param.tOnRespWin, [], 'mean');
tonsetRespAmp_b = characteriseResp(singleOnsetResp, ...
    winSamps, param.baseWin, [], 'mean');

%prefDir = 180/pi*circ_mean(tgtDir'/180*pi,tonsetRespAmp_p(:,psthIdx)); %NG
[fitPars, fitErr, r] = fitResponse(tgtDir, tonsetRespAmp(:,psthIdx), param.cardinalDir);
prefDir = fitPars(1); %[deg]
[~, prefDirIdx] = min(abs(circ_dist(pi/180*param.cardinalDir, pi/180*prefDir)));

%% check if a unit is visually responsive
if allTr
    theseTrials = 1:numel(tgtDir);
else
    theseTrials = find(tgtDir == param.cardinalDir(prefDirIdx));%trials with cell's preferred direction
%     for tt = 1:numel(tgtDir)
%         [~,minDirIdx(tt)] = min(abs(circ_dist(tgtDir(tt), pi/180*param.cardinalDir))); %from getTgtDirMtx 23/7/22
%     end
%     theseTrials = find(minDirIdx == prefDirIdx);
end

%plot(tgtDir'/180*pi, tonsetRespAmp_p(:,psthIdx),'.');
%hold on
%plot(tgtDir'/180*pi, tonsetRespAmp_b(:,psthIdx),'r.');


%         PtonsetResp = [];
%         for tt = 1:numel(respTidx)
%             PtonsetResp(tt) = dirDotProdTest(tgtDir'/180*pi, ...
%                 singleOnsetResp(:,psthIdx,respTidx(tt)), tonsetRespAmp_b(:,psthIdx));
%         end
PtonsetResp = dirDotProdTest(tgtDir'/180*pi, ...
    tonsetRespAmp_p(:,psthIdx), tonsetRespAmp_b(:,psthIdx));

%TODO: cluster-base permutation test to determine the P-value
%d-prime rather than the dot product test?

mtOnsetResp = squeeze(mean(singleOnsetResp(theseTrials,:,:)));%avg response to preferred direction
setOnsetResp = 1/sqrt(numel(theseTrials))*squeeze(std(singleOnsetResp(theseTrials,:,:)));%avg response to preferred direction

%% mean response amplitude across trials
%mtOnsetRespAmp_p_c = abs(mean(tonsetRespAmp_p .* exp(1i*tgtDir'/180*pi)));
%mtOnsetRespAmp_b_c = abs(mean(tonsetRespAmp_b .* exp(1i*tgtDir'/180*pi)));
mtOnsetRespAmp = mean(tonsetRespAmp(theseTrials,[psthIdx visionIdx eyevelIdx]),1);
%mtOnsetRespAmp = mtOnsetRespAmp_p_c - mtOnsetRespAmp_b_c;

if ~isempty(theseTrials)
    [PtonsetResp_paired] = signrank(tonsetRespAmp(theseTrials,visionIdx), ...
        tonsetRespAmp(theseTrials,eyevelIdx));
else 
    PtonsetResp_paired = nan;
end


%% triggered by saccade onsets (outside of the task)
[avgSaccResp, ~, singleSaccResp, sortedSaccLabels, uniqueSaccLabels] ...
    = eventLockedAvg(y_r',t_r, startSaccNoTask, saccDirNoTask, figTWin);

%use the same saccade direction to the one used for tOnset
if allTr
    theseSaccTrials = 1:numel(saccDirNoTask);
else
    theseSaccTrials = (saccDirNoTask == param.cardinalDir(prefDirIdx));
end
msaccResp = squeeze(nanmean(singleSaccResp(theseSaccTrials,:,:)));%avg response to preferred direction
sesaccResp = 1/sqrt(numel(theseSaccTrials))*squeeze(nanstd(singleSaccResp(theseSaccTrials,:,:),1));

saccRespAmp = characteriseResp(singleSaccResp, ...
    winSamps, param.tOnRespWin, param.baseWin, 'mean');
saccRespAmp_p = characteriseResp(singleSaccResp, ...
    winSamps, param.tOnRespWin, [], 'mean');
saccRespAmp_b = characteriseResp(singleSaccResp, ...
    winSamps, param.baseWin, [], 'mean');
msaccRespAmp = mean(saccRespAmp(theseSaccTrials,psthIdx),1);

nonanEvents = intersect(find(~isnan(saccRespAmp_p(:,1))), find(~isnan(saccRespAmp_b(:,1))));
if isempty(nonanEvents)
    return;
end
PsaccResp = dirDotProdTest(saccDirNoTask(nonanEvents)'/180*pi, ...
    saccRespAmp_p(nonanEvents,psthIdx), saccRespAmp_b(nonanEvents,psthIdx));


%% triggered by tOnset without saccade
validEvents = intersect(find(~isnan(catEvTimes.tOnset)), find(dd.successTrials==0));
%if ~isempty(validEvents)
    onsetTimes = catEvTimes.tOnset(validEvents);
    tgtDir_v = getTgtDir(dd.targetloc(validEvents), param.cardinalDir);
    
    [avgOnsetResp_v, winSamps, singleOnsetResp_v, ...
        sortedOnsetLabels_v, uniqueOnsetLabels_v] ...
        = eventLockedAvg(y_r', t_r, onsetTimes, tgtDir_v, param.figTWin);
    
    if allTr
        theseTrials_v = 1:numel(tgtDir_v);
    else
        theseTrials_v = find(tgtDir_v == param.cardinalDir(prefDirIdx));%trials with cell's preferred direction
    end
    mtOnsetResp_v = squeeze(mean(singleOnsetResp_v(theseTrials_v,:,:)));%avg response to preferred direction
    setOnsetResp_v = 1/sqrt(numel(theseTrials_v))*squeeze(std(singleOnsetResp_v(theseTrials_v,:,:)));%avg response to preferred direction
%end

%% categorise the cell
%unitClass = getCellClass2(PtonsetResp, PsaccResp, param.Pth);
unitClass = getCellClass(PtonsetResp, PtonsetResp_paired, mtOnsetRespAmp, PsaccResp, param.Pth);


%% save results
cellclassInfo.unitClass = unitClass;
cellclassInfo.PtonsetResp = PtonsetResp;
%cellclassInfo.PtonsetResp_paired = PtonsetResp_paired;
cellclassInfo.mtOnsetRespAmp = mtOnsetRespAmp;
cellclassInfo.mtOnsetResp = mtOnsetResp;
cellclassInfo.npreftonsetTrials = numel(theseTrials);
cellclassInfo.PsaccResp = PsaccResp;
cellclassInfo.msaccResp = msaccResp;
cellclassInfo.msaccRespAmp = msaccRespAmp;
cellclassInfo.nprefSaccTrials = numel(theseSaccTrials);
cellclassInfo.winSamps = winSamps;
%cellclassInfo.datech = datech;

%save(saveName, 'cellclassInfo');

medianSaccDelay = nanmedian(catEvTimes.cOnset-catEvTimes.tOnset);


%% visualize the result
f=figure('position',[0 0 800 1000]);
ax(1) = subplot(311);
%plot(winSamps, mtOnsetResp([allMdlIdx visionIdx eyevelIdx],:));hold on;
boundedline(winSamps, mtOnsetResp(psthIdx,:), setOnsetResp(psthIdx,:),'k', 'linewidth',2);
hold on;
boundedline(winSamps, mtOnsetResp(allMdlIdx,:), setOnsetResp(allMdlIdx,:),'b', 'transparency', 0.5);
boundedline(winSamps, mtOnsetResp(visionIdx,:), setOnsetResp(visionIdx,:),'m', 'transparency', 0.5);
boundedline(winSamps, mtOnsetResp(eyevelIdx,:), setOnsetResp(eyevelIdx,:),'c', 'transparency', 0.5);
%boundedline(winSamps, mtOnsetResp(5,:), setOnsetResp(5,:),'g', 'transparency', 0.5);
%vbox(param.baseWin(1), param.baseWin(2))
%vbox(param.tOnRespWin(1), param.tOnRespWin(2),[],[.7 1 .7]);
xlim([-0.1 0.5]);
ylabel(['tOnset (success), n=' num2str(numel(theseTrials))]);
%title([datech ' ' num2str(unitClass)]);
title(['resp to' num2str(param.cardinalDir(prefDirIdx)) 'deg']);

ax(2) = subplot(312);
boundedline(winSamps, mtOnsetResp_v(psthIdx,:), setOnsetResp_v(psthIdx,:),'k', 'linewidth',2);
hold on;
boundedline(winSamps, mtOnsetResp_v(allMdlIdx,:), setOnsetResp_v(allMdlIdx,:),'b', 'transparency', 0.5);
boundedline(winSamps, mtOnsetResp_v(visionIdx,:), setOnsetResp_v(visionIdx,:),'m', 'transparency', 0.5);
boundedline(winSamps, mtOnsetResp_v(eyevelIdx,:), setOnsetResp_v(eyevelIdx,:),'c', 'transparency', 0.5);
%vbox(param.baseWin(1), param.baseWin(2))
%vbox(param.tOnRespWin(1), param.tOnRespWin(2),[],[.7 1 .7]);
xlim([-0.1 0.5]);
ylabel(['tOnset (fail), n=' num2str(numel(theseTrials_v))] );

ax(3) = subplot(313);
boundedline(winSamps, msaccResp(psthIdx,:), sesaccResp(psthIdx,:),'k', 'linewidth',2);
hold on;
boundedline(winSamps, msaccResp(allMdlIdx,:), sesaccResp(allMdlIdx,:),'b', 'transparency', 0.5);
boundedline(winSamps, msaccResp(visionIdx,:), sesaccResp(visionIdx,:),'m', 'transparency', 0.5);
boundedline(winSamps, msaccResp(eyevelIdx,:), sesaccResp(eyevelIdx,:),'c', 'transparency', 0.5);

linkaxes(ax(:),'y');
vline(medianSaccDelay, ax(1));
vline(medianSaccDelay, ax(2));
vline(0, ax(3));

%vbox(param.baseWin(1), param.baseWin(2))
%vbox(param.tOnRespWin(1), param.tOnRespWin(2),[],[.7 1 .7]);
xlim([-0.1-medianSaccDelay winSamps(end)-medianSaccDelay])
ylabel(['saccade(outside task), n=' num2str(numel(theseSaccTrials))]);
legend('observed','all mdl','vision','eye velocity','location','northwest')


% saveas(gcf,saveFigName);
% screen2png(saveFigName(1:end-4));
% close;


end

