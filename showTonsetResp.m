function [f, cellclassInfo] = showTonsetResp(t_r, y_r, catEvTimes, dd, psthNames, ...
    startSaccNoTask, saccDirNoTask, param, allTr, includeTrials)
%[f, cellclassInfo] = showTonsetResp(t_r, y_r, catEvTimes, dd, psthNames, ...
%    startSaccNoTask, saccDirNoTask, param, figTWin)
%
% cellclassInfo:
%     cellclassInfo.unitClass: (PROVISIONAL)
%     cellclassInfo.PtonsetResp: p-value if the target response is deviated from (0,) in 2D space (dirDotProdTest.m)
%     cellclassInfo.mtOnsetRespAmp:
%     cellclassInfo.mtOnsetResp:
%     cellclassInfo.npreftonsetTrials:
%     cellclassInfo.PsaccResp: p-value if the saccade response is significant
%     cellclassInfo.msaccResp:
%     cellclassInfo.msaccRespAmp:
%     cellclassInfo.nprefSaccTrials:
%     cellclassInfo.winSamps:

if nargin < 10
    includeTrials = 1:dd.numTrials;
end
%created from classifyUnits_circ & resp_cueConditions.m
if nargin < 9 || isempty(allTr)
    allTr = 0; %only show response to preferred direction
end

figTWin = param.figTWin; %figure temporal window
compTWin = [figTWin(1) - 0.2 figTWin(2)]; %triggered response is computed

psthIdx = find(strcmp(psthNames, 'psth'));
allMdlIdx = find(strcmp(psthNames, 'predicted_all'));
visionIdx = find(strcmp(psthNames, 'vision'));
eyevelIdx = find(strcmp(psthNames, 'eyespeed'));
eyeposIdx = find(strcmp(psthNames, 'eyeposition'));

%hesitant = getChoiceOutcome_hesitant(catEvTimes, dd);
[choiceOutcome] = getChoiceOutcome(dd);


%% triggered by tOnsets
successEvents = intersect(find((choiceOutcome==1)), includeTrials);
successEvents = intersect(successEvents, find(catEvTimes.tOnset + compTWin(2) < max(t_r)));
%only use trials when the choices were registered.
%this is a temporary fix as my current algorithm assumes stimuli were NOT
%presented, causing no visual response in the model

onsetTimes = catEvTimes.tOnset(successEvents);
tgtDir = getTgtDir(dd.targetloc(successEvents), param.cardinalDir);

[~,dirIdx] = intersect(param.cardinalDir, unique(tgtDir));

[~, winSamps, singleOnsetResp, ...
    sortedOnsetLabels, uniqueOnsetLabels] ...
    = eventLockedAvg(y_r', t_r, onsetTimes, tgtDir, compTWin);
tonsetRespAmp = characteriseResp(singleOnsetResp, ...
    winSamps, param.tOnRespWin, param.baseWin, 'mean');

% prefDir = getPrefDir(y_r(:,1), t_r, onsetTimes, tgtDir, param);
prefDir = getPrefDir_wrapper(y_r, t_r, dd, catEvTimes, param, includeTrials);

tonsetRespAmp_p = characteriseResp(singleOnsetResp, ...
    winSamps, param.tOnRespWin, [], 'mean');
tonsetRespAmp_b = characteriseResp(singleOnsetResp, ...
    winSamps, param.baseWin, [], 'mean');


%% select trials
if allTr
    theseTrials = 1:numel(tgtDir);
else
%    theseTrials = prefDirTrials;%trials with cell's preferred direction
    theseTrials = find(tgtDir == prefDir);
end

%plot(tgtDir'/180*pi, tonsetRespAmp_p(:,psthIdx),'.');
%hold on
%plot(tgtDir'/180*pi, tonsetRespAmp_b(:,psthIdx),'r.');

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

%% triggered by cOnsets not tOnsets
onsetTimes = catEvTimes.cOnset(successEvents);
[avgOnsetResp, winSamps, singleOnsetResp, ...
    sortedOnsetLabels, uniqueOnsetLabels] ...
    = eventLockedAvg(y_r', t_r, onsetTimes, tgtDir, compTWin);
mcOnsetResp = squeeze(mean(singleOnsetResp(theseTrials,:,:)));%avg response to preferred direction
secOnsetResp = 1/sqrt(numel(theseTrials))*squeeze(std(singleOnsetResp(theseTrials,:,:)));%avg response to preferred direction


%% triggered by saccade onsets (outside of the task)

[avgSaccResp, ~, singleSaccResp, sortedSaccLabels, uniqueSaccLabels] ...
    = eventLockedAvg(y_r',t_r, startSaccNoTask, saccDirNoTask, compTWin);

%use the same saccade direction to the one used for tOnset
if allTr
    theseSaccTrials = 1:numel(sortedSaccLabels);
else
    theseSaccTrials = find(sortedSaccLabels == prefDir);
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
quiescentEvents = intersect(find(choiceOutcome == 3), includeTrials);
onsetTimes = catEvTimes.tOnset(quiescentEvents);
tgtDir_v = getTgtDir(dd.targetloc(quiescentEvents), param.cardinalDir);

[avgOnsetResp_v, winSamps, singleOnsetResp_v, ...
    sortedOnsetLabels_v, uniqueOnsetLabels_v] ...
    = eventLockedAvg(y_r', t_r, onsetTimes, tgtDir_v, compTWin);

if allTr
    theseTrials_v = 1:size(singleOnsetResp_v,1);
else
    theseTrials_v = find(sortedOnsetLabels_v == prefDir);%trials with cell's preferred direction
end
mtOnsetResp_v = squeeze(mean(singleOnsetResp_v(theseTrials_v,:,:)));%avg response to preferred direction
setOnsetResp_v = 1/sqrt(numel(theseTrials_v))*squeeze(std(singleOnsetResp_v(theseTrials_v,:,:)));%avg response to preferred direction


%% triggered by tOnset with wrong succade direction
wrongdirEvents = find(choiceOutcome == 2);
onsetTimes = catEvTimes.tOnset(wrongdirEvents);
tgtDir_f = getTgtDir(dd.targetloc(wrongdirEvents), param.cardinalDir);

[avgOnsetResp_f, winSamps, singleOnsetResp_f, ...
    sortedOnsetLabels_f, uniqueOnsetLabels_f] ...
    = eventLockedAvg(y_r', t_r, onsetTimes, tgtDir_f, compTWin);

if allTr
    theseTrials_f = 1:size(singleOnsetResp_f,1);
else
    theseTrials_f = find(sortedOnsetLabels_f == prefDir);%trials with cell's preferred direction
end
mtOnsetResp_f = squeeze(mean(singleOnsetResp_f(theseTrials_f,:,:)));%avg response to preferred direction
setOnsetResp_f = 1/sqrt(numel(theseTrials_f))*squeeze(std(singleOnsetResp_f(theseTrials_f,:,:)));%avg response to preferred direction


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

%% delays only use successful trials
medianSaccDelay = nanmedian(catEvTimes.cOnset(successEvents) - catEvTimes.tOnset(successEvents));
medianFixationDelay = nanmedian(catEvTimes.fOnset(successEvents) - catEvTimes.tOnset(successEvents));

%% visualize the result
f=figure('position',[0 0 800 1000]);
ax(1) = subplot(511);
%plot(winSamps, mtOnsetResp([allMdlIdx visionIdx eyevelIdx],:));hold on;
boundedline(winSamps, mtOnsetResp(psthIdx,:), setOnsetResp(psthIdx,:),'k', 'linewidth',2);
hold on;
boundedline(winSamps, mtOnsetResp(allMdlIdx,:), setOnsetResp(allMdlIdx,:),'b', 'transparency', 0.5);
boundedline(winSamps, mtOnsetResp(visionIdx,:), setOnsetResp(visionIdx,:),'m', 'transparency', 0.5);
if ~isempty(eyevelIdx)
    boundedline(winSamps, mtOnsetResp(eyevelIdx,:), setOnsetResp(eyevelIdx,:),'c', 'transparency', 0.5);
end
boundedline(winSamps, mtOnsetResp(eyeposIdx,:), setOnsetResp(eyeposIdx,:),'g', 'transparency', 0.5);
%boundedline(winSamps, mtOnsetResp(5,:), setOnsetResp(5,:),'g', 'transparency', 0.5);
%vbox(param.baseWin(1), param.baseWin(2))
%vbox(param.tOnRespWin(1), param.tOnRespWin(2),[],[.7 1 .7]);
set(gca,'tickdir','out');
ylabel(sprintf('tOnset (success)\n n=%d',numel(theseTrials)));
[~,ch]=fileparts(dd.path);
tname = [dd.subject ' ' dd.date ' ' ch];
if allTr
    title(['resp to all directions']);
else
    title(['resp to' num2str(prefDir) 'deg']);
end
axis tight;
xlim([-0.1 0.5]);

ax(2) = subplot(512);
boundedline(winSamps, mtOnsetResp_v(psthIdx,:), setOnsetResp_v(psthIdx,:),'k', 'linewidth',2);
hold on;
boundedline(winSamps, mtOnsetResp_v(allMdlIdx,:), setOnsetResp_v(allMdlIdx,:),'b', 'transparency', 0.5);
boundedline(winSamps, mtOnsetResp_v(visionIdx,:), setOnsetResp_v(visionIdx,:),'m', 'transparency', 0.5);
if ~isempty(eyevelIdx)
    boundedline(winSamps, mtOnsetResp_v(eyevelIdx,:), setOnsetResp_v(eyevelIdx,:),'c', 'transparency', 0.5);
end
boundedline(winSamps, mtOnsetResp_v(eyeposIdx,:), setOnsetResp_v(eyeposIdx,:),'g', 'transparency', 0.5);
%vbox(param.baseWin(1), param.baseWin(2))
%vbox(param.tOnRespWin(1), param.tOnRespWin(2),[],[.7 1 .7]);
axis tight;
ylabel(sprintf('tOnset (quiescent) \nn=%d', numel(theseTrials_v)));
set(gca,'tickdir','out');

ax(3) = subplot(513);
boundedline(winSamps, mtOnsetResp_f(psthIdx,:), setOnsetResp_f(psthIdx,:),'k', 'linewidth',2);
hold on;
boundedline(winSamps, mtOnsetResp_f(allMdlIdx,:), setOnsetResp_f(allMdlIdx,:),'b', 'transparency', 0.5);
boundedline(winSamps, mtOnsetResp_f(visionIdx,:), setOnsetResp_f(visionIdx,:),'m', 'transparency', 0.5);
if ~isempty(eyevelIdx)
    boundedline(winSamps, mtOnsetResp_f(eyevelIdx,:), setOnsetResp_f(eyevelIdx,:),'c', 'transparency', 0.5);
end
boundedline(winSamps, mtOnsetResp_f(eyeposIdx,:), setOnsetResp_f(eyeposIdx,:),'g', 'transparency', 0.5);
%vbox(param.baseWin(1), param.baseWin(2))
%vbox(param.tOnRespWin(1), param.tOnRespWin(2),[],[.7 1 .7]);
axis tight;
ylabel(sprintf('tOnset (fail) \nn=%d', numel(theseTrials_f)) );
set(gca,'tickdir','out');
linkaxes(ax(1:3),'x');
xlim(figTWin);

ax(4) = subplot(514);
boundedline(winSamps, msaccResp(psthIdx,:), sesaccResp(psthIdx,:),'k', 'linewidth',2);
hold on;
boundedline(winSamps, msaccResp(allMdlIdx,:), sesaccResp(allMdlIdx,:),'b', 'transparency', 0.5);
boundedline(winSamps, msaccResp(visionIdx,:), sesaccResp(visionIdx,:),'m', 'transparency', 0.5);
if ~isempty(eyevelIdx)
    boundedline(winSamps, msaccResp(eyevelIdx,:), sesaccResp(eyevelIdx,:),'c', 'transparency', 0.5);
end
boundedline(winSamps, msaccResp(eyeposIdx,:), sesaccResp(eyeposIdx,:),'g', 'transparency', 0.5);
axis tight;
xlim([figTWin(1)-medianSaccDelay figTWin(end)-medianSaccDelay])
ylabel(sprintf('saccade(outside task) \nn=%d', numel(theseSaccTrials)));

ax(5) = subplot(515);
boundedline(winSamps, mcOnsetResp(psthIdx,:), secOnsetResp(psthIdx,:),'k', 'linewidth',2);
hold on;
boundedline(winSamps, mcOnsetResp(allMdlIdx,:), secOnsetResp(allMdlIdx,:),'b', 'transparency', 0.5);
boundedline(winSamps, mcOnsetResp(visionIdx,:), secOnsetResp(visionIdx,:),'m', 'transparency', 0.5);
if ~isempty(eyevelIdx)
    boundedline(winSamps, mcOnsetResp(eyevelIdx,:), secOnsetResp(eyevelIdx,:),'c', 'transparency', 0.5);
end
boundedline(winSamps, mcOnsetResp(eyeposIdx,:), secOnsetResp(eyeposIdx,:),'g', 'transparency', 0.5);
axis tight;
xlim([figTWin(1)-medianSaccDelay figTWin(end)-medianSaccDelay])

linkaxes(ax(:),'y');
vline(medianSaccDelay, ax(1));
vline(medianSaccDelay, ax(2));
vline(medianSaccDelay, ax(3));
if medianFixationDelay>winSamps(1)
    vline(medianFixationDelay, ax(1));
    vline(medianFixationDelay, ax(2));
    vline(medianFixationDelay, ax(3));
end
vline(0, ax(4));
vline(0, ax(5));
set(gca,'tickdir','out');

ylabel(sprintf('cOnset (success)\nn=%d',numel(theseTrials)));
legend('observed','all mdl','vision','eye velocity','eye position','location','northwest')
%legend('observed','all mdl','vision','eye velocity','location','northwest')

end

