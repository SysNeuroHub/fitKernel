%load('/mnt/syncitium/Daisuke/cuesaccade_data/2022/hugo/eyeCat_hugo09September_06.mat')
%load('/mnt/syncitium/Daisuke/cuesaccade_data/2022/hugo/hugo09September_06_10_linear_rReg.mat')

function [fig] = ...
    getTonsetByCue(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, param, dd, validEvents)
% created from getTgtNeuroLatency(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, Thresh, param, dd)
% see also. showTonsetByCue


%% from showTonsetByCue:
onsetTimes_t = catEvTimes.tOnset(validEvents);
tgtDir = getTgtDir(dd.targetloc(validEvents), param.cardinalDir);

[~,dirIdx]=intersect(param.cardinalDir, unique(tgtDir));
% triggered by target onset
[~, winSamps_t, singleResp_t] ...
    = eventLockedAvg(PSTH_f', t_r, onsetTimes_t, tgtDir, tWin_t);
singleResp_t = reshape(singleResp_t, size(singleResp_t,1), size(singleResp_t,3));
%singleResp: trials x times


%% behavioural latency for sorting trials
latency_bhv_srt = getTgtBhvLatency(onsets_cat, dd, validEvents, 1);
latency_bhv_cOn = getTgtBhvLatency(onsets_cat, dd, validEvents, 0);

%% delay between Target onset and cue Onset of each trial
diffCueFOnset = getDiffCueTgtOnset(onsets_cat, catEvTimes); %3/6/24

%% divide trials
[choiceOutcome] = getChoiceOutcome(dd);
success = (choiceOutcome(validEvents) == 1) .* (latency_bhv_cOn - latency_bhv_srt <= 0.1);
hesitant = (latency_bhv_cOn - latency_bhv_srt > 0.1); %1st saccade did NOT lead to choice
fail = choiceOutcome(validEvents) >= 2; %quiescent + wrong

[prefDir, prefDirTrials_c] = getPrefDir(PSTH_f, t_r, onsetTimes_t, tgtDir, param);
prefDirTrials = zeros(numel(onsetTimes_t),1);
prefDirTrials(prefDirTrials_c) = 1;

[prevDir, prevDirTrials_c] = getPrevDir(tgtDir, param);
prevDirTrials = zeros(numel(onsetTimes_t),1);
prevDirTrials(prevDirTrials_c) = 1;


%% visualization
 [~, latency_bhv_sortedTrIdx] = sort(latency_bhv_srt);

nDivs = 3;
fig = figure('position',[1003         219         nDivs*588         664]);

divIdx = 0;
for idiv = 1:nDivs
    divIdx = divIdx(end)+1:ceil(idiv*numel(prefDirTrials)/nDivs);

    %divTrials = zeros(numel(prefDirTrials),1);
    %divTrials(latency_bhv_sortedTrIdx(divIdx)) = 1;
   divTrialIdx = latency_bhv_sortedTrIdx(divIdx);
    % thesePrefDirTrials = prefDirTrials(find(divTrials));
    %[~, idx] = sort(diffCueFOnset(find(thesePrefDirTrials)));
    %theseTrials = thesePrefDirTrials(idx);
    %nTrials =  [sum(success.*stimtrials) sum(hesitant.*stimtrials) sum(fail.*stimtrials)];
    [~, idx] = sort(diffCueFOnset(divTrialIdx));
    theseTrials = divTrialIdx(idx);

    axes(idiv) = subplot(1,nDivs,idiv);
    if sum(theseTrials)==0
        continue;
    end
    imagesc(winSamps_t,1:numel(theseTrials), singleResp_t(theseTrials,:));
    hold on
    %plot(latency_neuro(theseTrials), 1:numel(theseTrials),'r.');
    plot(-diffCueFOnset(theseTrials), 1:numel(theseTrials),'g.');
    plot(latency_bhv_srt(theseTrials), 1:numel(theseTrials),'w.');
    vline(0);
    %hline(cumsum(nTrials)+.5,gca,'-','m');
    title(['tgt stim on preferred direction ' num2str(prefDir) ]);
    xlabel('time from target onset [s]');
    ylabel('fail/ hesitant / success');
end
linkcaxes(axes);mcolorbar(gca, .5);
