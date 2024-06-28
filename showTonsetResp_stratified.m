function [fig, stats] = showTonsetResp_stratified(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t,  ...
    param, dd, validEvents, nDiv)
% fig = showTonsetResp_stratified(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t,  param, dd, validEvents)

%only use successful trials without hesitation
%
% see also: getTgtNeuroLatency

if nargin < 9
    nDiv = 2;
end

%% from showTonsetByCue:
onsetTimes_t = catEvTimes.tOnset(validEvents);
tgtDir = getTgtDir(dd.targetloc(validEvents), param.cardinalDir);

[~,dirIdx]=intersect(param.cardinalDir, unique(tgtDir));
% triggered by target onset
[~, winSamps_t, singleResp_t] ...
    = eventLockedAvg(PSTH_f', t_r, onsetTimes_t, tgtDir, tWin_t);

%% behavioural latency for sorting trials
latency_bhv_srt = getTgtBhvLatency(onsets_cat, dd, validEvents, 1);

hesitant = getChoiceOutcome_hesitant(onsets_cat, dd, validEvents);

%% divide trials
[choiceOutcome] = getChoiceOutcome(dd);
success = (choiceOutcome(validEvents) == 1) .* ~hesitant;
%hesitant = (latency_bhv_cOn - latency_bhv_srt > 0.1); %1st saccade did NOT lead to choice
%fail = choiceOutcome(validEvents) >= 2; %quiescent + wrong

[prefDir, prefDirTrials_c] = getPrefDir(PSTH_f, t_r, onsetTimes_t, tgtDir, param);
prefDirTrials = zeros(numel(onsetTimes_t),1);
prefDirTrials(prefDirTrials_c) = 1;

[prevDir, prevDirTrials_c] = getPrevDir(tgtDir, param);
prevDirTrials = zeros(numel(onsetTimes_t),1);
prevDirTrials(prevDirTrials_c) = 1;


theseColors = cool(nDiv);


nCols = 4;
fig = figure('position',[1003         219         nCols*588         1200]);

stats = [];
dt = []; latency_neuro = []; thresh_neuro = [];
for istimtype = 1:nCols
    if istimtype==1
        stimtrials = prefDirTrials;
    elseif istimtype == 2
        stimtrials = 1-prefDirTrials;
    elseif istimtype == 3
        stimtrials = prevDirTrials;
    elseif istimtype == 4
        stimtrials = ones(size(prefDirTrials));
    end

    theseTrials_c = find(success.*stimtrials);
    [~, idx] = sort(latency_bhv_srt(theseTrials_c));
    theseTrials = theseTrials_c(idx);

    if ceil(numel(theseTrials)/nDiv) <= 1
        continue;
    end

    divIdx = 0;
    for idiv = 1:nDiv
        divIdx = divIdx(end)+1:ceil(idiv*numel(theseTrials)/nDiv);
        theseTrials_div = theseTrials(divIdx);

        %% stats per division
        avgResp(idiv, istimtype, :) = median(singleResp_t(theseTrials_div,:),1);
        avgLatency_bhv(idiv, istimtype) = median(latency_bhv_srt(theseTrials_div),1);
        seResp_t(idiv, istimtype, :) = ste(singleResp_t(theseTrials_div,:), 1);
        seLatency_bhv(idiv, istimtype) = ste(latency_bhv_srt(theseTrials_div), 1);


        %% visualization
        axes(istimtype, idiv) = subplot(nDiv+2, nCols, nCols*(idiv-1)+istimtype);
        if sum(theseTrials_div)==0
            continue;
        end
        boundedline(winSamps_t, squeeze(avgResp(idiv,istimtype, :)),  squeeze(seResp_t(idiv, istimtype, :)), ...
            'cmap',theseColors(idiv,:));
        vline(avgLatency_bhv(idiv, istimtype),  axes(istimtype,idiv), ':', 'k');

        if istimtype==1
            ylabel([num2str(idiv) '/' num2str(nDiv)]);
        end
        if idiv==1
            if istimtype==1
                title(['pref direction ' num2str(prefDir)]);
            elseif istimtype==2
                title('other directions');
            elseif istimtype==3
                title(['most frequent directions ' num2str(prevDir)]);
             elseif istimtype==4
                title('all directions');
            end
        end
    end % end of idiv

    %% stats across divisions
    [latency_neuro(:,istimtype), thresh_neuro(:,istimtype)] = ...
        getSingleTrialLatency(avgResp(:,istimtype,:), winSamps_t, tWin_t, param.threshParam);
  
    for idiv = 1:nDiv
        hline(thresh_neuro(idiv,istimtype), axes(istimtype, idiv), ':', theseColors(idiv,:));
    end

    %% summary across divisions
    axes(istimtype, nDiv+1) = subplot(nDiv+2, nCols, nCols*nDiv+istimtype);
    for idiv = 1:nDiv
        plot(winSamps_t,  squeeze(avgResp(idiv, istimtype, :)), 'Color', theseColors(idiv,:)); hold on;
    end
    for idiv = 1:nDiv
        hold on;        vline(avgLatency_bhv(idiv, istimtype), axes(istimtype,nDiv+1),':',theseColors(idiv,:));
    end
    %hline(thresh_neuro(:,istimtype));
    if istimtype==1
        ylabel('all divisions');
    end

    axes(istimtype, nDiv+2) = subplot(nDiv+2, 4, 4*(nDiv+1)+istimtype);
    scatter(avgLatency_bhv(:,istimtype), latency_neuro(:,istimtype), [],theseColors)
    nonan = ~isnan(latency_neuro(:,istimtype));
    if sum(nonan)==0
        avgCorr(istimtype) = nan; pCorr(istimtype) = nan;
    else
        [avgCorr(istimtype), pCorr(istimtype)] = corr(avgLatency_bhv(nonan,istimtype), latency_neuro(nonan,istimtype), 'type','Spearman');
    end
    title(['rank corr: ' num2str(avgCorr(istimtype)) ', p: ' num2str(pCorr(istimtype))]);
    squareplot; axis padded;


end
linkaxes(axes(:,1:nDiv));
linkaxes(axes(:,nDiv+1));

%% summary stats 
stats.latency = latency_neuro;
stats.thresh = thresh_neuro;
stats.avglatency_bhv = avgLatency_bhv;
stats.avgResp = avgResp;
stats.avgCorr = avgCorr;
stats.pCorr = pCorr;
