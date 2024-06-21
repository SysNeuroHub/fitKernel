function [fig, stats] = showTonsetResp_stratified(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t,  param, dd, validEvents)
% fig = showTonsetResp_stratified(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t,  param, dd, validEvents)

%only use successful trials without hesitation
%
% see also: getTgtNeuroLatency

nDiv = 2;

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
zeroDirTrials = zeros(numel(onsetTimes_t),1);
zeroDirTrials(tgtDir == 0) = 1;

theseColors = jet(nDiv);


fig = figure('position',[1003         219         3*588         680]);

stats = [];
dt = []; latency_neuro = []; thresh_neuro = [];
for istimtype = 1:3
    if istimtype==1
        stimtrials = prefDirTrials;
    elseif istimtype == 2
        stimtrials = 1-prefDirTrials;
    % elseif istimtype == 3
    %     stimtrials = zeroDirTrials;
    elseif istimtype == 3
        stimtrials = ones(size(prefDirTrials));
    end

    theseTrials_c = find(success.*stimtrials);
    [~, idx] = sort(latency_bhv_srt(theseTrials_c));
    theseTrials = theseTrials_c(idx);

    divIdx = 0;
    for idiv = 1:nDiv
        divIdx = divIdx(end)+1:ceil(idiv*numel(theseTrials)/nDiv);
        theseTrials_div = theseTrials(divIdx);

        %% stats per division
        avgResp_t(idiv,:) = mean(singleResp_t(theseTrials_div,:));
        avgLatency_bhv(idiv) = mean(latency_bhv_srt(theseTrials_div));

        %% visualization
        axes(istimtype, idiv) = subplot(nDiv+1, 3, 3*(idiv-1)+istimtype);
        if sum(theseTrials_div)==0
            continue;
        end
        plot(winSamps_t, singleResp_t(theseTrials_div,:), 'color',theseColors(idiv,:)); hold on;
        plot(winSamps_t, mean(singleResp_t(theseTrials_div,:)), 'color','k','LineWidth',2);
        vline(latency_bhv_srt(theseTrials_div), axes(istimtype,idiv),':',theseColors(idiv,:));
        vline(mean(latency_bhv_srt(theseTrials_div)), axes(istimtype,idiv),':','k');
        if istimtype==1
            ylabel([num2str(idiv) '/' num2str(nDiv)]);
        end
        if idiv==1
            if istimtype==1
                title('pref direction');
            elseif istimtype==2
                title('other directions');
             elseif istimtype==3
                title('all directions');
            end
        end
    end % end of idiv

    %% stats across divisions
    % currently only for nDiv == 2
    [latency_neuro(:,istimtype), thresh_neuro(:,istimtype)] = getSingleTrialLatency(avgResp_t, winSamps_t, tWin_t, 'uniform', 4);
    dt(istimtype) = diff(latency_neuro(:,istimtype));

    %% summary across divisions
    axes(istimtype, nDiv+1) = subplot(nDiv+1, 3, 3*nDiv+istimtype);
    for idiv = 1:nDiv
        plot(winSamps_t, avgResp_t(idiv,:), 'Color', theseColors(idiv,:)); hold on;
    end
    for idiv = 1:nDiv
        hold on;        vline(avgLatency_bhv(idiv), axes(istimtype,nDiv+1),':',theseColors(idiv,:));
    end
    hline(thresh_neuro(:,istimtype));
    title(['dt: ' num2str(dt(istimtype))]);

    if istimtype==1
        ylabel('all divisions');
    end
end
linkaxes(axes(:,1:nDiv));
linkaxes(axes(:,nDiv+1));

stats.latency = latency_neuro;
stats.thresh = thresh_neuro;
stats.dt = dt;
