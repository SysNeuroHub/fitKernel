function [latency_neuro, thresh_neuro, tgtDir, fig, nlatencyTrials_pref_success, latencyStats] = ...
    getTgtNeuroLatency(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, ThreshParam, param, dd, validEvents)
% [latency_neuro, validEvents, thresh_neuro, tgtDir] = getTgtNeuroLatency(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, Thresh, param, dd)
% create figure of single-trials sorted by behavioural latency of
% individual neurons

%16/1/25 added nlatencyTrials_pref_success and latencyStats 

markerSize = 10;
figPosition = [0 0 163 2*163];

if isempty(ThreshParam)
    %threshOption = 'uniform'; %uniform threshold across all trials or inidividual threshold of each trial (India)
    ThreshParam.option = 'individual'; %OR  define by onset of behavioural fixation as 'fonset'; turned out to produce too short latency
    ThreshParam.dur = 0.04;%s
    ThreshParam.Thresh = 4;
end


%% target response
%baseline = 'tonset';
%tWin_f = [0  min(onsets_cat.cueOnset-onsets_cat.fOnset)];

%% delay between Target onset and cue Onset of each trial
diffCueFOnset = getDiffCueTgtOnset(onsets_cat, catEvTimes); %3/6/24
diffCueFOnset = diffCueFOnset(validEvents);

%% from showTonsetByCue:
onsetTimes_t = catEvTimes.tOnset(validEvents);
tgtDir = getTgtDir(dd.targetloc(validEvents), param.cardinalDir);
% onsetTimes_t = catEvTimes.tOnset; %13/12/24
% tgtDir = getTgtDir(dd.targetloc, param.cardinalDir); %13/12/24

[~,dirIdx]=intersect(param.cardinalDir, unique(tgtDir));
% triggered by target onset
[~, winSamps_t, singleResp_t] ...
    = eventLockedAvg(PSTH_f', t_r, onsetTimes_t, tgtDir, tWin_t);
singleResp_t = reshape(singleResp_t, size(singleResp_t,1), size(singleResp_t,3));
singleResp = nan(numel(validEvents), numel(winSamps_t));
singleResp(~isnan(onsetTimes_t),:) = singleResp_t; %16/1/25
%singleResp: trials x times (including events w onsetTimes_t == NAN)



%% neural latency
[latency_neuro, thresh_neuro] = getSingleTrialLatency(singleResp, winSamps_t, tWin_t, ThreshParam);
% from     latencyV1MT/secondary/findSingleTrialLatency.m

%% behavioural latency for sorting trials
latency_bhv_srt = getTgtBhvLatency(onsets_cat, dd, validEvents, 1);
latency_bhv_cOn = getTgtBhvLatency(onsets_cat, dd, validEvents, 0);

%% divide trials
[choiceOutcome] = getChoiceOutcome(dd);
success = (choiceOutcome(validEvents) == 1) .* (latency_bhv_cOn - latency_bhv_srt <= 0.1); %1st saccade lead to choice
fail = choiceOutcome(validEvents) >= 2; %quiescent + wrong
hesitant = (latency_bhv_cOn - latency_bhv_srt > 0.1); %1st saccade did NOT lead to choice

[prefDir, prefDirTrials_c] = getPrefDir(PSTH_f, t_r, onsetTimes_t, tgtDir, param); %same as showTonsetResp
prefDirTrials = zeros(numel(onsetTimes_t),1);
prefDirTrials(prefDirTrials_c) = 1;
clear prefDirTrials_c


%% visualization
if nargout>3
    % fig = figure('position',[1003         300        482         664]);
    fig = figure('position', figPosition);

    stimtrials = prefDirTrials;

    theseTrials = []; nTrials = []; difflatency = nan;
    for ievtype = 1%:3
        switch ievtype
            case 1
                theseTrials_c = find(success.*stimtrials);
            case 2
                theseTrials_c = find(hesitant.*stimtrials);%find(quiescent.*stimtrials);
            case 3
                theseTrials_c = find(fail.*stimtrials);
        end

        theseTrials_c = intersect(theseTrials_c, validEvents);
        bothidx = find(~isnan(latency_bhv_srt(theseTrials_c)).*~isnan(latency_neuro(theseTrials_c)));
        theseTrials_c = theseTrials_c(bothidx);

        [~, idx] = sort(latency_bhv_srt(theseTrials_c));
        theseTrials = [theseTrials; theseTrials_c(idx)];
        nTrials(ievtype) = numel(theseTrials_c);
    
        %from getTgtlatencyCorr.m
        corrOption = 'Spearman'; %'Pearson'
        try
            [r(ievtype), p(ievtype)] = corr(latency_bhv_srt(theseTrials_c), latency_neuro(theseTrials_c), 'type', corrOption);
        catch err
            r(ievtype) = nan;
            p(ievtype) = nan;
        end

        axes(2*ievtype-1) = subplot(3, 1, 1);
        if sum(theseTrials)==0
            continue;
        end
        imagesc(winSamps_t, 1:numel(theseTrials_c), singleResp(theseTrials_c(idx),:));
        colormap(1-gray);
        hold on
        plot(latency_bhv_srt(theseTrials_c(idx)), 1:numel(theseTrials_c(idx)),'m.','MarkerSize',markerSize);
        plot(latency_neuro(theseTrials_c(idx)), 1:numel(theseTrials_c(idx)),'c.','MarkerSize', markerSize);
        legend('Behavioural','Neural','Location','northwest');
        %plot(-diffCueFOnset(theseTrials_c(idx)), 1:numel(theseTrials_c(idx)),'g.');
        %vline(0);

        difflatency = median(latency_neuro(theseTrials_c(idx)) - latency_bhv_srt(theseTrials_c(idx)));
        if ievtype == 1
            title(['tgt stim on preferred direction ' num2str(prefDir) ]);
            xlabel('Time from target onset [s]')
            %ylabel('success'); 
            ylabel('Trial ID (sorted)');
        elseif ievtype == 2
            ylabel('hesitent');
        elseif ievtype == 3
            ylabel('failure');
        end
        set(axes(2*ievtype-1),'tickdir','out');
        data_c = singleResp(theseTrials_c(idx),:);
        clim([max(0, min(data_c(:))) 0.1*round(10*max(data_c(:)))]);
        mcolorbar(gca, .5);

        %% scatter plot
        axes(2*ievtype) = subplot(3,1, [2 3]);
        scatter(latency_bhv_srt(theseTrials_c(idx)), latency_neuro(theseTrials_c(idx)),markerSize,'MarkerEdgeColor','k'); hold on;
        xlabel('Behavioural');
        ylabel('Neural');
        ax = gca;
        ax.XColor = 'm'; ax.YColor = 'c';
        tname = sprintf('corr: %.2f (p=%.2f)', r(ievtype), p(ievtype));
        title(tname); %title(['r=' num2str(r) ', r(success)=' num2str(r_success)]);
        squareplot;
        box off;
        axis padded; set(gca,'TickDir','out');

    end
end
%linkcaxes(axes([1 3 5]));

nlatencyTrials_pref_success = nTrials(1);
latencyStats.latency_r = r(1);
latencyStats.latency_p = p(1);
latencyStats.difflatency = difflatency;

