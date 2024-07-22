function [fig, stats] = showTonsetResp_stratifiedByCue(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t,  ...
    param, dd, validEvents)
% fig = showTonsetResp_stratified(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t,  param, dd, validEvents)

%only use successful trials without hesitation
%
% see also: getTgtNeuroLatency

%% delay between Target onset and cue Onset of each trial
diffCueFOnset = getDiffCueTgtOnset(onsets_cat, catEvTimes); %3/6/24
diffCueFOnset = diffCueFOnset(validEvents);

%% from showTonsetByCue:
onsetTimes_t = catEvTimes.tOnset(validEvents);
tgtDir = getTgtDir(dd.targetloc(validEvents), param.cardinalDir);

[~,dirIdx]=intersect(param.cardinalDir, unique(tgtDir));
% triggered by target onset
[~, winSamps_t, singleResp_t] ...
    = eventLockedAvg(PSTH_f', t_r, onsetTimes_t, tgtDir, tWin_t);

hesitant = getChoiceOutcome_hesitant(onsets_cat, dd, validEvents);

%% divide trials
[choiceOutcome] = getChoiceOutcome(dd);
success = (choiceOutcome(validEvents) == 1) .* ~hesitant;

% [prefDir, prefDirTrials_c] = getPrefDir(PSTH_f, t_r, onsetTimes_t, tgtDir, param);
% prefDirTrials = zeros(numel(onsetTimes_t),1);
% prefDirTrials(prefDirTrials_c) = 1;
% 
% [prevDir, prevDirTrials_c] = getPrevDir(tgtDir, param);
% prevDirTrials = zeros(numel(onsetTimes_t),1);
% prevDirTrials(prevDirTrials_c) = 1;

%% decide response polarity
bWin =  intersect(find(winSamps_t >= param.baseWin(1)), find(winSamps_t <= param.baseWin(2)));
rWin =  intersect(find(winSamps_t >= param.respWin(1)), find(winSamps_t <= param.respWin(2)));
bSignal = mean(singleResp_t(tgtDir == 0, bWin),2);
rSignal = mean(singleResp_t(tgtDir == 0, rWin),2);
p_visResp = signrank(bSignal, rSignal); %is there a modulation after target stimulation?
%respPolarity = sign(mean(rSignal) - mean(bSignal));

theseColors = cool(2); %wcue v wocue

nCols = 1;%2;
fig = figure('position',[1003         219         nCols*588         1200]);

stats = [];
dt = []; latency_neuro = []; thresh_neuro = [];
for istimtype = 1:nCols
    if istimtype==1
        %stimtrials = prefDirTrials;
        stimtrials = find(tgtDir == 0);
      end

    % theseTrials_c = find(success.*stimtrials);
    % [~, idx] = sort(latency_bhv_srt(theseTrials_c));
    % theseTrials = theseTrials_c(idx);
    
    if ceil(numel(stimtrials)) <= 1
        continue;
    end

    %divIdx = 0;
    ampResp = []; factor = []; 
    mampResp = nan(2,1);
    for idiv = 1:2
        %divIdx = divIdx(end)+1:ceil(idiv*numel(theseTrials)/nDiv);
        switch idiv
            % case 1
            %     divTrials = find(~isinf(diffCueFOnset));
            %     thisLabel = 'w cue';
            case 1
               divTrials = intersect(find(dd.cuedLoc(validEvents)==1), find(~isinf(diffCueFOnset)));
               thisLabel = 'congruent cue >.6';
            case 2
                divTrials = find(isinf(diffCueFOnset));
                thisLabel = 'wo cue';
        end
        
        theseTrials_div = intersect(intersect(find(success), stimtrials), divTrials);

        %% stats per division
        avgResp(idiv, istimtype, :) = median(singleResp_t(theseTrials_div,:),1);
        seResp_t(idiv, istimtype, :) = ste(singleResp_t(theseTrials_div,:), 1);
        
        thisAmpResp = mean(singleResp_t(theseTrials_div, rWin),2);
        ampResp = [ampResp; thisAmpResp];
        factor = [factor; idiv*ones(numel(theseTrials_div),1)];
        mampResp(idiv) = median(thisAmpResp);
    
        %% visualization
        axes(istimtype, idiv) = subplot(1+3, nCols, nCols*(idiv-1)+istimtype, 'Parent', fig);
         
        if istimtype==1&&idiv==1
           title(['stimDir = 0deg, p visResp: ' num2str(p_visResp)]);
        end
        ylabel(thisLabel);
        
        if sum(theseTrials_div)==0
            continue;
        end
        boundedline(winSamps_t, squeeze(avgResp(idiv,istimtype, :)),  squeeze(seResp_t(idiv, istimtype, :)), ...
            'cmap',theseColors(idiv,:));
        %vline(avgLatency_bhv(idiv, istimtype),  axes(istimtype,idiv), ':', 'k');
        xlabel(num2str(numel(theseTrials_div)));
        xlim(tWin_t);       
    end % end of idiv

    %% stats across divisions
    % [latency_neuro(:,istimtype), thresh_neuro(:,istimtype)] = ...
    %     getSingleTrialLatency(avgResp(:,istimtype,:), winSamps_t, tWin_t, param.threshParam);
    % 
    % for idiv = 1:nDiv
    %     hline(thresh_neuro(idiv,istimtype), axes(istimtype, idiv), ':', theseColors(idiv,:));
    % end

    %% summary across divisions
    axes(istimtype, 3) = subplot(1+3, nCols, nCols*2+istimtype, 'Parent', fig);
    for idiv = 1:1+1
        plot(winSamps_t,  squeeze(avgResp(idiv, istimtype, :)), 'Color', theseColors(idiv,:)); hold on;
    end
    xlim(tWin_t);

    % for idiv = 1:nDiv+1
    %     hold on;        vline(avgLatency_bhv(idiv, istimtype), axes(istimtype,nDiv+2),':',theseColors(idiv,:));
    % end
    set(gca,'box','off')
    %hline(thresh_neuro(:,istimtype));
    if istimtype==1
        ylabel('all divisions');
    end

    %% ANOVA comparing different cue conditions (across rows)
    % [p_anova, ~, stats_anova] = anova1(ampResp, factor, "off");
    % [multcomp_anova, ~,~, gnames] = multcompare(stats_anova,"Display","off");
    % 
    % multcomp_anova_stimtype{istimtype} = multcomp_anova;
    % multcomp_anova_stimtype{istimtype}(:,1) = changem(multcomp_anova_stimtype{istimtype}(:,1), str2num(cell2mat(gnames))', 1:numel(gnames));
    % multcomp_anova_stimtype{istimtype}(:,2) = changem(multcomp_anova_stimtype{istimtype}(:,2), str2num(cell2mat(gnames))', 1:numel(gnames));
    if  sum(factor==1)*sum(factor==2) ~= 0
        p_cueMod = ranksum(ampResp(factor == 1), ampResp(factor==2));
    else
        p_cueMod = NaN;
    end
    c = zeros(1,6);
    c(1:2) = [1 2];
    c(6) = p_cueMod;

    axes(istimtype, 4) = subplot(1+3, nCols, nCols*3+istimtype, 'Parent', fig);
    simple_violin_scatter(factor, ampResp, 10, .5, {2, theseColors(factor,:)})
    title(sprintf('p cueMod: %f', p_cueMod));
    addSignStar(axes(istimtype, 1+3), c);

end
linkaxes(axes(:),'y');
%linkaxes(axes(:,nDiv+1));

%% summary stats 
%stats.latency = latency_neuro;
%stats.thresh = thresh_neuro;
%stats.avglatency_bhv = avgLatency_bhv;
stats.mampResp = mampResp;
stats.p_visResp = p_visResp;
stats.p_cueModulation = p_cueMod;

