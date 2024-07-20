function [fig, stats] = showTonsetResp_stratifiedByCue(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t,  ...
    param, dd, validEvents)
% fig = showTonsetResp_stratified(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t,  param, dd, validEvents)

%only use successful trials without hesitation
%
% see also: getTgtNeuroLatency

nDiv = 2;%4;


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

%% behavioural latency for sorting trials
latency_bhv_srt = getTgtBhvLatency(onsets_cat, dd, validEvents, 1);

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
%p = signrank(bSignal, rSignal);
respPolarity = sign(mean(rSignal) - mean(bSignal));

theseColors = cool(nDiv+1);

nCols = 1;%2;
fig = figure('position',[1003         219         nCols*588         1200]);

stats = [];
dt = []; latency_neuro = []; thresh_neuro = [];
for istimtype = 1:nCols
    if istimtype==1
        %stimtrials = prefDirTrials;
        stimtrials = find(tgtDir == 0);
    elseif istimtype == 2
        %stimtrials = 1-prefDirTrials;
        stimtrials = find(tgtDir == 180);
    elseif istimtype == 3
        %stimtrials = prevDirTrials;
    elseif istimtype == 4
        %stimtrials = ones(size(prefDirTrials));
    end

    % theseTrials_c = find(success.*stimtrials);
    % [~, idx] = sort(latency_bhv_srt(theseTrials_c));
    % theseTrials = theseTrials_c(idx);
    
    if ceil(numel(stimtrials)/nDiv) <= 1
        continue;
    end

    %divIdx = 0;
    ampResp = []; factor = []; 
    multcomp_anova_stimtype = cell(nCols,1);
    for idiv = 1:nDiv+1
        %divIdx = divIdx(end)+1:ceil(idiv*numel(theseTrials)/nDiv);
        switch idiv
            % case 1
            %     divTrials = find(~isinf(diffCueFOnset));
            %     thisLabel = 'w cue';
            case 1
               divTrials = intersect(find(dd.cuedLoc(validEvents)==1), intersect(find(diffCueFOnset > .6), find(~isinf(diffCueFOnset))));
               thisLabel = 'congruent cue >.6';
            case 2
                divTrials = intersect(find(dd.cuedLoc(validEvents)==1), intersect(find(diffCueFOnset <= .6), find(~isinf(diffCueFOnset))));
                thisLabel = 'congruent cue <.6';
            % case 3
            %     divTrials = intersect(find(dd.cuedLoc(validEvents)==0), find(~isinf(diffCueFOnset)));
            %     thisLabel = 'incongruent';
            case nDiv+1
                divTrials = find(isinf(diffCueFOnset));
                thisLabel = 'wo cue';
        end
        
        theseTrials_div = intersect(intersect(find(success), stimtrials), divTrials);

        %% stats per division
        avgResp(idiv, istimtype, :) = median(singleResp_t(theseTrials_div,:),1);
        %avgLatency_bhv(idiv, istimtype) = median(latency_bhv_srt(theseTrials_div),1);
        seResp_t(idiv, istimtype, :) = ste(singleResp_t(theseTrials_div,:), 1);
        %seLatency_bhv(idiv, istimtype) = ste(latency_bhv_srt(theseTrials_div), 1);
        %avgAmpResp(idiv,istimtype,:) = median(mean(singleResp_t(theseTrials_div, rWin),2), 1);
        ampResp = [ampResp; mean(singleResp_t(theseTrials_div, rWin),2)];
        factor = [factor; idiv*ones(numel(theseTrials_div),1)];
    
        %% visualization
        axes(istimtype, idiv) = subplot(nDiv+3, nCols, nCols*(idiv-1)+istimtype, 'Parent', fig);
         
        if idiv==1
            if istimtype==1
                %title(['pref direction ' num2str(prefDir)]);
                title('stimDir = 0deg')
            elseif istimtype==2
                %title('other directions');
                title('stimDir = 180deg')
            elseif istimtype==3
                %title(['most frequent directions ' num2str(prevDir)]);
             elseif istimtype==4
                %title('all directions');
            end
         end

          if istimtype==1
            ylabel(thisLabel);
          end

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
    axes(istimtype, nDiv+2) = subplot(nDiv+3, nCols, nCols*(nDiv+1)+istimtype, 'Parent', fig);
    for idiv = 1:nDiv+1
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
    [p_anova, ~, stats_anova] = anova1(ampResp, factor, "off");
    [multcomp_anova, ~,~, gnames] = multcompare(stats_anova,"Display","off");

    multcomp_anova_stimtype{istimtype} = multcomp_anova;
    multcomp_anova_stimtype{istimtype}(:,1) = changem(multcomp_anova_stimtype{istimtype}(:,1), str2num(cell2mat(gnames))', 1:numel(gnames));
    multcomp_anova_stimtype{istimtype}(:,2) = changem(multcomp_anova_stimtype{istimtype}(:,2), str2num(cell2mat(gnames))', 1:numel(gnames));

    axes(istimtype, nDiv+3) = subplot(nDiv+3, nCols, nCols*(nDiv+2)+istimtype, 'Parent', fig);
    %violin(ampResp_c,'facecolor',theseColors, 'plotlegend', 0);
        simple_violin_scatter(factor, ampResp, 10, .5, {2, theseColors(factor,:)})
        title(sprintf('ANOVA p: %f', p_anova));
        addSignStar(axes(istimtype, nDiv+3), multcomp_anova_stimtype{istimtype});
   



end
linkaxes(axes(:),'y');
%linkaxes(axes(:,nDiv+1));

%% summary stats 
stats.latency = latency_neuro;
stats.thresh = thresh_neuro;
%stats.avglatency_bhv = avgLatency_bhv;
stats.avgResp = avgResp;
% stats.p_anova = p_anova;
stats.multcomp_anova = multcomp_anova_stimtype;

%stats.avgCorr = avgCorr;
%stats.pCorr = pCorr;
