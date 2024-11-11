function [fig, avgAmp, p_hm] = showTonsetResp_hm(y_r, t_r, onsets_cat, catEvTimes, tWin_t,  ...
    param, dd, validEvents)
% fig = showTonsetResp_hm(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t,  param, dd, validEvents)
% 
%

%% from showTonsetByCue:
onsetTimes_t = catEvTimes.tOnset(validEvents);
tgtDir = getTgtDir(dd.targetloc(validEvents), param.cardinalDir);

[~,dirIdx]=intersect(param.cardinalDir, unique(tgtDir));
% triggered by target onset
[~, winSamps_t, singleResp_t] ...
    = eventLockedAvg(y_r', t_r, onsetTimes_t, tgtDir, tWin_t);

hesitant = getChoiceOutcome_hesitant(onsets_cat, dd, validEvents);

%% divide trials
[choiceOutcome] = getChoiceOutcome(dd);
success = (choiceOutcome(validEvents) == 1);% .* ~hesitant;
failure = (choiceOutcome(validEvents) ~= 1);

[prefDir, prefDirTrials_c] = getPrefDir(y_r(:,1), t_r, onsetTimes_t, tgtDir, param);
% prefDirTrials = zeros(numel(onsetTimes_t),1);
% prefDirTrials(prefDirTrials_c) = 1;

%% decide response polarity
% bWin =  intersect(find(winSamps_t >= param.baseWin(1)), find(winSamps_t <= param.baseWin(2)));
rWin =  intersect(find(winSamps_t >= param.respWin(1)), find(winSamps_t <= param.respWin(2)));
% bSignal = mean(singleResp_t(tgtDir == 0, bWin),2);
% rSignal = mean(singleResp_t(tgtDir == 0, rWin),2);
% p_visResp = signrank(bSignal, rSignal); %is there a modulation after target stimulation?
%respPolarity = sign(mean(rSignal) - mean(bSignal));

theseColors = [1 0 0; 0 0 1]; %h v m

fig = figure('position',[1003         219         2*588         1200]);

avgResp = nan(numel(param.cardinalDir), 2,2, numel(winSamps_t));
seResp = nan(numel(param.cardinalDir), 2,2, numel(winSamps_t));
avgAmp = nan(numel(param.cardinalDir), 2,2);
singleAmp = cell(numel(param.cardinalDir), 2,2);
for iregress = 1:2
    for ihm = 1:2
        switch ihm
            case 1
                stimTrials_c = find(success); %intersect(intersect(find(success), stimtrials), divTrials);
            case 2
                stimTrials_c = find(failure); %intersect(intersect(find(failure), stimtrials), divTrials);
        end

        for itgtDir = 1:numel(param.cardinalDir)
            theseTrials = intersect(find(tgtDir == param.cardinalDir(itgtDir)), stimTrials_c);

            if ceil(numel(theseTrials)) <= 1
                continue;
            end

            %% stats per division
            if iregress == 1 %observed
                theseData = squeeze(singleResp_t(theseTrials,1,:));
                tname = 'observed';
            elseif iregress == 2 %subtracted away by the full model
                theseData = -squeeze(diff(singleResp_t(theseTrials,:,:),1,2));
                tname = 'after subtraction by full mdl';
            end
                avgResp(itgtDir, ihm, iregress, :) = mean(theseData,1);
                seResp(itgtDir, ihm, iregress, :) = ste(theseData, 1);

                avgAmp(itgtDir, ihm, iregress) = mean(mean(theseData(:,rWin),1));
                singleAmp{itgtDir, ihm, iregress} = mean(theseData(:,rWin),2);
        end % end of ihm

        %% visualization
        %% time course at preferred direction
        axes(1,iregress) = subplot(2, 2, iregress, 'Parent', fig);

        prefDirIdx = find(param.cardinalDir == prefDir);
        boundedline(winSamps_t, squeeze(avgResp(prefDirIdx, ihm,iregress, :)),  ...
            squeeze(seResp(prefDirIdx, ihm, iregress, :)), 'color', theseColors(ihm,:),'transparency', 0.5); hold on;
        %vbox(param.respWin(1), param.respWin(2));
        %vline(param.respWin, ':', 'k');
        xlabel(num2str(numel(theseTrials)));
        xlim(tWin_t);
        title(tname);
        if iregress == 2
            legend('success','failure');
        end

        %% circular distribution
        axes(2,iregress) = subplot(2, 2, iregress+2, 'Parent', fig);
        polarplot(pi/180*param.cardinalDir, squeeze(avgAmp(:,ihm,iregress))', 'color', theseColors(ihm,:)); hold on
        polarscatter(pi/180*param.cardinalDir(prefDirIdx), squeeze(avgAmp(prefDirIdx,ihm,iregress)),'color', theseColors(ihm,:));
    end

    %% hit v miss stats
    if ~isempty(singleAmp{prefDirIdx, 1, iregress}) && ~isempty(singleAmp{prefDirIdx, 2, iregress})
        p_hm(iregress) = ranksum(singleAmp{prefDirIdx, 1, iregress}, singleAmp{prefDirIdx, 2, iregress});
    else
        p_hm(iregress) = NaN;
    end
    title(['p = ' num2str(p_hm(iregress))])

end
linkaxes(axes(1,:));

%linkaxes(axes(2,:));%NG as axes inputs must be Cartesian.
rlim(axes(2,:), [0 max(avgAmp(:))]);



%% summary stats
% stats.mampResp = mampResp;
% stats.p_visResp = p_visResp;
% stats.p_cueModulation = p_cueMod;

