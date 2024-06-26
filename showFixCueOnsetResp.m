function [f, avgfOnsetResp, avgCueResp, winSamps] = showFixCueOnsetResp(t_r, ...
    y_r, catEvTimes, dd, psthNames, figTWin, diff_at_zero)
% f = showFixCueOnsetResp(t_r, y_r, catEvTimes, dd, predictorNames, respWin)
% respWin = [-0.5 1];

if nargin < 7
    diff_at_zero = 0;
end

nRespTypes = size(y_r,2);

%% response to fixation w/wo cue
%cueOnset = catEvTimes.cueOnset;
fOnset = catEvTimes.fOnset;
validTrials = find(~isnan(fOnset));

[avgfOnsetResp, winSamps,~,~,uniqueLabels] ...
    = eventLockedAvg(y_r', t_r, fOnset(validTrials), dd.cueOn(validTrials), ...
    figTWin);%param.figTWin);

if numel(uniqueLabels)==1
    for icue = 1:2
        if icue == uniqueLabels+1
            avgfOnsetResp_c(icue,:,:) = avgfOnsetResp;
        else
            avgfOnsetResp_c(icue,:,:) = nan(size(avgfOnsetResp));
        end
    end
    avgfOnsetResp = avgfOnsetResp_c;
end

if diff_at_zero
    [~,zeroidx] = min(abs(winSamps - 0));
    avgfOnsetResp = avgfOnsetResp - avgfOnsetResp(:,:,zeroidx);
end

crange = prctile(avgfOnsetResp(:),[0 100]);
f = figure('position',[0 0 1000 500]);
for icue = 1:2
    switch icue
        case 1
            condName = 'wo cue';
        case 2
            condName = 'w cue';
    end
    subplot(1,3,icue);
    if nRespTypes ==1
        plot(winSamps, squeeze(avgfOnsetResp(icue,:,:)));
        ylabel(psthNames);
    else
        imagesc(winSamps, 1:nRespTypes, squeeze(avgfOnsetResp(icue,:,:)));
        if icue == 1
            set(gca,'yticklabel',psthNames);%[{'observed'},{'predicted_all'},predictorNames(:)']);
        end
    end
    vline(0);
    xlabel('Time from fixation behaviour Onset [s]');
    title(['fixOn ' condName ': ' num2str(sum(dd.cueOn(validTrials)==icue-1)) ' trials']);
    if nRespTypes == 1
        ylim(gca,crange);
    else
        caxis(gca,crange);
    end
end
% mcolorbar;
% screen2png(fullfile(saveFigFolder,['fOnset_' saveSuffix]));
% close;


%% response to cueOnset
cueOnset = catEvTimes.cueOnset;
fOnset = catEvTimes.fOnset;
validEvents = intersect(find(dd.cueOn==1), find(~isnan(fOnset)));
cueOnset = cueOnset(validEvents);
%cuedLoc = dd.cuedLoc(validEvents);

if ~isempty(validEvents)
    [avgCueResp, winSamps] ...
        = eventLockedAvg(y_r', t_r, cueOnset, ones(numel(cueOnset),1), figTWin);%param.figTWin);
    
    if diff_at_zero
        [~,zeroidx] = min(abs(winSamps - 0));
        avgCueResp = avgCueResp - avgCueResp(:,:,zeroidx);
    end

    subplot(1,3,3);
    if nRespTypes ==1
        plot(winSamps, squeeze(avgCueResp));
    else
        imagesc(winSamps, 1:nRespTypes, squeeze(avgCueResp));
    end
    vline(0);
    if nRespTypes == 1
        ylim(gca,crange);
    else
        caxis(gca,crange);
        mcolorbar(gca,.5);
    end
    xlabel('Time from cueOnset [s]');
    title(['cueOn']);
else
    avgCueResp = nan(size(avgfOnsetResp));
end