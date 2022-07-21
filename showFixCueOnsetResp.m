function [f, avgfOnsetResp, avgCueResp, winSamps] = showFixCueOnsetResp(t_r, ...
    y_r, catEvTimes, dd, psthNames, figTWin)
% f = showFixCueOnsetResp(t_r, y_r, catEvTimes, dd, predictorNames, respWin)
% respWin = [-0.5 1];

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
    imagesc(winSamps, 1:nRespTypes, squeeze(avgfOnsetResp(icue,:,:)));
    if icue == 1
    set(gca,'yticklabel',psthNames);%[{'observed'},{'predicted_all'},predictorNames(:)']);
    end
    vline(0);
    xlabel('Time from fOnset [s]');
    title(['fixOn ' condName ': ' num2str(sum(dd.cueOn(validTrials)==icue-1)) ' trials']);
    caxis(gca,crange);
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
    
    subplot(1,3,3);
    imagesc(winSamps, 1:nRespTypes, squeeze(avgCueResp));
    %set(gca,'yticklabel',[{'observed'},{'predicted_all'},param.predictorNames(:)']);
    vline(0);
    caxis(gca,crange);
    xlabel('Time from cueOnset [s]');
    title(['cueOn']);
    mcolorbar(gca,.5);
    %             screen2png(fullfile(saveFigFolder,['cueOn_' saveSuffix]));
    %             close;
else
    avgCueResp = nan(size(avgfOnsetResp));
end