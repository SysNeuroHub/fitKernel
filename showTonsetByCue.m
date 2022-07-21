function [f, avgTonsetByCue, winSamps_tonsetByCue] = showTonsetByCue(t_r, ...
    y_r, cardinalDir, catEvTimes, dd, psthNames, figTWin)
% [f, avgTOnsetByCue, winSamps_tonsetByCue] = showTonsetByCue(t_r, ...
%     y_r, cardinalDir, catEvTimes, dd, psthNames, figTWin)

% created from resp_cueConditions
onset = catEvTimes.tOnset;

% y_r = cat(1,PSTH_f',PSTH_ff', predicted_all, predicted, pdiam_r)

for icue = 1:2
    validEvents = intersect(find(~isnan(onset)), find(dd.cueOn==icue-1));
    %< this condition only includes sucessful saccades
    
    onsetTimes = onset(validEvents);
    tgtDir = getTgtDir(dd.targetloc(validEvents), cardinalDir);
    
    [~,dirIdx]=intersect(cardinalDir, unique(tgtDir));
    [avgTonsetByCue(dirIdx,:,:,icue), winSamps_tonsetByCue] ...
        = eventLockedAvg(y_r', t_r, onsetTimes, tgtDir, figTWin);
end

%% w vs wo cue
colormap('parula');
f = figure('position',[0 0 1400 1000]);
for icue = 1:3
    nvars = size(avgTonsetByCue,2);
    
    switch icue
        case {1,2}
            crange = prctile(avgTonsetByCue(:),[1 99]);
        case 3
            cache = abs(diff(avgTonsetByCue,1,4));
            crange = [-prctile(cache(:),99) prctile(cache(:),99)];
    end
    for ivar = 1:nvars
        if icue <3
            thisData = squeeze(avgTonsetByCue(:,ivar,:,icue));
        elseif icue ==3
            thisData = squeeze(diff(avgTonsetByCue(:,ivar,:,:),1,4));
        end
        subplot(nvars, 3, 3*(ivar-1)+icue);
        imagesc(winSamps_tonsetByCue, cardinalDir, thisData);
        set(gca, 'ytick',cardinalDir);
        if icue==1
            ylabel(psthNames{ivar});
        end
        caxis(crange);
        %mcolorbar(gca,.5);
        if ivar == 1
            if icue==1
                title('wo cue');
            elseif icue ==2
                title('w cue');
            elseif icue == 3
                title('wcue-wocue');
            end
        end
    end
    mcolorbar;
end
xlabel(['time from tOnset [s]']);
