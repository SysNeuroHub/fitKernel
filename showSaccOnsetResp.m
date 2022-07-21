function [f, avgSaccResp, winSamps_sacc, singleSaccResp, sortedSaccLabels] = ...
    showSaccOnsetResp(t_r, y_r, cardinalDir, psthNames, ...
    startSaccNoTask, saccDirNoTask, figTWin)
%f = showSaccOnsetResp(t_r, y_r, catEvTimes, cardinalDir, psthNames, blinks, outliers, figTWin)

% tOnset = catEvTimes.tOnset;
% cOnset = catEvTimes.cOnset;
% validEvents = intersect(find(~isnan(tOnset)), find(~isnan(cOnset)));
% tOnset = tOnset(validEvents);
% cOnset = cOnset(validEvents);
% 
% tcOnset_trace = event2Trace(t_cat, [tOnset; cOnset], 2*0.5);
% excEventT_cat = (tcOnset_trace + blinks + outliers > 0); %28/1/22
% 
% [startSaccNoTask, endSaccNoTask] = selectSaccades(catEvTimes.saccadeStartTimes, ...
%     catEvTimes.saccadeEndTimes, t_cat, excEventT_cat);%param.minSaccInterval);
% 
% [saccDirNoTask, dirIndexNoTask] = getSaccDir(startSaccNoTask, endSaccNoTask, ...
%     eyeData_rmotl_cat, param.cardinalDir);

[avgSaccResp, winSamps_sacc, singleSaccResp, sortedSaccLabels] ...
    = eventLockedAvg(y_r', t_r, startSaccNoTask, saccDirNoTask, figTWin);

crange = prctile(avgSaccResp(:),[1 99]);
nvars = size(avgSaccResp,2);
f = figure('position',[0 0 400 1000]);
for ivar = 1:nvars
    subplot(nvars, 1, ivar);
    imagesc(winSamps_sacc, cardinalDir, squeeze(avgSaccResp(:,ivar,:)));
    set(gca, 'ytick',cardinalDir);
    xlabel('time from saccade onset [s]');
    ylabel(psthNames{ivar});
    caxis(gca,crange);
end
mcolorbar(gca,.5);

%         screen2png(fullfile(saveFigFolder,['saccOn_' saveSuffix]));
%         close;
