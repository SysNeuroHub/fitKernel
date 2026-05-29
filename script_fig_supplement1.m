f = figure('position', [0 0 168 168]);
histogram('BinCounts', mean(saccDirNoTask_spkOkUCue_hist_pop,2), 'BinEdges', binEdges); hold on;
errorbar(param.cardinalDir, mean(saccDirNoTask_spkOkUCue_hist_pop,2), 1/sqrt(size(saccDirNoTask_spkOkUCue_hist_pop,2))*std(saccDirNoTask_spkOkUCue_hist_pop,[],2),'k.');
set(gca,'tickdir','out','XTick',param.cardinalDir);
ylabel('# saccades');xlabel('Direction [deg]');
box(gca,'off');