function f = function_fig_supplement2_1(corr_tgt_avg_pop, animalid_pop, param)
f = figure('position', [0 0 200 7*100]);
for imodality = 1:6
    ax(imodality) = subplot(6,1,imodality);
    subplot(6,1,imodality);
    histogram(corr_tgt_avg_pop(imodality,animalid_pop==1),linspace(-1,1,20),'FaceColor',[.5 .5 .5]); hold on;
    histogram(corr_tgt_avg_pop(imodality,animalid_pop==2),linspace(-1,1,20),'FaceColor',[.1 .1 .1]);
    axis tight
    box off
    set(gca,'tickdir','out');
    if imodality ==6
        xlabel('correlation'); ylabel('# units');
    elseif imodality == 1
        legend('M1','M2','Location','northwest');
        title('Full model');
    end
    if imodality > 1
        title(param.predictorNames{imodality-1});
    end
    %title(['mean: ' num2str(nanmean(corr_tgt_avg_pop(imodality, :))) ', s.d.:' num2str(nanstd(corr_tgt_avg_pop(imodality, :)))]);
    set(gca,'ytick',0:100:400);
end
linkaxes(ax);