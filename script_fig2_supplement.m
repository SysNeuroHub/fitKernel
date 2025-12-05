f = figure('position', [0 0 200 200]);
histogram(corr_tgt_avg_pop(1,animalid_pop==1),linspace(param.corr_tgtTh,1,20),'FaceColor',[.5 .5 .5]); hold on;
histogram(corr_tgt_avg_pop(1,animalid_pop==2),linspace(param.corr_tgtTh,1,20),'FaceColor',[.1 .1 .1])
axis tight square 
box off
set(gca,'tickdir','out')
xlabel('correlation of the full model'); ylabel('# units');
legend('M1','M2','Location','northwest');
title(['mean: ' num2str(mean(corr_tgt_avg_pop(1, :))) ', s.d.:' num2str(std(corr_tgt_avg_pop(1, :)))]);