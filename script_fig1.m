f = figure('position', [0 0 168 168]);
histogram(prefDir_pop(animalid_pop == 1), 0:5:359,'FaceColor',[.5 .5 .5]); hold on;
histogram(prefDir_pop(animalid_pop == 2), 0:5:359,'FaceColor',[.1 .1 .1]);
axis tight square; box off;
ylim([0 35])
set(gca,'tickdir','out','xtick',[0 90 180 270 360]);
title(['n=' num2str(numel(animalid_pop))]);
xlabel('Target stimulus direction [deg]'); ylabel('# Units');
legend('M1','M2');