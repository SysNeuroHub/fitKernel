%% how ROC analysis works:
hit = 10+randn(100,1);
miss = 11+randn(120,1);

labels = [zeros(1,numel(hit)) ones(1,numel(miss))];
scores = [hit; miss]';
[X,Y,~,auc] = perfcurve(labels, scores, 1);

f = figure('position',[0 0 2*168 168]);
subplot(121);
histogram(hit, 'binLimits',prctile(scores(:),[0 100]),'NumBins',15);hold on;
histogram(miss, 'binLimits',prctile(scores(:),[0 100]),'NumBins',15);hold on;
axis square;
box off;
set(gca,'tickdir','out');
ylabel('# trials'); xlabel('Firing rate [spikes/s]');
legend('Correct','False');

subplot(122);
plot(X,Y);squareplot;
xlabel('False positive rate')
ylabel('True positive rate')
title(['AUC:' num2str(auc)])
%title('ROC for Classification by Logistic Regression');
% rocObj = rocmetrics(species(51:end,:),scores,'virginica');
% plot(rocObj)
savePaperFigure(f,'figure_supplement3_1');
