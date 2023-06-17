
yaxis = linspace(0, 315, 8);

subplot(18,1,[1 8]);
imagesc(t_r(1:1e3), yaxis,predictorInfo.predictors_r(1:8,1:1e3));
set(gca,'ytick',yaxis,'xtick',[],'tickdir','out','color','m');
xlim([3 18]);

subplot(18,1,[9 16]);
imagesc(t_r(1:1e3), yaxis,predictorInfo.predictors_r(9:16,1:1e3));
set(gca,'ytick',yaxis,'xtick',[],'tickdir','out');
xlim([3 18]);

subplot(18,1,[17]);
imagesc(t_r(1:1e3), 1,predictorInfo.predictors_r(17,1:1e3));
set(gca,'xtick',[],'ytick',[]);
xlim([3 18]);

subplot(18,1,[18]);
imagesc(t_r(1:1e3), 1,predictorInfo.predictors_r(18,1:1e3));
set(gca,'tickdir','out','ytick',[],'xtick',[]);
hold on
line([16 17],[1 1],'linewidth',2)
xlabel('Time [s]');
xlim([3 18]);


colormap(1-gray);
