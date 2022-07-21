function f = showKernel(t_r, y_r, kernelInfo, cardinalDir)
%f = showKernel(t_r, y_r, kernelInfo, cardinalDir)

f = figure('position',[0 0 1000 500]);
subplot(1,2,1);
plot(t_r, y_r(:,1), 'color',[.5 .5 .5]);hold on
plot(t_r, y_r(:,2), 'linewidth',2);
%xlim([1935 1966]);
xlim([100 200])
legend('recorded','fitted');
xlabel('time [s]'); ylabel('firing rate [Hz]');

title(['expval: ' num2str(kernelInfo.expval), ', R: ' num2str(kernelInfo.corrcoef)]);

a2=subplot(3,2,2);
thisIm = kernelInfo.kernel{1}';
crange = prctile(abs(thisIm(:)),99);
%crange = prctile(thisIm(:),[1 99]);
imagesc(kernelInfo.tlags{1}(:,1), cardinalDir, thisIm);
caxis([-crange crange]);
set(gca,'ytick',cardinalDir);
xlabel('time from targetOnset [s]');
mcolorbar(a2,.5);

a3=subplot(3,2,4);
thisIm = kernelInfo.kernel{2}';
crange = prctile(abs(thisIm(:)),99);
%crange = prctile(thisIm(:),[1 99]);
imagesc(kernelInfo.tlags{2}(:,1),cardinalDir, thisIm);
caxis([-crange crange]);
set(gca,'ytick',cardinalDir);
xlabel('time from eye movement [s]');
mcolorbar(a3,.5);

a4=subplot(3,2,6);
plot(kernelInfo.tlags{3}, kernelInfo.kernel{3}');hold on
plot(kernelInfo.tlags{4}, kernelInfo.kernel{4}');hold on
xlabel('time from pupil dilation/blink [s]');
axis tight;