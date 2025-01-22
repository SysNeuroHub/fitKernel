function f = showKernel(t_r, y_r, kernelInfo, cardinalDir)
%f = showKernel(t_r, y_r, kernelInfo, cardinalDir)
%
% 5/5/23 now compatible w both eye speed and position

useSameYrange = true;
if useSameYrange
    absMax = max(cellfun(@(x)abs(max(x(:))), kernelInfo.kernel));
end

f = figure('position',[0 0 1000 500]);
subplot(1,2,1);
plot(t_r, y_r(:,1), 'color',[.5 .5 .5]);hold on
plot(t_r, y_r(:,2), 'linewidth',2);
xlim([100 200])
legend('recorded','fitted');
xlabel('time [s]'); ylabel('firing rate [Hz]');

title(['expval: ' num2str(kernelInfo.expval), ', R: ' num2str(kernelInfo.corrcoef)]);

a2=subplot(4,2,2);
thisIm = kernelInfo.kernel{1}';
crange = prctile(abs(thisIm(:)),99);
if useSameYrange
    crange = absMax;
end
%crange = prctile(thisIm(:),[1 99]);
imagesc(kernelInfo.tlags{1}(:,1), cardinalDir, thisIm);
caxis([-crange crange]);
set(gca,'ytick',cardinalDir);
xlabel('time from targetOnset [s]');
mcolorbar(a2,.5);

a3=subplot(4,2,4);
thisIm = kernelInfo.kernel{2}';
crange = prctile(abs(thisIm(:)),99);
if useSameYrange
    crange = absMax;
end
imagesc(kernelInfo.tlags{2}(:,1),cardinalDir, thisIm);
caxis([-crange crange]);
set(gca,'ytick',cardinalDir);
xlabel('time from eye movement [s]');
mcolorbar(a3,.5);

if size(kernelInfo.kernel{3},2)>1
    a4=subplot(4,2,6);
    thisIm = kernelInfo.kernel{3}';
    crange = prctile(abs(thisIm(:)),99);
    if useSameYrange
        crange = absMax;
    end
    imagesc(kernelInfo.tlags{3}(:,1),cardinalDir, thisIm);
    caxis([-crange crange]);
    set(gca,'ytick',cardinalDir);
    xlabel('time from eye movement [s]');
    mcolorbar(a4,.5);

    a5=subplot(4,2,8);
    plot(kernelInfo.tlags{4}, kernelInfo.kernel{4}');hold on
    plot(kernelInfo.tlags{5}, kernelInfo.kernel{5}');hold on
    xlabel('time from pupil dilation/blink [s]');
    axis tight;
    if useSameYrange
        ylim([-absMax absMax]);
    end
    linkaxes([a2 a3 a4 a5],'x');
else
    a4=subplot(4,2,8);
    plot(kernelInfo.tlags{3}, kernelInfo.kernel{3}');hold on
    plot(kernelInfo.tlags{4}, kernelInfo.kernel{4}');hold on
    xlabel('time from pupil dilation/blink [s]');
    axis tight;
    linkaxes([a2 a3 a4],'x');
end