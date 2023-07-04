function fig = showKernel2(kernel_avg, kernel_se, tlags, cardinalDir)
%
% created rom showKernel
fig = figure('position',[0 0 500 1000]);
a2=subplot(4,1,1);
thisIm = kernel_avg{1}';
crange = prctile(abs(thisIm(:)),99);
%crange = prctile(thisIm(:),[1 99]);
imagesc(tlags{1}(:,1), cardinalDir, thisIm);
vline(0);
caxis([-crange crange]);
set(gca,'ytick',cardinalDir);
xlabel('time from targetOnset [s]');
mcolorbar(a2,.5);

a3=subplot(4,1,2);
thisIm = kernel_avg{2}';
crange = prctile(abs(thisIm(:)),99);
%crange = prctile(thisIm(:),[1 99]);
imagesc(tlags{2}(:,1),cardinalDir, thisIm);
vline(0);
caxis([-crange crange]);
set(gca,'ytick',cardinalDir);
xlabel('time from eye movement [s]');
mcolorbar(a3,.5);

if size(kernel_avg{3},2)>1
    a4=subplot(4,1,3);
    thisIm = kernel_avg{3}';
    crange = prctile(abs(thisIm(:)),99);
    %crange = prctile(thisIm(:),[1 99]);
    imagesc(tlags{3}(:,1),cardinalDir, thisIm);
    vline(0);
    caxis([-crange crange]);
    set(gca,'ytick',cardinalDir);
    xlabel('time from eye movement [s]');
    mcolorbar(a4,.5);

    a5=subplot(4,1,4);
    yyaxis left
    boundedline(tlags{4}, kernel_avg{4}', kernel_se{4}' , 'b');hold on
    yyaxis right
    boundedline(tlags{5}, kernel_avg{5}', kernel_se{5}', 'k');hold on
    vline(0);
    xlabel('time from pupil dilation/blink [s]');
    axis tight;
    linkaxes([a2 a3 a4 a5],'x');
else
    a4=subplot(4,1,4);
    boundedline(tlags{3}, kernel_avg{3}', kernel_se{3}');hold on
    boundedline(tlags{4}, kernel_avg{4}', kernel_se{4}');hold on
    vline(0);
    xlabel('time from pupil dilation/blink [s]');
    axis tight;
    linkaxes([a2 a3 a4],'x');
end