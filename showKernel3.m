function [fig, kernel_avg,prefdir, prefdirPval] = showKernel3(kernel_pop, ...
    tlags, cardinalDir, centering, tgtRange)
% fig = showKernel3(kernel_pop, tlags, cardinalDir, centering)

prefDirOption = 0;%
%very slow in option=1
%something is not right in option=2 

[ncol, nrow] = size(kernel_pop);
nUnits = nrow;%size(kernel_pop{1,1},2);
kernel_avg = [];prefdir=[];prefdirPval=[];
for col = 1:ncol %predictor type
    if centering && col<=3
        allmatrix = reshape(zeros(size(kernel_pop{col,1})),[],1);
        for row = 1:nrow %unit
            allmatrix(:,row) = reshape(kernel_pop{col,row},1,[]);
        end
        orisize = size(kernel_pop{col,1});
        allmatrix = reshape(allmatrix, orisize(1), orisize(2),[]);
        
        tgtTimes = intersect(find(tlags{col}(:,1)>tgtRange(col,1)), ...
            find(tlags{col}(:,1)<tgtRange(col,2)));
        [directions, allmatrix_c, prefdir{col}, prefdirPval{col}]  = ...
            alignMtxDir(allmatrix, tgtTimes, cardinalDir, prefDirOption );
        kernel_avg{col} = squeeze(mean(allmatrix_c,3));        
    else
        if ~centering
            directions = cardinalDir;
        end
        allmatrix = reshape(zeros(size(kernel_pop{col,1})),[],1);
        for row = 1:nrow
            allmatrix(:,row) = reshape(kernel_pop{col,row},1,[]);
        end
        kernel_avg{col} = reshape(mean(allmatrix,2),size(kernel_pop{col,row}));
    end
end

% created rom showKernel
fig = figure('position',[1          41        1440         783]);
a2=subplot(3,2,1);
thisIm = kernel_avg{1}';
crange = prctile(abs(thisIm(:)),99);
%crange = prctile(thisIm(:),[1 99]);
imagesc(tlags{1}(:,1), directions, thisIm);
vline(0);
if centering
    hline(0);
end
caxis([-crange crange]);
set(gca,'ytick',directions);
xlabel('time from targetOnset [s]');
mcolorbar(a2,.5);
title(['n=' num2str(nUnits)]);
set(gca,'tickdir','out');
   
a3=subplot(3,2,3);
thisIm = kernel_avg{2}';
crange = prctile(abs(thisIm(:)),99);
%crange = prctile(thisIm(:),[1 99]);
imagesc(tlags{2}(:,1),directions, thisIm);
vline(0);
if centering
    hline(0);
end
caxis([-crange crange]);
set(gca,'ytick',directions);
xlabel('time from eye movement [s]');
mcolorbar(a3,.5);
set(gca,'tickdir','out');

if size(kernel_avg{3},2)>1
    a4=subplot(3,2,5);
    thisIm = kernel_avg{3}';
    crange = prctile(abs(thisIm(:)),99);
    %crange = prctile(thisIm(:),[1 99]);
    imagesc(tlags{3}(:,1),directions, thisIm);
    vline(0);
    if centering
        hline(0);
    end
    caxis([-crange crange]);
    set(gca,'ytick',directions);
    set(gca,'tickdir','out');
    xlabel('time from eye movement [s]');
    mcolorbar(a4,.5);

    a5=subplot(3,2,2);
    plot(tlags{4}, [kernel_pop{4,:}],'color',[.7 .7 .7]);
    hold on
    plot(tlags{4}, kernel_avg{4}','k','linewidth',2);
    vline(0);
    xlabel('time from pupil dilation [s]');
    set(gca,'tickdir','out');
    axis tight;

    a6=subplot(3,2,4);
    plot(tlags{5}, [kernel_pop{5,:}],'color',[.7 .7 .7]);
    hold on
    plot(tlags{5}, kernel_avg{5}','k','linewidth',2);
    vline(0);
    xlabel('time from blink [s]');
    set(gca,'tickdir','out');
    axis tight;

    % a7=subplot(3,2,6);
    % plot(tlags{6}, [kernel_pop{6,:}],'color',[.7 .7 .7]);
    % hold on
    % plot(tlags{6}, kernel_avg{6}','k','linewidth',2);
    % vline(0);
    % xlabel('time from fixation onset [s]');
    % axis tight;

    linkaxes([a2 a3 a4 a5 a6],'x');
else
    %NOT YET
    a4=subplot(5,1,4);
    boundedline(tlags{3}, kernel_avg{3}', kernel_se{3}');hold on
    boundedline(tlags{4}, kernel_avg{4}', kernel_se{4}');hold on
    vline(0);
    xlabel('time from pupil dilation/blink [s]');
    axis tight;
    linkaxes([a2 a3 a4],'x');
end