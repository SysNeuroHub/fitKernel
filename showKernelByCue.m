function [fig, kernel_avg,prefdir, prefdirPval] = showKernelByCue(kernel_pop, ...
    tlags, cardinalDir, centering, tgtRange)
% created from showKernel3(kernel_pop, tlags, cardinalDir, centering)

prefDirOption = 0;%
%very slow in option=1
%something is not right in option=2

[ncol, nrow] = size(kernel_pop); %modality x units
nTrials = size(kernel_pop{1,1},2);
kernel_avg = [];prefdir=[];prefdirPval=[];
for icue = 1:2
    for col = 1:3;%5
        if centering %&& col<=3
            allmatrix = reshape(zeros(size(kernel_pop{col+5*(icue-1),1})),[],1);
            for row = 1:nrow
                allmatrix(:,row) = reshape(kernel_pop{col+5*(icue-1),row},1,[]);
            end
            orisize = size(kernel_pop{col+5*(icue-1),1});
            allmatrix = reshape(allmatrix, orisize(1), orisize(2),[]);

            tgtTimes = intersect(find(tlags{col+5*(icue-1)}(:,1)>tgtRange(col,1)), ...
                find(tlags{col+5*(icue-1)}(:,1)<tgtRange(col,2)));
            if icue == 1
                [directions, allmatrix_c, prefdir{col}, prefdirPval{col}]  = ...
                    alignMtxDir(allmatrix, tgtTimes, cardinalDir, prefDirOption );
            elseif icue == 2
               [directions, allmatrix_c, prefdir{col}, prefdirPval{col}]  = ...
                    alignMtxDir(allmatrix, tgtTimes, cardinalDir, prefDirOption, prefdir{col});
            end
            kernel_avg{col, icue} = squeeze(mean(allmatrix_c,3));
        else
            if ~centering
                directions = cardinalDir;
            end
            allmatrix = reshape(zeros(size(kernel_pop{col+5*(icue-1),1})),[],1);
            for row = 1:nrow
                allmatrix(:,row) = reshape(kernel_pop{col+5*(icue-1),row},1,[]);
            end
            kernel_avg{col, icue} = reshape(mean(allmatrix,2),size(kernel_pop{col+5*(icue-1),row}));
        end
    end
end

% created rom showKernel
fig = figure('position',[1          41        1440         783]);
for icue = 1:2
    a2(icue)=subplot(3,2,icue);
    thisIm = kernel_avg{1,icue}';
    crange = prctile(abs(thisIm(:)),99);
    %crange = prctile(thisIm(:),[1 99]);
    imagesc(tlags{1}(:,1), directions, thisIm);
    vline(0);
    if centering
        hline(0);
    end
    %caxis([-crange crange]);
    set(gca,'ytick',directions);
    xlabel('time from targetOnset [s]');
   colorbar;
    title(['n=' num2str(nTrials)]);

    a3(icue)=subplot(3,2,2+icue);
    thisIm = kernel_avg{2,icue}';
    crange = prctile(abs(thisIm(:)),99);
    %crange = prctile(thisIm(:),[1 99]);
    imagesc(tlags{2}(:,1),directions, thisIm);
    vline(0);
    if centering
        hline(0);
    end
    %caxis([-crange crange]);
    set(gca,'ytick',directions);
    xlabel('time from eye movement [s]');
   colorbar;

    a4(icue)=subplot(3,2,4+icue);
    thisIm = kernel_avg{3,icue}';
    crange = prctile(abs(thisIm(:)),99);
    %crange = prctile(thisIm(:),[1 99]);
    imagesc(tlags{3}(:,1),directions, thisIm);
    vline(0);
    if centering
        hline(0);
    end
    %caxis([-crange crange]);
    set(gca,'ytick',directions);
    xlabel('time from eye movement [s]');
   colorbar;
end
linkcaxes(a2,'mirror');%mcolorbar(a2(2),.5);
linkcaxes(a3,'mirror');%mcolorbar(a3(2),.5);
linkcaxes(a4,'mirror');%mcolorbar(a4(2),.5);
linkaxes([a2 a3 a4],'x');
