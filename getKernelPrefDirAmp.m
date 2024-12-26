function [prefDir, amp] = getKernelPrefDirAmp(kernel_pop, tlags, tgtRange, cardinalDir, prefDirOption)
%[prefDir, amp] = getKernelPrefDirAmp(kernel_pop, tlags, tgtRange)
%Input:
% kernel_pop: {modality x units}
%
%Output:
% prefDir: preferred direction of the kernel (units x modality) in [deg]
% amp: kernel amplitude at the preferred direction (units x modality)
%
% created from fitPSTH_pop.m & alignMtxDir
if nargin < 5
    prefDirOption = 2;
end

prefDir = [];
amp = [];
for col = 1:3 %modality
    nrow = size(kernel_pop,2);
    allmatrix = reshape(zeros(size(kernel_pop{col,1})),[],1);
    for row = 1:nrow
        allmatrix(:,row) = reshape(kernel_pop{col,row},1,[]);
    end
    orisize = size(kernel_pop{col,1});
    allmatrix = reshape(allmatrix, orisize(1), orisize(2),[]);

    tgtTimes = intersect(find(tlags{col}(:,1)>tgtRange(col,1)), ...
        find(tlags{col}(:,1)<tgtRange(col,2)));

    for idata = 1:size(allmatrix,3)
        resp = mean(allmatrix(tgtTimes,:,idata),1)';
         switch prefDirOption %from alignMtxDir
             case 0
                  [~,prefBin] = max(resp);
                prefDir(idata,col) = cardinalDir(prefBin);
             case 1
                 %NOT YET
             case 2
                 prefDir(idata, col) = 180/pi*circ_mean(cardinalDir'*pi/180, resp);
         end
        amp(idata, col) =  circ_r(cardinalDir'*pi/180, resp);
    end
end