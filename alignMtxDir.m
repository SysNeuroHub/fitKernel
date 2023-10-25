function [centeredDir, centeredMtx, prefDir, pvalue] = ...
    alignMtxDir(mtx, tgtTidx, cardinalDir, prefDirOption)
%[centeredDir, centeredMtx] = alignMtx(mtx, tgtTidx, cardinalDir)
% returns a matrix which 2nd dimension is rotated based on its activity
% during tgtTidx in the 1st dimension
%
% INPUTS:
%     mtx: [time x directions (x channels)]
%     tgtTidx: time index to compute preferred direction
% OUTPUTS:
%     prefDir: direction in [deg]

if nargin < 4
    prefDirOption = 1;
end

% if nargin < 4
%     prefDirs = [];
% end

centralBin = round(0.5*length(cardinalDir));
centralDir = cardinalDir(centralBin);

nTimes = size(mtx,1);
assert(length(tgtTidx) == length(intersect(1:nTimes, tgtTidx)));

nDir = size(mtx,2);
assert(nDir == length(cardinalDir));

nElements = size(mtx);
dims = length(nElements);

if dims>=3
    mtx = reshape(mtx, [nTimes nDir prod(nElements(3:end))]);
end

centeredDir = 180/pi*circ_dist(pi/180*cardinalDir, pi/180*centralDir);
%tgtTimes = intersect(find(kerneltlags>0.03), find(kerneltlags<0.25));

centeredMtx = zeros(size(mtx));
prefDir = nan(1,size(mtx,3));
pvalue = nan(1,size(mtx,3));
for idata = 1:size(mtx,3)
    resp = mean(mtx(tgtTidx,:,idata),1);
    switch prefDirOption
        case 0
            [~,prefBin] = max(resp);
            prefDir(idata) = cardinalDir(prefBin);
        case 1
            [fitPars, fitErr] ...
                = fitoriWrapped(cardinalDir, resp,...
                [], [nan nan 0 min(resp) nan],'',20, []);
            prefDir(idata) = fitPars(1);
            [~,prefBin] = min(abs(prefDir(idata) - cardinalDir));
        case 2
             prefDir(idata) = 180/pi*circ_mean(cardinalDir'*pi/180, resp');
             [~,prefBin] = min(abs(prefDir(idata) - cardinalDir));
    end
    centeredMtx(:,:,idata) = circshift(mtx(:,:,idata), centralBin - prefBin, 2);
    
    if nargout>=4
        pvalue(idata) =  dirDotProdTest(cardinalDir'*pi/180, resp');
    end
end
if dims>=3
    centeredMtx = reshape(centeredMtx, [nTimes nDir nElements(3:end)]);
end
