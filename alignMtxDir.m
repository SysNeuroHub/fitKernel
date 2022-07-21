function [centeredDir, centeredMtx] = alignMtxDir(mtx, tgtTidx, cardinalDir, prefDirs)
%[centeredDir, centeredMtx] = alignMtx(mtx, tgtTidx, cardinalDir)
% returns a matrix which 2nd dimension is rotated based on its activity
% during tgtTidx in the 1st dimension
%
% INPUTS:
%     mtx: [time x directions (x channels)]
%     tgtTidx: time index to compute preferred direction

if nargin < 4
    prefDirs = [];
end

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
for idata = 1:size(mtx,3)
    
    if ~isempty(prefDirs)
        prefDir = prefDirs(idata);
    else
        [~,prefDir] = max(mean(mtx(tgtTidx,:,idata),1));
    end
    centeredMtx(:,:,idata) = circshift(mtx(:,:,idata), centralBin - prefDir, 2);
end
if dims>=3
    centeredMtx = reshape(centeredMtx, [nTimes nDir nElements(3:end)]);
end