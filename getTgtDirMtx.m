function dirMtx = getTgtDirMtx(dd, t_r, onsets_cat, cardinalDir, eyeData_cat)
% dirMtx = getTgtDirMtx(dd, t_r, tOnset_cat, cardinalDir)
%
% INPUT:
% dd
% t_r: time axis of concatenated events
% onsets_cat: time of target onsets
% eyeData_cat:
% cardinalDir: list of directions of targets in deg
%
% OUTPUT:
% dirMtx: matrix [cardinalDir x t_r]


tgtR = 5; %[deg] according to Mau
if nargin < 5
    eyeData_cat = [];
end

if nargin < 4
    cardinalDir = unique(dd.targetloc);
end

% [~,minDirIdx] = arrayfun(@(x)(min(abs(circ_dist(pi/180*x, pi/180*cardinalDir)))), dd.targetloc);

dirMtx = zeros(length(cardinalDir), length(t_r));
for itr = 1:dd.numTrials
    tOnset = onsets_cat.tOnset(itr);
    cOnset = onsets_cat.cOnset(itr);
    
    if isnan(tOnset) || isnan(cOnset)
        continue;
    end
        
    [~, onIdx] = min(abs(t_r - tOnset));
    [~, offIdx] = min(abs(t_r - cOnset));
    
    tgtX = tgtR*cos(pi/180*dd.targetloc(itr));
    tgtY = tgtR*sin(pi/180*dd.targetloc(itr));
    
    tidx = onIdx:offIdx;
    if ~isempty(eyeData_cat)
        eyeX = eyeData_cat.x(tidx);
        eyeY = eyeData_cat.y(tidx);
        
        tgtX_r = tgtX - eyeX;
        tgtY_r = tgtY - eyeY;
        targetloc_r = atan2(tgtY_r, tgtX_r); %[-pi pi]
    else
        targetloc_r = repmat(atan2(tgtY, tgtX), [numel(tidx),1]);
    end
    
    minDirIdx = [];
    for tt = 1:numel(tidx)
        [~,minDirIdx(tt)] = min(abs(circ_dist(targetloc_r(tt), pi/180*cardinalDir)));
        dirMtx(minDirIdx(tt), tidx(tt)) = 1;
    end
end

