function dirMtx = getTgtDirMtx(dd, t_r, onsets_cat, cardinalDir)
% dirMtx = getTgtDirMtx(dd, t_r, tOnset_cat, cardinalDir)
%
% INPUT:
% dd
% t_r: time axis of concatenated events
% tOnset_cat: time of target onsets
% cardinalDir: list of directions of targets in deg
%
% OUTPUT:
% dirMtx: matrix [cardinalDir x t_r]

if nargin < 4
    cardinalDir = unique(dd.targetloc);
end

minDirIdx = zeros(dd.numTrials);
for itr = 1:dd.numTrials
    [~, minDirIdx(itr)] = min(abs(circ_dist(pi/180*dd.targetloc(itr), pi/180*cardinalDir)));
end

dirMtx = zeros(length(cardinalDir), length(t_r));
for itr = 1:dd.numTrials
    tOnset = onsets_cat.tOnset(itr);
    cOnset = onsets_cat.cOnset(itr);
    [~, onIdx] = min(abs(t_r - tOnset));
    [~, offIdx] = min(abs(t_r - cOnset));
    dirMtx(minDirIdx(itr), onIdx:offIdx) = 1;
end