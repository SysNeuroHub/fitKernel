function dirMtx = getDirMtx(dd, cardinalDir)

if nargin < 2
    cardinalDir = unique(dd.targetloc);
else
    cardinalDir = linspace(0, 360, 9);
    cardinalDir = cardinalDir(1:end-1);
end

minDirIdx = zeros(dd.numTrials);
for itr = 1:dd.numTrials
    [~, minDirIdx(itr)] = min(abs(circ_dist(pi/180*dd.targetloc(itr), pi/180*cardinalDir)));
end

dirMtx = zeros(length(cardinalDir), length(t_r));
for itr = 1:dd.numTrials
    tOnset = onsets_cat.tOnset(itr);
    [~,theseTimeIdx] = min(abs(t_r - tOnset));
    dirMtx(minDirIdx(itr), theseTimeIdx) = 1;
end