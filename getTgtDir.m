function [saccDir, dirIndex] = getTgtDir(targetloc, cardinalDir)
% [saccDir, dirIndex] = getTgtDir(targetloc, cardinalDir)
% 19/2/22 created from getSaccDir

[~, dirIndex] = arrayfun(@(x)(min(abs(circ_dist(x, pi/180*cardinalDir)))), pi/180*targetloc);
    
saccDir = cardinalDir(dirIndex);