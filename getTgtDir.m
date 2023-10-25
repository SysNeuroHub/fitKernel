function [inputDir_q, dirIndex] = getTgtDir(inputDir, cardinalDir)%, targetDir)
% [saccDir, dirIndex] = getTgtDir(targetloc, cardinalDir)
% 19/2/22 created from getSaccDir
%
%INPUT:
% inputDir: original directions
% cardinalDir: quantized directions
% targetDir: if not empty, distance between inputDir and targetDir is registered
%
%OUTPUT:
% dirIndex: index of input directions in cardinalDir
% inputDir_q: quantized input direction

% inputDir = circ_dist(pi/180*inputDir, pi/180*targetDir);

[~, dirIndex] = arrayfun(@(x)(min(abs(circ_dist(x, pi/180*cardinalDir)))), pi/180*inputDir);
    
inputDir_q = cardinalDir(dirIndex);