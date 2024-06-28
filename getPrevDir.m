function [prevDir_quantized, prevDirTrials] = getPrevDir(tgtDir, param)
%[prevDir_quantized, prevDirTrials] = getPrevDir(tgtDir, param)
% return the most prevalent stimulus direction
% Note this direction is irrelevant of neural response
% see also: getPrefDir

prevDir = mode(tgtDir);
[~, prevDirIdx] = min(abs(circ_dist(pi/180*param.cardinalDir, pi/180*prevDir)));
prevDir_quantized = param.cardinalDir(prevDirIdx);
prevDirTrials =  find(tgtDir == prevDir_quantized);%trials with cell's preferred direction
