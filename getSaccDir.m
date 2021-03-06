function [saccDir, dirIndex] = getSaccDir(startSacc, endSacc, eyeData_rmotl_cat, cardinalDir)
% [saccDir, dirIndex] = getSaccDir(startSacc, endSacc, eyeData_rmotl_cat, cardinalDir)
%
% TODO: make this function faster...

x = eyeData_rmotl_cat.x;
y = eyeData_rmotl_cat.y;
t = eyeData_rmotl_cat.t;


dirIndex = zeros(length(startSacc),1);
parfor isacc = 1:length(startSacc)
    tsnippet = intersect(find(t>=startSacc(isacc)), find(t<=endSacc(isacc)));
    %if isempty(tsnippet)
    %    continue;
    %end
    eyeRad = atan2(y(tsnippet)-y(tsnippet(1)), x(tsnippet)-x(tsnippet(1)));
    [~, minDirIdx] = arrayfun(@(x)(min(abs(circ_dist(x, pi/180*cardinalDir)))), eyeRad);
    dirIndex(isacc) = mode(minDirIdx);
end

saccDir = cardinalDir(dirIndex);