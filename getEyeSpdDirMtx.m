function dirMtx = getEyeSpdDirMtx(dd, t_r, eyeData_rmblk_cat, cardinalDir, blinks)
% dirMtx = getEyeSpdDirMtx(dd, t_r, tOnset_cat, cardinalDir)
%
% INPUT:
% dd
% t_r: time axis of concatenated events
% cardinalDir: list of directions of targets in deg
%
% OUTPUT:
% dirMtx: matrix [cardinalDir x t_r]

%WARNING
%this function will miss saccades if resampling rate is too low (>100Hz)

if nargin < 5
    blinks = zeros(length(eyeData_rmblk_cat.t),1);
end
if nargin < 4
    cardinalDir = unique(dd.targetloc);
end

dt = median(diff(eyeData_rmblk_cat.t));

order = 3;
framelen = 11;
x_f = sgolayfilt(eyeData_rmblk_cat.x, order, framelen);
y_f = sgolayfilt(eyeData_rmblk_cat.y, order, framelen);

dx_f = [0; diff(x_f)]/dt;
dy_f = [0; diff(y_f)]/dt;

dx_f(blinks) = 0;
dy_f(blinks) = 0;

%dx_r = interp1(eyeData_rmblk_cat.t, dx_f, t_r);
%dy_r = interp1(eyeData_rmblk_cat.t, dy_f, t_r);


eyeRad = atan2(dy_f, dx_f); %[-pi pi]
dist = sqrt(dy_f.^2+dx_f.^2);

eyeRad_r = interp1(eyeData_rmblk_cat.t, eyeRad, t_r);
dist_r = interp1(eyeData_rmblk_cat.t, dist, t_r);

[~, minDirIdx] = arrayfun(@(x)(min(abs(circ_dist(x, pi/180*cardinalDir)))), eyeRad_r);

dirMtx = zeros(length(cardinalDir), length(t_r));
for tt = 1:length(t_r)
    dirMtx(minDirIdx(tt),tt) = dist_r(tt);
end
