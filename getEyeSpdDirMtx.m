function dirMtx = getEyeSpdDirMtx(dd, t_r, eyeData_rmblk_cat, cardinalDir, excTimes)
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
    excTimes = [];
end
if nargin < 4
    cardinalDir = unique(dd.targetloc);
end

dt = median(diff(eyeData_rmblk_cat.t));
dt_r = median(diff(t_r));

order = 3;
framelen = 11;
x_f = sgolayfilt(eyeData_rmblk_cat.x, order, framelen);
y_f = sgolayfilt(eyeData_rmblk_cat.y, order, framelen);

x_r = interp1(eyeData_rmblk_cat.t, x_f, t_r);
y_r = interp1(eyeData_rmblk_cat.t, y_f, t_r);

blinks = (event2Trace(eyeData_rmblk_cat.t, excTimes)>0);
%blinks_r = (event2Trace(t_r, excTimes)>0);%NG
blinks_r = logical(interp1(eyeData_rmblk_cat.t, single(blinks), t_r, 'nearest'));

dx_r = [0; diff(x_r)]/dt_r;
dy_r = [0; diff(y_r)]/dt_r;

dx_r(blinks_r) = 0; 
dy_r(blinks_r) = 0; 

eyeRad_r = atan2(dy_r, dx_r); %[-pi pi]
dist_r = sqrt(dy_r.^2+dx_r.^2);

%test: without resampling
% blinks = (event2Trace(eyeData_rmblk_cat.t, excTimes)>0);
% dx = [0; diff(x_f)]/dt;
% dy = [0; diff(y_f)]/dt;
% 
% dx(blinks) = 0;
% dy(blinks) = 0;
% 
% eyeRad = atan2(dy, dx); %[-pi pi]
% dist = sqrt(dy.^2+dx.^2);
% 
% plot(eyeData_rmblk_cat.t, dist, t_r, dist_r);

[~, minDirIdx] = arrayfun(@(x)(min(abs(circ_dist(x, pi/180*cardinalDir)))), eyeRad_r);

dirMtx = zeros(length(cardinalDir), length(t_r));
for tt = 1:length(t_r)
    dirMtx(minDirIdx(tt),tt) = dist_r(tt);
end
