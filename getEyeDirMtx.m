function dirMtx = getEyeDirMtx(dd, t_r, eyeData_rmblk_cat, cardinalDir)
% dirMtx = getTgtDirMtx(dd, t_r, tOnset_cat, cardinalDir)
%
% INPUT:
% dd
% t_r: time axis of concatenated events
% cardinalDir: list of directions of targets in deg
%
% OUTPUT:
% dirMtx: matrix [cardinalDir x t_r]

%TODO
% use cOnset as a proxy of saccade onsets?


if nargin < 4
    cardinalDir = unique(dd.targetloc);
end

x_r = interp1(eyeData_rmblk_cat.t, eyeData_rmblk_cat.x, t_r);
y_r = interp1(eyeData_rmblk_cat.t, eyeData_rmblk_cat.y, t_r);
eyeRad = atan2(y_r, x_r); %[-pi pi]
dist = sqrt(y_r.^2+x_r.^2);

minDirIdx = zeros(length(t_r),1);
for tt = 1:length(t_r)
    [~,minDirIdx(tt)] = min(abs(circ_dist(eyeRad(tt), pi/180*cardinalDir)));
end

dirMtx = zeros(length(cardinalDir), length(t_r));
for tt = 1:length(t_r)
    dirMtx(minDirIdx(tt),tt) = dist(tt);
end
