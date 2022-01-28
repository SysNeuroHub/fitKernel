function PSTH = getPSTH(spk, t)
%PSTH = getPSTH(spk, t)

if size(t,1) < size(t,2)
    t = t';
end

dt = median(diff(t));

tbins = cat(1, t, t(end)+dt) - 0.5*dt;
% hh = histogram(spk, tbins);
% PSTH = hh.Values'/dt; %better to do temporal smoothing??
N = histcounts(spk, tbins);
PSTH = N'/dt;