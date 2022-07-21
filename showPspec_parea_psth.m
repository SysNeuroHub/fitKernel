function [f] = showPspec_parea_psth(PSTH_f)
%[f] = showPspec_parea_psth(PSTH_f)
% NOT YET FUNCTIONAL
disp('spectral analysis');
f = figure('position',[0 0 1000 500]);
subplot(121);
loglog(faxis_parea, pspec_parea);
grid on
axis tight
xlabel('frequency [Hz]'); ylabel('psd');
title('parea');

[pspec_psth,faxis_psth] = pmtm(PSTH_f, 10, length(PSTH_f), 1/param.dt_r);%slow
subplot(122);
semilogy(faxis_psth, pspec_psth);
grid on
axis tight
xlabel('frequency [Hz]'); ylabel('psd');
title('psth');
screen2png(['pspec_parea_psth' saveSuffix]);
close;
