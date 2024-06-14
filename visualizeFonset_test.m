
[~, pdiam] = getPupilDiameter(eyeData_rmotl_cat);
pdiam_r = interp1(eyeData_rmotl_cat.t, pdiam, t_r)';
pdiam_rf= hpFilt(pdiam_r', 1/param.dt_r, .5)'%param.cutoffFreq)';
pdiam_rf= lpFilt(pdiam_rf', 1/param.dt_r, 2)'; %test 20 5
%pdiam_rf = [diff(pdiam_rf) 0]/param.dt_r; %spd

%% longitudinal
subplot(211);
plot(eyeData_rmotl_cat.tsample, zscore(eyeData_rmotl_cat.parea));hold on;
plot(t_r, zscore(pdiam_rf));
vbox(catEvTimes.blinkStartTimes, catEvTimes.blinkEndTimes, gca, 'k');
vbox(catEvTimes.outlierStartTimes, catEvTimes.outlierEndTimes);
vbox(catEvTimes.saccadeStartTimes, catEvTimes.saccadeEndTimes,gca,'g');
vline(catEvTimes.fOnset, gca, '-','r');
vline(catEvTimes.tOnset, gca, '-','b');
xlim([25 40]);

%% event trigerred avg
diamRange = [-3 3];%[-15 15];
[avgfOnsetResp, winSamps, singlefOnsetResp,~,uniqueLabels] ...
    = eventLockedAvg((pdiam_rf), t_r, catEvTimes.fOnset, [], [-5 15]);%[-0.5 1.5]);

subplot(2,4, [5:7]);
plot(winSamps, squeeze(singlefOnsetResp)); hold on
plot(winSamps, squeeze(avgfOnsetResp), 'k','linewidth',3);
%ylim(prctile(singlefOnsetResp(:),[1 99]))
ylim(diamRange);
vline(0, gca, '-','r');
hline(min(avgfOnsetResp))
set(gca,'TickDir','out');

subplot(2,4,8);
histogram(squeeze(pdiam_rf),'binlimits',diamRange, 'Orientation', 'horizontal');
hold on;
histogram(squeeze(singlefOnsetResp),'binlimits',diamRange, 'Orientation', 'horizontal');
ylim(diamRange);
