function                 f = showKernel_prefDir( t_r, PSTH_f, dd, catEvTimes, kernelInfo, param, includeTrials)

prefDir = getPrefDir_wrapper(PSTH_f, t_r, dd, catEvTimes, param, includeTrials);
[~, prefDirIdx] = min(abs(prefDir - param.cardinalDir));

f = figure;
for ivar = 1:3
    plot(kernelInfo.tlags{ivar}(:,prefDirIdx), kernelInfo.kernel{ivar}(:, prefDirIdx));
    hold on;
end
set(gca,'tickdir','out');
legend(param.predictorNames(1:3));
xlabel('Time relative to covariate [s]');
title(['Pref dir: ' num2str(prefDir)])