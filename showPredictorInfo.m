function fig = showPredictorInfo(predictorInfo, param, trange)
%fig = showPredictorInfo(predictorInfo, param.predictorNames, trange)

if nargin<3
    trange = [0 20]; %[s]
end

[~,trangeIdx(1)] = min(abs(predictorInfo.t_r-trange(1)));
[~,trangeIdx(2)] = min(abs(predictorInfo.t_r-trange(2)));
tidx = trangeIdx(1):trangeIdx(2);

twoD = find(predictorInfo.npredVars>1);

nWindows = numel(twoD) + 1;
ax = [];
fig = figure('position',[0 0 500 1000]);
for widx = 1:numel(twoD)
    ax(widx)=subplot(nWindows, 1, widx);
    if widx==1
        iidx = 1:predictorInfo.npredVars(widx);
    else
        iidx = 1+sum(predictorInfo.npredVars(1:widx-1)):sum(predictorInfo.npredVars(1:widx));
    end
    imagesc(predictorInfo.t_r(tidx), param.cardinalDir, predictorInfo.predictors_r(iidx,tidx))
    
    mcolorbar;
    ylabel(param.predictorNames{widx});
    set(gca,'tickdir','out');
end
colormap(1-gray);

ax(nWindows)=subplot(nWindows, 1, nWindows);
%yyaxis left;
plot(predictorInfo.t_r(tidx), predictorInfo.predictors_r(end-1,tidx));
ylabel(param.predictorNames{end-1});
% yyaxis right;
% plot(predictorInfo.t_r(tidx), predictorInfo.predictors_r(end,tidx));
% ylabel(param.predictorNames{end});
tt=trace2Event(predictorInfo.predictors_r(end,tidx));
vbox(predictorInfo.t_r(tidx(tt(:,1))), predictorInfo.t_r(tidx(tt(:,2))));
xlim([min(predictorInfo.t_r(tidx)) max(predictorInfo.t_r(tidx))]);

linkaxes(ax(:),'x');
