function fig = showPredictorInfo(eyeData_rmotl_cat, predictorInfo, param, catEvTimes, dd, ...
    PSTH_f, predicted_all,trange)
%fig = showPredictorInfo(predictorInfo, param.predictorNames, trange)

if nargin<8
    trange = [0 20]; %[s]
end

%% temporal window to show
[~,trangeIdx(1)] = min(abs(predictorInfo.t_r-trange(1)));
[~,trangeIdx(2)] = min(abs(predictorInfo.t_r-trange(2)));
tidx = trangeIdx(1):trangeIdx(2);

[~,trangeIdx_full(1)] = min(abs(eyeData_rmotl_cat.tsample - trange(1)));
[~,trangeIdx_full(2)] = min(abs(eyeData_rmotl_cat.tsample - trange(2)));
tidx_full = trangeIdx_full(1):trangeIdx_full(2);

%% retrieve if trial was success or not
nTrials = numel(dd.started);
[rewardPunishTimes, successTimes_c, trialOutcome] = getRewardTimes(dd);

fixStart = nan(nTrials,1);
trialStart = nan(nTrials,1);
trialEnd = nan(nTrials,1);
for itrial = 1:numel(catEvTimes.fixOnset)
    fixStart(itrial) = catEvTimes.fixOnset(itrial); %fixation cue onset
    trialStart(itrial) = catEvTimes.tOnset(itrial); %target stim onset
    trialEnd(itrial) = rewardPunishTimes(itrial);
end


%% downsampling 
fixStartIdx = []; trialStartIdx = []; trialEndIdx = [];
for itrial = 1:nTrials
    [~, fixStartIdx(itrial)] = min(abs(predictorInfo.t_r - fixStart(itrial)));
    [~, trialStartIdx(itrial)] = min(abs(predictorInfo.t_r - trialStart(itrial)));
    [~, trialEndIdx(itrial)] = min(abs(predictorInfo.t_r - trialEnd(itrial)));
end
fixStart = predictorInfo.t_r(fixStartIdx);
trialStart = predictorInfo.t_r(trialStartIdx);
trialEnd = predictorInfo.t_r(trialEndIdx);



%% select only trials with target stimulus
theseTrials = (catEvTimes.fixOnset > trange(1)).*(catEvTimes.fixOnset < trange(2));

fixStart = fixStart( ~isnan(catEvTimes.tOnset) & theseTrials);
fixEnd = trialStart( ~isnan(catEvTimes.tOnset) & theseTrials);

trialStart_reward = trialStart(trialOutcome'==1 & theseTrials);
trialEnd_reward = trialEnd(trialOutcome'==1 & theseTrials);

%'incorrect' trials: trialOutcome'==-1 & ~isnan(catEvTimes.cOnset) <> wrong choice was made
trialStart_nreward = trialStart(trialOutcome'==-1 & ~isnan(catEvTimes.tOnset)  & theseTrials); %target stimulus was presented & unrewarded
trialEnd_nreward = trialEnd(trialOutcome'==-1 & ~isnan(catEvTimes.tOnset) & theseTrials);


twoD = find(predictorInfo.npredVars>1);

nWindows = numel(twoD) + 3;
ax = [];
fig = figure('position',[0 0 500 1000]);

ax(1)=subplot(nWindows, 1, 1);
plot(eyeData_rmotl_cat.t(tidx_full), eyeData_rmotl_cat.x(tidx_full),...
    eyeData_rmotl_cat.t(tidx_full), eyeData_rmotl_cat.y(tidx_full));
set(gca,'tickdir','out');
xlim(trange);
hline(0);
vbox(fixStart, fixEnd, gca, [.5 .5 1 .5]);
vbox(trialStart_reward, trialEnd_reward, gca, [.5 1 .5 .5]);
vbox(trialStart_nreward, trialEnd_nreward, gca, [1 .5 .5 .5]);

for widx = 1:numel(twoD)
    ax(widx+1)=subplot(nWindows, 1, widx+1);
    if widx==1
        iidx = 1:predictorInfo.npredVars(widx);
    else
        iidx = 1+sum(predictorInfo.npredVars(1:widx-1)):sum(predictorInfo.npredVars(1:widx));
    end
    imagesc(predictorInfo.t_r(tidx), param.cardinalDir, predictorInfo.predictors_r(iidx,tidx),...
       'alphadata',predictorInfo.predictors_r(iidx,tidx));
    %drawRectanglesFromMatrix(predictorInfo.predictors_r(iidx,tidx),[],gcf,[nWindows,1,widx+1]);
    
    vbox(fixStart, fixEnd, gca, [.5 .5 1 .5]);
    vbox(trialStart_reward, trialEnd_reward, gca, [.5 1 .5 .5]);
    vbox(trialStart_nreward, trialEnd_nreward, gca, [1 .5 .5 .5]);

    mcolorbar;
    ylabel(param.predictorNames{widx});
    set(gca,'tickdir','out');
end
colormap(1-gray);

ax(nWindows-1)=subplot(nWindows, 1, nWindows-1);
plot(predictorInfo.t_r(tidx), predictorInfo.predictors_r(end-1,tidx),'k');
ylabel(param.predictorNames{end-1});

vbox(fixStart, fixEnd, gca, [.5 .5 1 .5]);
vbox(trialStart_reward, trialEnd_reward, gca, [.5 1 .5 .5]);
vbox(trialStart_nreward, trialEnd_nreward, gca, [1 .5 .5 .5]);

%% show blink times
tt=trace2Event(logical(predictorInfo.predictors_r(end,tidx)));
vbox(predictorInfo.t_r(tidx(tt(:,1))), predictorInfo.t_r(tidx(tt(:,2))));
    set(gca,'tickdir','out');

%% show spikes 
ax(nWindows)=subplot(nWindows, 1, nWindows);
plot(predictorInfo.t_r(tidx), PSTH_f(tidx), 'k');hold on; 
plot(predictorInfo.t_r(tidx), predicted_all(tidx), 'color',[.9 .3 .9], 'linewidth',2);
set(gca,'tickdir','out');

linkaxes(ax(:),'x');
xlim([min(predictorInfo.t_r(tidx)) max(predictorInfo.t_r(tidx))]);
