function fig = showPredictorInfo(eyeData_rmotl_cat, predictorInfo, param, catEvTimes, dd, trange)
%fig = showPredictorInfo(predictorInfo, param.predictorNames, trange)

if nargin<3
    trange = [0 20]; %[s]
end

%% temporal window to show
[~,trangeIdx(1)] = min(abs(predictorInfo.t_r-trange(1)));
[~,trangeIdx(2)] = min(abs(predictorInfo.t_r-trange(2)));
tidx = trangeIdx(1):trangeIdx(2);

%% retrieve if trial was success or not
nTrials = dd.numTrials;
[rewardTimes_c, punishTimes_c, successTimes_c, trialOutcome] = getRewardTimes(dd);

rewardPunishTimes = nan(nTrials,1);
rewardPunishTimes(trialOutcome==1) = rewardTimes_c;%successTimes_c;%
rewardPunishTimes(trialOutcome==-1) = punishTimes_c;

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
for itrial = 1:numel(catEvTimes.fixOnset)
    [~, fixStartIdx(itrial)] = min(abs(predictorInfo.t_r - fixStart(itrial)));
    [~, trialStartIdx(itrial)] = min(abs(predictorInfo.t_r - trialStart(itrial)));
    [~, trialEndIdx(itrial)] = min(abs(predictorInfo.t_r - trialEnd(itrial)));
end
fixStart = predictorInfo.t_r(fixStartIdx);
trialStart = predictorInfo.t_r(trialStartIdx);
trialEnd = predictorInfo.t_r(trialEndIdx);


%% select only trials with target stimulus
fixStart = fixStart( ~isnan(catEvTimes.tOnset));
fixEnd = trialStart( ~isnan(catEvTimes.tOnset));

trialStart_reward = trialStart(trialOutcome==1);
trialEnd_reward = trialEnd(trialOutcome==1);

%'incorrect' trials: trialOutcome'==-1 & ~isnan(catEvTimes.cOnset) <> wrong choice was made
trialStart_nreward = trialStart(trialOutcome'==-1 & ~isnan(catEvTimes.tOnset)); %target stimulus was presented & unrewarded
trialEnd_nreward = trialEnd(trialOutcome'==-1 & ~isnan(catEvTimes.tOnset));


twoD = find(predictorInfo.npredVars>1);

nWindows = numel(twoD) + 2;
ax = [];
fig = figure('position',[0 0 500 1000]);

ax(1)=subplot(nWindows, 1, 1);
plot(eyeData_rmotl_cat.t, eyeData_rmotl_cat.x,...
    eyeData_rmotl_cat.t, eyeData_rmotl_cat.y);
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

ax(nWindows)=subplot(nWindows, 1, nWindows);
plot(predictorInfo.t_r(tidx), predictorInfo.predictors_r(end-1,tidx));
ylabel(param.predictorNames{end-1});

vbox(fixStart, fixEnd, gca, [.5 .5 1 .5]);
vbox(trialStart_reward, trialEnd_reward, gca, [.5 1 .5 .5]);
vbox(trialStart_nreward, trialEnd_nreward, gca, [1 .5 .5 .5]);

%% show blink times
tt=trace2Event(predictorInfo.predictors_r(end,tidx));
vbox(predictorInfo.t_r(tidx(tt(:,1))), predictorInfo.t_r(tidx(tt(:,2))));


xlim([min(predictorInfo.t_r(tidx)) max(predictorInfo.t_r(tidx))]);

linkaxes(ax(:),'x');
