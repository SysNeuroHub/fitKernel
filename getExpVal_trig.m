function [expval_trig, corr_trig, winSamps] = getExpVal_trig(PSTH_f, predicted, t_r, eventTimes, tWin)
%[expval_trig, corr_trig] = getExpVal_trig(PSTH_f, predicted, t_r, eventTimes, tWin)
% returns explained variance of triggered events, each moment in time from the event
%
% output:
% expval_trig: [kernel type x time]
% corr_trig: [kernel type x time]
% 13/6/24 created from getExpVal_tgt

%% event-triggered traces
[~, ~, resp] = eventLockedAvg(PSTH_f', t_r, eventTimes, [], tWin);
[~, winSamps, resp_predicted] = eventLockedAvg(predicted', t_r, eventTimes, [], tWin);

%% explained variance for each moment in time
expval_trig = zeros(size(predicted,2),numel(winSamps));
corr_trig = zeros(size(predicted,2),numel(winSamps));
for ivar = 1:size(predicted,2)
    for it = 1:numel(winSamps)
        [expval_trig(ivar, it), ~, corr_trig(ivar, it)] = getExpVal(resp(:,1,it), resp_predicted(:,ivar,it));
    end
end