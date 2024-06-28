function [trialEndTimes, successTimes, trialOutcome] = getRewardTimes(dd)
%[trialEndTimes, successTimes, trialOutcome] = getRewardTimes(dd)
%returns times of reward delivery detected as item2delivered
%
%cf. concatenate_eye.m

eyeData = dd.eye;

t_cat = [];
nTrials = length(eyeData);
trialEndTimes = nan(1,nTrials);
successTimes = nan(1,nTrials);
trialOutcome = nan(1, nTrials);
for itr = 1:nTrials
    if isempty(eyeData(itr).t)
        continue;
    end

    if isempty(t_cat)
        t0 = eyeData(itr).t(1);
    else
        t0 =  max(t_cat)-eyeData(itr).t(1)+eyeData(itr).dt;
    end
    t_cat = cat(1, t_cat, eyeData(itr).t+t0);
    
    %end of trial to initiate reward delivery
    [time_s,trialInfo,~,state] = dd.meta.choice.state('trial',itr); % time stamps of the trial 'states'
    success_ind = strcmpi(state, 'SUCCESS');
    fail_ind = strcmp(state, 'FAIL');
        
    %actual delivery time of reward
    [time_d, trialInfo, frame, data] = dd.meta.newera.item2delivered('trial',itr);
    deliveryTime = time_d(end) - dd.meta.cic.firstFrame('trial',itr).time + t0; %<multiple entries to "time"
    if sum(success_ind)>0 %dd.successTrials(itr) 
        successTimes(itr) = time_s(success_ind) - dd.meta.cic.firstFrame('trial',itr).time + t0;
        trialEndTimes(itr) = deliveryTime;
        trialOutcome(itr) = 1;
    elseif sum(fail_ind)>0
        trialEndTimes(itr) = time_s(fail_ind) - dd.meta.cic.firstFrame('trial',itr).time + t0;
        trialOutcome(itr) = -1;
    end

end
%trialOutcome(~isnan(rewardTimes))=1;
%trialOutcome(~isnan(punishTimes))=-1;
% punishTimes = punishTimes(trialOutcome==-1);
% rewardTimes = rewardTimes(trialOutcome==1);
% successTimes = successTimes(~isnan(successTimes));


