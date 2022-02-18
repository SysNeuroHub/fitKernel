function [rewardTimes, punishTimes, successTimes] = getRewardTimes(dd)
%[rewardTimes, successTimes] = getRewardTimes(dd)
%returns times of reward delivery detected as item2delivered
%
%cf. concatenate_eye.m

eyeData = dd.eye;

t_cat = [];
nTrials = length(eyeData);
rewardTimes = nan(1,nTrials);
punishTimes = nan(1,nTrials);
successTimes = nan(1,nTrials);
for itr = 1:nTrials
    
    if isempty(t_cat)
        t0 = eyeData(itr).t(1);
    else
        t0 =  max(t_cat)-eyeData(itr).t(1)+eyeData(itr).dt;
    end
    t_cat = cat(1, t_cat, eyeData(itr).t+t0);
    
    %end of trial to initiate reward delivery
    [time_s,trialInfo,~,state] = dd.meta.choice.state('trial',itr); % time stamps of the trial 'states'
    success_ind = strcmpi(state, 'SUCCESS');
        
    %actual delivery time of reward
    [time_d, trialInfo, frame, data] = dd.meta.newera.item2delivered('trial',itr);
    deliveryTime = time_d(end) - dd.meta.cic.firstFrame('trial',itr).time + t0; %<multiple entries to "time"
    if sum(success_ind)>0 %dd.successTrials(itr) 
        successTimes(itr) = time_s(success_ind) - dd.meta.cic.firstFrame('trial',itr).time + t0;
        rewardTimes(itr) = deliveryTime;
    else
        punishTimes(itr) = deliveryTime;
    end
end
punishTimes = punishTimes(~isnan(punishTimes));
rewardTimes = rewardTimes(~isnan(rewardTimes));
successTimes = successTimes(~isnan(successTimes));


