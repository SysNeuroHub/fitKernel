function [choiceOutcome,succLoc, succDist] = getChoiceOutcome(dd, distTh)
% [choiceOutcome,succLoc, succDist] = getChoiceOutcome(dd, cardinalDir, distTh)
% choiceOutcome
% succLoc: location of the 1st succade after target onset [deg]
% succDist: distance of the 1st succade

if nargin < 2
    distTh = 3; %[deg]
    %Jo: always 5deg
end

[outcomes, cOnsetTimes] = getChoice(dd);

nTrials = numel(dd.started);
succLoc = nan(nTrials,1);
succDist = nan(nTrials,1);
for itrial = 1:nTrials
   
    %cOnsetTime = dd.cOnset(itrial); %only success trials
    cOnsetTime = cOnsetTimes(itrial);
    
    if ~isnan(cOnsetTime)
        [~,tidx] = min(abs(dd.eye(itrial).t - cOnsetTime));
        
        eyeRad = atan2(dd.eye(itrial).y(tidx), dd.eye(itrial).x(tidx)); %[-pi pi]
        %[~,minDirIdx] = min(abs(circ_dist(eyeRad, pi/180*cardinalDir)));
        succLoc(itrial) = eyeRad * 180/pi;%cardinalDir(minDirIdx); %direction of the saccade
        succDist(itrial) = sqrt(dd.eye(itrial).y(tidx).^2+dd.eye(itrial).x(tidx).^2); %saccade distance 
    end
end
     
successTr = find(dd.successTrials==1);
wrongSuccTr = find(dd.successTrials==0&~isnan(dd.tOnset)&succDist>=distTh);
quiescentTr = find(dd.successTrials==0&~isnan(dd.tOnset)&succDist<distTh);

choiceOutcome = nan(nTrials,1);
choiceOutcome(successTr) = 1; %succade to target location
choiceOutcome(wrongSuccTr) = 2; %succade to wrong location
choiceOutcome(quiescentTr) = 3; %did not succade

%% sanity check
% ax(1)=subplot(211)
% plot(dd.targetloc,'*');
% hold on
% plot(succLoc,'o');
% 
% successEvents=trace2Event(dd.successTrials==1);
% vbox(successEvents(:,1),successEvents(:,2),gca,[.5 1 .5]);
% 
% %saccade to wrong direction
% wrongSaccTr=trace2Event(dd.successTrials==0&outcomes==0&succDist>=distTh);
% vbox(wrongSaccTr(:,1),wrongSaccTr(:,2),gca,[1 .5 .5]);
% 
% %did not move eye
% quiescentTr=trace2Event(dd.successTrials==0&outcomes==0&succDist<distTh);
% vbox(quiescentTr(:,1), quiescentTr(:,2),gca,[.5 .5 .5]);
% 
% axis padded
% 
% ax(2)=subplot(212);
% plot(succDist,'sk');
% hline(distTh);
% linkaxes(ax(:),'x');
% 
