function [t_tr, catEvTimes, validTrials] = trimInvalids(t_tr, catEvTimes)
% [t_tr, catEvTimes, validTrials] = trimInvalids(t_tr, catEvTimes)
validTrials = find(cellfun(@numel,t_tr) > 1);
t_tr = t_tr(validTrials);
catEvTimes.tOnset = catEvTimes.tOnset(validTrials);
catEvTimes.fixOnset = catEvTimes.fixOnset(validTrials);
catEvTimes.cOnset = catEvTimes.cOnset(validTrials);
catEvTimes.cueOnset = catEvTimes.cueOnset(validTrials);
catEvTimes.fOnset = catEvTimes.fOnset(validTrials);

% CANNOT MODIFY DD
% dd.eye = dd.eye(validTrials);
% dd.numTrials = numel(validTrials);
% dd.started = dd.started(validTrials);
% dd.complete = dd.complete(validTrials);
% dd.successTrials = dd.successTrials(validTrials);
% dd.fixationdt = dd.fixationdt(validTrials);
% dd.theta = dd.theta(validTrials);
% dd.theataoff=dd.thetaoff(validTrials);
% dd.dtheta=dd.dtheta(validTrials);
% dd.pCued=dd.pCued(validTrials);
% dd.cuedLoc=dd.cuedLoc(validTrials);
% dd.cueOn=dd.cueOn(validTrials);
% dd.tOnset=dd.tOnset(validTrials);
% dd.cOnset=dd.cOnset(validTrials);
% dd.fixOnset=dd.fixOnset(validTrials);
% dd.fOnset=dd.fOnset(validTrials);
% dd.targetloc=dd.targetloc(validTrials);
