function dds = getCueSaccadeSmall(dd)
% dds = getCueSaccadeSmall(dd)

% % %showTonsetResp
% % dd.successTrials
% % dd.targetloc
% % 
% % %showTonsetByCue
% % dd.cOnset
% % dd.targetloc
% % 
% % %showFixCueOnsetResp
% % dd.cueOn


dds=[];
dds.started = dd.started;
dds.complete = dd.complete;%: [843×1 logical]
dds.successTrials = dd.successTrials;%: [843×1 double]
dds.fixationdt = dd.fixationdt;%: [843×1 double]
dds.theta = dd.theta;%: [843×1 double]
dds.thetaoff = dd.thetaoff;%: [843×1 double]
dds.dtheta=dd.dtheta;%: [843×1 double]
dds.pCued=dd.pCued;%: [843×1 double]
dds.cuedLoc=dd.cuedLoc;%: [843×1 double]
dds.cueOn=dd.cueOn;%: [843×1 logical]
dds.tOnset=dd.tOnset;%: [843×1 double]
dds.cOnset=dd.cOnset;%: [843×1 double]
dds.fixOnset=dd.fixOnset;%: [843×1 double]
dds.fOnset=dd.fOnset;%: [843×1 double]
dds.srt=dd.srt;%: [843×1 double]
dds.targetloc=dd.targetloc;%: [843×1 double]