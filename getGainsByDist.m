function [gain_dir, gain_dist, winSamps_tonsetByCue, ...
    avgTonsetByCue] = getGainsByDist(t_r, ...
    y_r, cardinalDir, catEvTimes, dd, figTWin, onlySuccess, prefDir)

% y_r: time x 2 (observed or modelled)

% created from resp_cueConditions
onset = catEvTimes.tOnset;

cardinalDist = [0 45 90 135 180];
    

for icue = 1:2
    % icue==1 > without cue
    % icue==2 > with cue
    
    validEvents = intersect(find(~isnan(onset)), find(dd.cueOn==icue-1));
    %< this condition only includes all trials irrespective of the trial outcome

    if onlySuccess
        validEvents = intersect(validEvents, find(~isnan(dd.cOnset)));
    end
    
    
    onsetTimes = onset(validEvents);
    tgtDir = getTgtDir(dd.targetloc(validEvents), cardinalDir);
    
    %% tgt direction
    [~,dirIdx]=intersect(cardinalDir, unique(tgtDir));
    [avgTonsetByCue(dirIdx,:,:,icue), winSamps_tonsetByCue] ...
        = eventLockedAvg(y_r', t_r, onsetTimes, tgtDir, figTWin);
    
    %% distance to cue 
    distToCueDir = 180/pi*abs(circ_dist(pi/180*dd.targetloc(validEvents), pi/180*prefDir));
    [distToCueDir_q] = getTgtDir(distToCueDir, cardinalDist);
    [~,dirIdx2]=intersect(cardinalDist, unique(distToCueDir_q));
    [avgTonsetByCue_distCue(dirIdx2,:,:,icue), winSamps_tonsetByCue] ...
        = eventLockedAvg(y_r', t_r, onsetTimes, distToCueDir_q, figTWin);
     
end

%% get gain, by way of dividing observed signal with predicted by kernels
gain_dir = squeeze(avgTonsetByCue(:,1,:,:)./avgTonsetByCue(:,2,:,:));
gain_dist = squeeze(avgTonsetByCue_distCue(:,1,:,:)./avgTonsetByCue_distCue(:,2,:,:));

