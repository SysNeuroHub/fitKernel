function [respByCue, thisAxis] = dealRespByCue(respByCue, centering, showDistance, prefDir, cardinalDir)
%INPUT:
%respByCue: direction x time x wo/w cue
%
%OUTPUT:
%respByCue: 
%thisAxis: direction or distance from prefDir

%created from showGainInfo_pop.m
% gainInfo = gainInfo_pop(idata);
%     avgTonsetByCue = squeeze(gainInfo.avgTonsetByCue(:,1,:,:));

% if sum(sum(avgTonsetByCue(:,:,2))) ~= 0
%     okgain = [okgain idata];
% end

if centering
    centralBin = round(0.5*length(cardinalDir));
    centralDir = cardinalDir(centralBin);
    centeredDir = 180/pi*circ_dist(pi/180*cardinalDir, pi/180*centralDir );
    dirAxis = centeredDir;
    distAxis = unique(abs(round(dirAxis)));
else
    dirAxis = cardinalDir;
     centeredDir = 180/pi*circ_dist(pi/180*cardinalDir,0);
     distAxis = unique(abs(round(centeredDir)));
end


if centering %alignMtxDir
    [~,prefBin] = min(abs(prefDir - cardinalDir));
    centralBin = round(0.5*length(cardinalDir));
    respByCue = circshift(respByCue, centralBin - prefBin, 1);
end

respByCue(respByCue==0) = nan;

if showDistance
    avgTonsetByCue_c = [];
    for idist =1:numel(distAxis)
        theseDirs = find(abs(round(dirAxis)) == distAxis(idist));
        avgTonsetByCue_c(idist,:,:) = mean(respByCue(theseDirs,:,:),1);
    end
    respByCue = avgTonsetByCue_c;
end


if showDistance
    thisAxis = distAxis;
else
    thisAxis = dirAxis;
end