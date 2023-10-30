load('/mnt/syncitium/Daisuke/cuesaccade_data/figPSTH_pop20231026hugo/fitPSTH_pop20231026hugo.mat')


prefDirOption = 0;%

param.mfiringRateTh = 5;
param.expvalTh = 3;%0;
param.ntargetTrTh = 200;
param.ptonsetRespTh = 0.05;
param.ampTh = 0.5;%1; %for detecting preffered direction on 18/7/23

[okunits, mfiringRateOK, expvalOK, ntargetTrOK, ptonsetRespOK] ...
    = inclusionCriteria(mFiringRate_pop, expval_ind_pop(1,:), ntargetTrials_pop, PtonsetResp_pop, param);

%% hack exclude redundant data
[~, mfiringRate_u] = unique(mFiringRate_pop);
[~, expval_u] = unique(expval_ind_pop(1,:));
okunits_u = intersect(mfiringRate_u, expval_u);
okunits = intersect(okunits, okunits_u);

%% only retain ok data
gainInfo_pop = gainInfo_pop(okunits);
id_pop = id_pop(okunits);

thisID = find(strcmp(id_pop, 'hugo/2022/07July/26/19'));
[ncol, nrow] = size(gainInfo_pop);

winSamps = gainInfo_pop(1).winSamps;
cardinalDir = gainInfo_pop(1).cardinalDir;

centering = 0;
okgain = [];
gain_pop = [];
 avgTonsetByCue_pop = [];
for idata = 1:nrow
    gainInfo = gainInfo_pop(idata);

    avgTonsetByCue = squeeze(gainInfo.avgTonsetByCue(:,1,:,:));

    if sum(sum(avgTonsetByCue(:,:,2))) ~= 0
        okgain = [okgain idata];
    end

    if centering
        [~,prefBin] = min(abs(gainInfo.prefDir - gainInfo.cardinalDir));
        centralBin = round(0.5*length(gainInfo.cardinalDir));
        avgTonsetByCue = circshift(avgTonsetByCue, centralBin - prefBin, 1);
        centeredDir = 180/pi*circ_dist(pi/180*gainInfo.cardinalDir, pi/180*gainInfo.prefDir );
    end

    avgTonsetByCue(avgTonsetByCue==0) = nan;
    avgTonsetByCue_pop(:,:,:,idata) = avgTonsetByCue;

    gain = squeeze(avgTonsetByCue(:,:,2) ./ avgTonsetByCue(:,:,1));
    gain(gain==0) = nan;

        % gain = squeeze(log(avgTonsetByCue(:,1,:,2) ./ avgTonsetByCue(:,1,:,1)));
    % gain(isinf(gain))=nan;
    gain_pop(:,:,idata) = gain;

end

if centering
    yaxis = centeredDir;
else
    yaxis = cardinalDir;
end
ax(1)=subplot(211);
imagesc(winSamps, yaxis, squeeze(nanmean(avgTonsetByCue_pop(:,:,1,okgain(1:600)), 4)));
ax(2)=subplot(212);
imagesc(winSamps, yaxis, squeeze(nanmean(avgTonsetByCue_pop(:,:,2,okgain(1:600)), 4)));
linkcaxes(ax);





