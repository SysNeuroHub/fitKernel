

%% recorded data
animal = 'hugo';
rootFolder = '\\storage.erc.monash.edu.au\shares\R-MNHS-Physio\SysNeuroData\Monash Data\Joanita\2021/cuesaccade_data/';
saveFolder = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\Daisuke\cuesaccade_data';

excDur = 0.25; %[s]
respWin = [0 0.2]; %[s]

[loadNames, months, dates, channels] = getMonthDateCh(animal, rootFolder);

%idata = find(1-cellfun(@isempty, regexp(loadNames, regexptranslate('wildcard','06June\18\*_ch17.mat'))));

jdata = 1;
for idata = 1:length(channels) %1061;%865;%
datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
disp(datech);

saveSuffix = [animal replace(datech,'/','_')];

thisDate = [months{idata} '_' dates{idata}];

    
    saveName = fullfile(saveFolder, [saveSuffix '.mat']);
    
    if ~exist(saveName,'file')
        continue;
    end
    
    load(saveName, 'mFiringRate', 'param');%,'eyeData_rmotl_cat');
    
    %load(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');
    %load(fullfile(saveFolder,['eyeCat_' thisDate '.mat']), ...
    %    'onsets_cat','meta_cat','blinks','outliers','t_tr');
    
    load(saveName, 'mDir_sacc', 'mDir_sacc_pred', 'avgSaccResp', 'winSamps_sacc');

    mDir_sacc_pop(jdata) = mDir_sacc;
    mDir_sacc_pred_pop(jdata) = mDir_sacc_pred;
    avgSaccResp_pop(:,:,:,jdata) = avgSaccResp;%dir x variable x time
    jdata = jdata+1;
end

plot(mDir_sacc_pop, mDir_sacc_pred,'.');
marginplot;
squareplot;
xlabel('recorded');
ylabel('predicted by eyeposition kernel');
title('preferred direction of saccade response [rad]');
screen2png('saccResp_pop');
saveas(gcf,'saccResp_pop');



%% align to preferred direction
% cf fitPSTH_pop.m
centerBin = 4;
centeredDir = 180/pi*circ_dist(pi/180*param.cardinalDir, pi/180*param.cardinalDir(centerBin)); 
tgtTimes = intersect(find(winSamps_sacc>respWin(1)), find(winSamps_sacc<respWin(2)));
centered_avgSaccResp_pop = zeros(size(avgSaccResp_pop));
for jdata = 1:size(avgSaccResp_pop,4)
    thisData = squeeze(avgSaccResp_pop(:,:,:,jdata));
    [~,prefDir] = max(mean(thisData(:,1,tgtTimes),3));
    centered_avgSaccResp_pop(:,:,:,jdata) = circshift(thisData, centerBin - prefDir, 1);
end

nVars = size(avgSaccResp_pop,2);
for ivar = 1:nVars
   subplot(nVars,1,ivar); 
    imagesc(winSamps_sacc, centeredDir, squeeze(mean(centered_avgSaccResp_pop(:,ivar,:,:),4)));
    mcolorbar(gca,.5);
end
saveas(gcf,'saccResp_pop2');
screen2png('saccResp_pop2');