saveServer = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\Daisuke\cuesaccade_data';%'E:/tmp/cuesaccade_data';
load('fitPSTH_pop20230704hugo.mat')
load(fullfile(saveServer,'param20230405.mat'),'param');
param.mfiringRateTh = 5;
param.expvalTh = 3;
param.ntargetTrTh = 200;
param.ptonsetRespTh = 0.05;
param.ampTh = 0.5;%1;

[okunits, mfiringRateOK, expvalOK, ntargetTrOK, ptonsetRespOK] ...
    = inclusionCriteria(mFiringRate_pop, expval_ind_pop(1,:), ntargetTrials_pop, PtonsetResp_pop, param);

kernel_pop = kernel_pop(:,okunits);
expval_tgt_pop = expval_tgt_pop(:,okunits);
id_pop = id_pop(:,okunits);

animal = 'hugo';
theseIDs = {'hugo/2021/09September/01/25',... %vision
    'hugo/2022/03March/10/20',... %eye speed
    'hugo/2022/07July/29/19'}; %eye position
[~, selectedIDs] = intersect(id_pop, theseIDs);



tgtRange = [0.05 0.15; 0.03 0.25; -0.1 0.1];

%fitPars = [];
prefDir = [];
amp = [];
for col = 1:3
    nrow = size(kernel_pop,2);
    allmatrix = reshape(zeros(size(kernel_pop{col,1})),[],1);
    for row = 1:nrow
        allmatrix(:,row) = reshape(kernel_pop{col,row},1,[]);
    end
    orisize = size(kernel_pop{col,1});
    allmatrix = reshape(allmatrix, orisize(1), orisize(2),[]);
    
    tgtTimes = intersect(find(tlags{col}(:,1)>tgtRange(col,1)), ...
        find(tlags{col}(:,1)<tgtRange(col,2)));
    
    for idata = 1:size(allmatrix,3)
        resp = mean(allmatrix(tgtTimes,:,idata),1)';
        prefDir(idata, col) = 180/pi*circ_mean(param.cardinalDir'*pi/180, resp);
        amp(idata, col) =  circ_r(param.cardinalDir'*pi/180, resp);
    end
end
tuned = amp>param.ampTh;


%% computation of preferred direction and significance in tuning
% tuning = abs(squeeze(fitPars(2,:,:)./fitPars(4,:,:)));%looks good!
% for col=1:3
%     tuned(:,col) = (tuning(:,col)>ampTh);
% end
% allTuned = find(sum(tuned,2)==3);
% 
% %% very few are tuned for all the three modalities(?)
% loglog(abs(tuning(:,1)), abs(tuning(:,2)), '.');
% 
% prefDir = squeeze(fitPars(1,:,:));


%% figures for explained variance
fig = showScatterTriplets(expval_tgt_pop(2:4,:), ...
    param.predictorNames, [], selectedIDs);
screen2png(['expval_tgt_a_' animal]);close;

fig = showScatterTriplets(100*expval_tgt_pop(2:4,:)./expval_tgt_pop(1,:), ...
    param.predictorNames, [-100 200], selectedIDs);
screen2png(['expval_tgt_r_' animal]);close;





