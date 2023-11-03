%% pick up X number of functional units based on Rsqadj

nUnits = 40;

dataDir = '/mnt/syncitium/Daisuke/cuesaccade_data/figPSTH_pop20231026hugo/';

load(fullfile(dataDir, 'fitPSTH_pop20231026hugo.mat'), ...
    'Rsqadj_pop','mFiringRate_pop','expval_ind_pop','ntargetTrials_pop','PtonsetResp_pop',...
    'id_pop','gainInfo_pop');

load('/mnt/syncitium/Daisuke/cuesaccade_data/param20231102.mat', 'param');

%okunits = 1:numel(mFiringRate_pop);
[okunits, mfiringRateOK, expvalOK, ntargetTrOK, ptonsetRespOK] ...
    = inclusionCriteria(mFiringRate_pop, expval_ind_pop(1,:), ntargetTrials_pop, PtonsetResp_pop, ...
    param);

%% [HACK] exclude redundant data
[~, mfiringRate_u] = unique(mFiringRate_pop);
[~, expval_u] = unique(expval_ind_pop(1,:));
okunits_u = intersect(mfiringRate_u, expval_u);
okunits = intersect(okunits, okunits_u);

%% [HACK] exclude data without cue trials (from showGainInfo_pop)
%should include this in inclusionCriteria.m
okunits_gain = [];
for idata = 1:numel(gainInfo_pop)
    gainInfo = gainInfo_pop(idata);
    avgTonsetByCue = squeeze(gainInfo.avgTonsetByCue(:,1,:,:));

    if sum(sum(avgTonsetByCue(:,:,2))) ~= 0
        okunits_gain = [okunits_gain idata];
    end
end
okunits = intersect(okunits, okunits_gain);

%% retain ok units
id_pop = id_pop(okunits);
Rsqadj_pop = Rsqadj_pop(:,okunits);

Rsqadj_pop_a = Rsqadj_pop([2 3 4],:);
Rsqadj_pop_r = Rsqadj_pop([2 3 4],:)./Rsqadj_pop([1],:);

[~, idx_v_p] = sort(Rsqadj_pop_r(1,:));
[~, idx_es_p] = sort(Rsqadj_pop_r(2,:));
[~, idx_ep_p] = sort(Rsqadj_pop_r(3,:));
[~, idx_i_p] = sort(Rsqadj_pop_r(1,:)+Rsqadj_pop_r(2,:)+Rsqadj_pop_r(3,:));


idx_v_c = idx_v_p(1:nUnits);
idx_es_c = idx_es_p(1:nUnits);
idx_ep_c = idx_ep_p(1:nUnits);
idx_i_c = idx_i_p(1:nUnits);


idx_v = setdiff(idx_v_c, [idx_es_c idx_ep_c]);
idx_es = setdiff(idx_es_c, [idx_v_c idx_ep_c ]);
idx_ep = setdiff(idx_ep_c, [idx_v_c idx_es_c]);
idx_i = setdiff(idx_i_c, [idx_v idx_es idx_ep]);


funcClass.id_v = id_pop(idx_v);
funcClass.id_es = id_pop(idx_es);
funcClass.id_ep = id_pop(idx_ep);
funcClass.id_i = id_pop(idx_i);
funcClass.id_all = id_pop;

save(fullfile(dataDir,'pickUnitsByClass.mat'),"funcClass",'nUnits');


showScatterTriplets(1-Rsqadj_pop_r,{'vision','eye speed','eye position'},[-0.1 0.8],idx_v);
screen2png('pickedUnits_vision');

showScatterTriplets(1-Rsqadj_pop_r,{'vision','eye speed','eye position'},[-0.1 0.8],idx_es);
screen2png('pickedUnits_eyeSpeed');

showScatterTriplets(1-Rsqadj_pop_r,{'vision','eye speed','eye position'},[-0.1 0.8],idx_ep);
screen2png('pickedUnits_eyePosition');

showScatterTriplets(1-Rsqadj_pop_r,{'vision','eye speed','eye position'},[-0.1 0.8],idx_i);
screen2png('pickedUnits_integrator');

%id_pop([idx_v idx_es idx_ep])'

