%16/12/24 created from fitPSTH_pop.m

[saveServer, rootFolder] = getReady();

animal = 'ollie';
saveSuffix_p = ['20241212'];

n=load(fullfile(saveServer,'param20241206.mat'),'param');
param =n.param;
     
%% load result of mainscript_assemble.m
mFiringRate_pop = [];
expval_ind_pop = [];
ntargetTrials_pop = [];
PtonsetResp_pop = [];
id_pop = [];
kernel_pop = [];
expval_tgt_pop = [];
expval_tgt_rel_pop = [];
corr_tgt_pop = [];
corr_tgt_rel_pop = [];
latency_bhv_pop = [];
latency_neuro_pop = [];
latency_r_pop = [];
stats_stratified_pop = [];
p_hm_pop = [];
spkNGRate_pop = [];
CueTrRate_pop = [];
nLatencyTrials_pop = [];
nLatencyTrials_pref_pop = [];

for yy = 1:3
    switch yy
        case 1
            year = '2021';
        case 2
            year = '2022';
        case 3
            year = '2023';
    end
    saveFolder = fullfile(saveServer, year,animal);%17/6/23

    assemblyData = fullfile(saveFolder, ['assembly' saveSuffix_p '.mat']);

    if exist(assemblyData, 'file')
        assembly = load(assemblyData);
    else
        continue;
    end
    
    entries =  find(1-cellfun(@isempty, assembly.id_pop));

    %% inclusion criteria
    mFiringRate_pop = [mFiringRate_pop; [assembly.mFiringRate_pop{entries}]'];
    ntargetTrials_pop = [ntargetTrials_pop; [assembly.ntargetTrials_pop{entries}]'];
    PtonsetResp_pop = [PtonsetResp_pop; [assembly.PtonsetResp_pop{entries}]'];

    % unit stats
    id_pop = [id_pop; [assembly.id_pop(entries)]];
   
    kernel_pop = cat(2,kernel_pop, assembly.kernel_pop(:,entries));
    expval_tgt_pop = [expval_tgt_pop; assembly.expval_tgt_pop{entries}];
    corr_tgt_tmp = [assembly.corr_tgt_pop{entries}];
    corr_tgt_pop = [corr_tgt_pop; corr_tgt_tmp'];
    corr_tgt_rel_tmp = [assembly.corr_tgt_rel_pop{entries}]';
    corr_tgt_rel_pop = [corr_tgt_rel_pop; corr_tgt_rel_tmp];
   %avgAmp_hm_pop = [avgAmp_hm_pop; [assembly.avgAmp_hm_pop{entries}]']; unncessary???
   p_hm_pop = [p_hm_pop; [assembly.p_hm_pop(entries)]'];
   latency_r_pop = [latency_r_pop; [assembly.latency_r_pop(entries)]'];
   spkNGRate_pop = [spkNGRate_pop; [assembly.spkNGRate_pop{entries}]'];
   CueTrRate_pop = [CueTrRate_pop; [assembly.CueTrRate_pop{entries}]'];
   nLatencyTrials_pop = [nLatencyTrials_pop; assembly.nLatencyTrials_pop{entries}];
   nLatencyTrials_pref_pop = [nLatencyTrials_pref_pop; assembly.nLatencyTrials_pref_pop{entries}];

   %load just once
   tlags = assembly.tlags_pop{entries(1)};
end

%% apply inclusion critetia
% param.corr_tgtTh = 0.05;
% param.ptonsetRespTh=.1; 
% param.nLatencyTr_pref_th = 5;
[okunits, corr_tgtOK, ntargetTrOK, ptonsetRespOK] ...
    = inclusionCriteria(corr_tgt_pop(:,1), ntargetTrials_pop, PtonsetResp_pop, param);

    % omit data with the following criteria??
    % no saccade response
    % low spontaneous firing
    % low number of successful trials

%% retain included units
% units x 1
kernel_pop = kernel_pop(:,okunits); 
corr_tgt_pop = corr_tgt_pop(okunits,:);
corr_tgt_rel_pop = corr_tgt_rel_pop(okunits,:);
id_pop = id_pop(okunits);
mFiringRate_pop = mFiringRate_pop(okunits);
latency_r_pop = latency_r_pop(okunits);
p_hm_pop = p_hm_pop(okunits);
nLatencyTrials_pop = nLatencyTrials_pop(okunits);
nLatencyTrials_pref_pop = nLatencyTrials_pref_pop(okunits);
%stats_stratified_pop = stats_stratified_pop(okunits);
%stats_stratifiedByCue_pop = stats_stratifiedByCue_pop(okunits);
%expval_avgtgt_pop = expval_avgtgt_pop(:,okunits);
%corr_avgtgt_pop = corr_avgtgt_pop(:,okunits);
%Rsqadj_pop = Rsqadj_pop(okunits);


%% convert from cell to matrix
% thisTrialType = 1; %successful trials to preferred stimulus direction
% stratified_latency_pop_c = cellfun(@(a)a.latency(:,thisTrialType), stats_stratified_pop, 'UniformOutput', false);
% stratified_latency_pop = [stratified_latency_pop_c{:}];
% stratified_avgCorr_pop = cellfun(@(a)a.avgCorr(thisTrialType), stats_stratified_pop);
% stratified_pCorr_pop = cellfun(@(a)a.pCorr(thisTrialType), stats_stratified_pop);
 
latency_r_nb_pop = cellfun(@(a)a.r_success_pref, latency_r_pop);
latency_p_nb_pop = cellfun(@(a)a.p_success_pref, latency_r_pop);
p_hm_pop_before = cellfun(@(a)a(1), p_hm_pop);
p_hm_pop_after = cellfun(@(a)a(2), p_hm_pop);
p_hm_pop = [p_hm_pop_before p_hm_pop_after];


%% selected units
theseIDs = {'hugo/2021/08August/25/27',... %vision
    'hugo/2022/07July/26/19',... %eye speed
    'hugo/2022/08August/05/2'}; %integrator OK
[~, selectedIDs] = intersect(id_pop, theseIDs);

% high latency correlation, also eye-position driven
theseIDs2 = {'hugo/2022/07July/29/32'};
[~, selectedIDs2] = intersect(id_pop, theseIDs2);


%% kernel avg across units (FIG2)
%average kernel before centering
[f, kernel_avg] = showKernel3(kernel_pop, tlags, param.cardinalDir, 0);
savePaperFigure(f,fullfile(saveServer,saveSuffix_p,['avgKernel_' animal]));

%centerring by preferred direction
tgtRange = [0.05 0.15; 0.03 0.25; -0.1 0.1];
[f, kernel_centered_avg] = showKernel3(kernel_pop, tlags, param.cardinalDir, 1, tgtRange);
savePaperFigure(f,fullfile(saveServer,saveSuffix_p,['avgKernel_centered_' animal]));

%preferred direction across 3 kernels
[prefDir, amp] = getKernelPrefDirAmp(kernel_pop, tlags, tgtRange, param.cardinalDir);
%tuned = amp>param.ampTh;
tuned = 1:numel(okunits); %TEMP

figure('position',[ 680         485        1181         493]);
for ii = 1:3
    switch ii
        case 1
            v = [1 2];
        case 2
            v = [1 3];
        case 3
            v = [2 3];
    end
    doubleTuned = find(tuned(:,v(1))+tuned(:,v(2))==2);
    subplot(2,3,ii);
    plot(prefDir(:,v(1)), prefDir(:,v(2)), '.','color',[.7 .7 .7]);hold on
    plot(prefDir(doubleTuned,v(1)), prefDir(doubleTuned,v(2)), 'b.');
    squareplot;
    [rho, pval] = circ_corrcc(prefDir(doubleTuned,v(1))*pi/180, prefDir(doubleTuned,v(2))*pi/180);
    title(['rho:' num2str(rho) ', pval:' num2str(pval)]);
    xlabel(param.predictorNames{v(1)});
    ylabel(param.predictorNames{v(2)});
    subplot(2,3,ii+3)
    histogram(prefDir(:,v(1)) - prefDir(:,v(2)),-180:5:180, 'facecolor',[.7 .7 .7]); hold on;
    histogram(prefDir(doubleTuned,v(1)) - prefDir(doubleTuned,v(2)), -180:5:180,  'facecolor', 'b');
    vline(0);
    xlabel([param.predictorNames{v(1)} '-' param.predictorNames{v(2)}]);
    pval = circ_medtest(prefDir(doubleTuned,v(1)) - prefDir(doubleTuned,v(2)),0);
    title(['pval:' num2str(pval)]);
end
savePaperFigure(gcf,fullfile(saveServer,saveSuffix_p,['tuning_pop_' animal]));


%% cell type distibution (FIG3)
% 2023 JNS / Osaka / Hokkaido
% fig = showScatterTriplets(expval_tgt_pop(:,2:4)', ...
%     param.predictorNames, [-10 35], selectedIDs);
% screen2png(fullfile(saveServer,saveSuffix_p,['expval_tgt_a_' animal]));close;


% correlation during target presentation
fig = showScatterTriplets(corr_tgt_pop', ...
    param.predictorNames, [-.2 1], selectedIDs);
squareplots(fig, [-.2 1]);
savePaperFigure(fig, fullfile(saveServer,saveSuffix_p,['corr_tgt_' animal]));close;

fig_rel = showScatterTriplets(corr_tgt_rel_pop', ...
    param.predictorNames, [-100 100], selectedIDs);
squareplots(fig_rel, [-100 100]);
savePaperFigure(fig_rel, fullfile(saveServer,saveSuffix_p,['corr_tgt_r_' animal]));close;


%% latency stats on correlation to tgt (FIG4)
figure('position',[0 0 500 500]);
%latencyOK = (sum(isnan(stratified_latency_pop),1) < 1);
units_latency = nLatencyTrials_pref_pop >= param.nLatencyTr_pref_th;

%single-trial latency correlation 2024June
scatter(corr_tgt_rel_pop(:,2),corr_tgt_rel_pop(:,1),15,[.5 .5 .5])
hold on
scatter(corr_tgt_rel_pop(units_latency,2),corr_tgt_rel_pop(units_latency,1),15,latency_r_nb_pop(units_latency),'filled')
scatter(corr_tgt_rel_pop(selectedIDs,2),corr_tgt_rel_pop(selectedIDs,1),30,[1 1 0]);%stratified_avgCorr_pop(selectedIDs));
%scatter(corr_tgt_rel_pop(selectedIDs2,2),corr_tgt_rel_pop(selectedIDs2,1),30,[1 0 1]);%stratified_avgCorr_pop(selectedIDs));
clim([-0 .6]);
colormap("winter");colorbar;
xlabel('eye speed'); ylabel('stimulus');
title(sprintf('corr to tgt, relative to full mdl\nlatency correlation (success pref)'));
squareplot(gca, [-50 150]);
screen2png(fullfile(saveServer,saveSuffix_p,'corr_tgt_rel_pop_latency_r_success_pref'));close;


%% hit v miss (FIG 5)
units_hm = ~isnan(p_hm_pop(:,1));

figure('position',[0 0 1000 500]);
for iregress = 1:2
    subplot(1,2,iregress);
    scatter(corr_tgt_rel_pop(:,2),corr_tgt_rel_pop(:,1),15,[.5 .5 .5])
    hold on
    scatter(corr_tgt_rel_pop(units_hm,2),corr_tgt_rel_pop(units_hm,1),15,-log(p_hm_pop(units_hm, iregress)),'filled')
    clim([0 4]);
    colormap("cool");hh=colorbar; hh.Label.String = '-log(p hit v miss)';
    xlabel('eye speed'); ylabel('stimulus');
    tname = sprintf('Hit v Miss at preferred direction');
    if iregress == 2
        tname = [tname ' after regression'];
    end
    title(tname);
    squareplot(gca, [-50 150]);
end
screen2png(fullfile(saveServer,saveSuffix_p,'p_hm'));close;

% [val, id] = sort(p_hm_pop(units_hm,2), 'ascend');
% units_hmID = find(units_hm);
% selectedUnits = units_hmID(id(1:100));
%  id_pop(selectedUnits)

