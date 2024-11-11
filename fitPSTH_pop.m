[saveServer, rootFolder] = getReady();

animal = 'hugo';
saveSuffix_p = ['20240715'];

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
%stats_stratifiedByCue_pop = [];
mampRespByCue_pop = [];
p_visRespByCue_pop = [];
p_cueModulationByCue_pop = [];
latency_r_cue_pop = [];
avgAmp_hm_pop = [];
p_hm_pop = [];
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
    load(fullfile(saveFolder, 'assembly.mat'),'assembly','param');
    tmp=load(fullfile(saveFolder, 'assembly_tmp.mat'),'assembly','param');
    assembly.avgAmp_hm_pop = tmp.assembly.avgAmp_hm_pop;
    assembly.p_hm_pop = tmp.assembly.p_hm_pop;
    
    save(fullfile(saveFolder, 'assembly.mat'),'assembly','param');
    
    entries =  find(1-cellfun(@isempty, assembly.id_pop));

    %% inclusion criteria
    mFiringRate_pop = [mFiringRate_pop; [assembly.mFiringRate_pop{entries}]'];
    expval_ind_pop = [expval_ind_pop; [assembly.expval_ind_pop{entries}]'];
    ntargetTrials_pop = [ntargetTrials_pop; [assembly.ntargetTrials_pop{entries}]'];
    PtonsetResp_pop = [PtonsetResp_pop; cellfun(@(a)a.PtonsetResp, assembly.cellclassInfo_pop(entries))];

    % unit stats
    id_pop = [id_pop; [assembly.id_pop(entries)]];
    kernel_pop = [kernel_pop; cellfun(@(a)a.kernel, assembly.kernelInfo_pop(entries), 'UniformOutput', false)];
    expval_tgt_pop = [expval_tgt_pop; assembly.expval_tgt_pop(entries)];
    expval_tgt_rel_pop = [expval_tgt_rel_pop; assembly.expval_tgt_rel_pop(entries)];
    corr_tgt_pop = [corr_tgt_pop; assembly.corr_tgt_pop(entries)];
    corr_tgt_rel_pop = [corr_tgt_rel_pop; assembly.corr_tgt_rel_pop(entries)];
   % hit v miss
   %avgAmp_hm_pop = [avgAmp_hm_pop; assembly.avgAmp_hm_pop()]
   p_hm_pop = [p_hm_pop; assembly.p_hm_pop(entries)];
   

    % latency 
    % latency_bhv_pop = [latency_bhv_pop; assembly.latency_bhv_pop{entries}'];
     %latency_neuro_pop = [latency_neuro_pop; [assembly.latency_neuro_pop{entries}]'];
     latency_r_pop = [latency_r_pop; assembly.latency_r_pop(entries)];
     stats_stratified_pop = [stats_stratified_pop; assembly.stats_stratified_pop(entries)];
     %stats_stratifiedByCue_pop = [stats_stratifiedByCue_pop; assembly.stats_stratifiedByCue_pop(entries)];
     mampRespByCue_pop = [mampRespByCue_pop; assembly.mampRespByCue_pop(entries)];
     p_visRespByCue_pop = [p_visRespByCue_pop; assembly.p_visRespByCue_pop(entries)];
     p_cueModulationByCue_pop = [p_cueModulationByCue_pop; assembly.p_cueModulationByCue_pop(entries)];
     latency_r_cue_pop = [latency_r_cue_pop; assembly.latency_r_cue_pop(entries)];
end

%% apply inclusion critetia
% [okunits, mfiringRateOK, expvalOK, ntargetTrOK, ptonsetRespOK] ...
%     = inclusionCriteria(mFiringRate_pop, expval_pop, ntargetTrials_pop, PtonsetResp_pop, param);
% [okunits, mfiringRateOK, expvalOK, ntargetTrOK, ptonsetRespOK] ...
%     = inclusionCriteria(mFiringRate_pop, expval_tgt_pop(1,:), ntargetTrials_pop, PtonsetResp_pop, param);
[okunits, mfiringRateOK, expvalOK, ntargetTrOK, ptonsetRespOK] ...
    = inclusionCriteria(mFiringRate_pop, expval_ind_pop(:,1), ntargetTrials_pop, PtonsetResp_pop, param);

    % omit data with the following criteria??
    % no saccade response
    % low spontaneous firing
    % low number of successful trials

%% retain included units
% units x 1
kernel_pop = kernel_pop(okunits); %FIXME
expval_ind_pop = expval_ind_pop(okunits,:);
expval_tgt_pop = [expval_tgt_pop{okunits}]';
expval_tgt_rel_pop = [expval_tgt_rel_pop{okunits}]';
corr_tgt_pop = [corr_tgt_pop{okunits}]';
corr_tgt_rel_pop = [corr_tgt_rel_pop{okunits}]';
id_pop = id_pop(okunits);
mFiringRate_pop = mFiringRate_pop(okunits);
latency_r_pop = latency_r_pop(okunits);
stats_stratified_pop = stats_stratified_pop(okunits);
%stats_stratifiedByCue_pop = stats_stratifiedByCue_pop(okunits);
%expval_avgtgt_pop = expval_avgtgt_pop(:,okunits);
%corr_avgtgt_pop = corr_avgtgt_pop(:,okunits);
%Rsqadj_pop = Rsqadj_pop(okunits);
latency_r_cue_pop = latency_r_cue_pop(okunits);
mampRespByCue_pop = mampRespByCue_pop(okunits);
p_visRespByCue_pop = p_visRespByCue_pop(okunits);
p_cueModulationByCue_pop = p_cueModulationByCue_pop(okunits);
p_hm_pop = p_hm_pop(okunits);


%% convert from cell to matrix
thisTrialType = 1; %successful trials to preferred stimulus direction
stratified_latency_pop_c = cellfun(@(a)a.latency(:,thisTrialType), stats_stratified_pop, 'UniformOutput', false);
stratified_latency_pop = [stratified_latency_pop_c{:}];
stratified_avgCorr_pop = cellfun(@(a)a.avgCorr(thisTrialType), stats_stratified_pop);
stratified_pCorr_pop = cellfun(@(a)a.pCorr(thisTrialType), stats_stratified_pop);

latency_r_nb_pop = cellfun(@(a)a.r_success_pref, latency_r_pop);
latency_p_nb_pop = cellfun(@(a)a.p_success_pref, latency_r_pop);

latency_r_cue_pop = [latency_r_cue_pop{:}];
mampRespByCue_pop = [mampRespByCue_pop{:}];
p_visRespByCue_pop = [p_visRespByCue_pop{:}];
p_cueModulationByCue_pop = [p_cueModulationByCue_pop{:}];
p_hm_pop_before = cellfun(@(a)a(1), p_hm_pop);
p_hm_pop_after = cellfun(@(a)a(2), p_hm_pop);
p_hm_pop = [p_hm_pop_before p_hm_pop_after];


%% time-resolved exp val
%expval_trig_pop = expval_trig_pop(:,:,:,okunits);
%mexpval_trig_pop = squeeze(mean(expval_trig_pop,4));
%
% for ievtype = 1:4
%     ax_expval_trig(ievtype) = subplot(4,1,ievtype);
%     plot(winSamps, squeeze(mexpval_trig_pop(:,:,ievtype))');
% end
% linkaxes(ax_expval_trig);
% ylim([-10 40])
% legend(['all',param.predictorNames],'location','west');
% screen2png(fullfile(saveServer,[saveSuffix_p animal '_expval_trig' limSuffix '.png']));


%% selected units
% theseIDs = {'hugo/2021/09September/01/25',...
%     'hugo/2021/11November/16/6',...
%     'hugo/2021/12December/14/13'};
% theseIDs = {'hugo/2021/09September/01/25',...
%     'hugo/2021/11November/02/18',...
%     'hugo/2022/07July/29/19'};
% theseIDs = {'hugo/2022/03March/10/20',... %eye speed driven
%     'hugo/2022/07July/29/19'}; %eye position driven
theseIDs = {'hugo/2021/08August/25/27',... %vision
    'hugo/2022/07July/26/19',... %eye speed
    'hugo/2022/08August/05/2'} %integrator OK
%    'hugo/2022/08August/15/5'}
%     'hugo/2021/12December/13/13'}; %integrator ... too law exp var
%    'hugo/2021/12December/09/8'}; %integrator ... too law exp var
%    'hugo/2021/12December/14/13'};%integrator ... too law exp var
%     'hugo/2021/09September/01/25',... %NG did not meet incl critaria

[~, selectedIDs] = intersect(id_pop, theseIDs);

% high latency correlation, also eye-position driven
theseIDs2 = {'hugo/2022/07July/29/32'};

[~, selectedIDs2] = intersect(id_pop, theseIDs2);

% %% histogram of individual explained variance
% figure('position',[ 1120         454         787         500]);
% for ii = 1:size(expval_ind_pop,1)
%     ax(ii) = subplot(size(expval_ind_pop,1),2,2*ii-1);
%     histogram(expval_ind_pop(ii,:),[-10:1:20]);
%     if ii==1
%         ylabel('all');
%         title('expval ind pop');
%     else
%         ylabel(param.predictorNames{ii-1});
%     end
% 
%     ax(ii) = subplot(size(expval_ind_pop,1),2,2*ii);
%     histogram(expval_tgt_pop(ii,:),[-10:1:30]);
%     if ii==1
%         title('expval tgt pop');
%     end
% end
% xlabel('Explained variance [%]')
% savePaperFigure(gcf,['expval_' animal]);
% 
% %% scatter plot of Rsquare adjusted
% % fig = showScatterTriplets(Rsqadj_pop([2 1 3],:), ...
% %     {'wo eye speed','full mdl','wo eye pos'}, [0 .5]);
% fig = showScatterTriplets(Rsqadj_pop([2 3 4],:), ...
%     {'wo vision','wo eye spd','wo eye pos'}, [0 .5]);
% squareplots;
% savePaperFigure(gcf,fullfile(saveServer,saveSuffix_p,['Rsqadj_' animal]));
% 
% fig = showScatterTriplets(Rsqadj_pop([2 3 4],:)./Rsqadj_pop([1],:), ...
%     {'wo vision','wo eye spd','wo eye pos'});
% squareplots;
% savePaperFigure(gcf,fullfile(saveServer,saveSuffix_p,['Rsqadj_' animal '_relativeToFullMdl']));
% 
% %% scatter plot of individual explained variances
% fig = showScatterTriplets(expval_ind_pop(2:4,:), ...
%     param.predictorNames, [-4 20], selectedIDs);
% screen2png(fullfile(saveServer,saveSuffix_p,['expval_all_a' animal]));
% 
% fig = showScatterTriplets(100*expval_ind_pop(2:4,:)./expval_ind_pop(1,:), ...
%     param.predictorNames, [-100 200], selectedIDs);
% screen2png(fullfile(saveServer,saveSuffix_p,['expval_all_r' animal]));

%% explained variance between kernels
% 2023 JNS / Osaka / Hokkaido
% fig = showScatterTriplets(expval_tgt_pop(:,2:4)', ...
%     param.predictorNames, [-10 35], selectedIDs);
% screen2png(fullfile(saveServer,saveSuffix_p,['expval_tgt_a_' animal]));close;

%expval_tgt_rel = 100*expval_tgt_pop(:,2:4)./expval_tgt_pop(:,1);
fig = showScatterTriplets(expval_tgt_rel_pop', ...
    param.predictorNames, [-100 200], selectedIDs);
squareplots
savePaperFigure(gcf,fullfile(saveServer,saveSuffix_p,['expval_tgt_r_' animal]));close;

%% correlation during target presentation
fig = showScatterTriplets(corr_tgt_pop', ...
    param.predictorNames, [-.2 1], selectedIDs);
squareplots
savePaperFigure(gcf,fullfile(saveServer,saveSuffix_p,['corr_tgt_' animal]));close;

fig = showScatterTriplets(corr_tgt_rel_pop', ...
    param.predictorNames, [-100 200], selectedIDs);
squareplots
savePaperFigure(gcf,fullfile(saveServer,saveSuffix_p,['corr_tgt_r_' animal]));close;

%% hit v miss
units_hm = ~isnan(p_hm_pop(:,1));

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
screen2png(fullfile(saveServer,saveSuffix_p,'p_hm'));

[val, id] = sort(p_hm_pop(units_hm,2), 'ascend');
units_hmID = find(units_hm);
selectedUnits = units_hmID(id(1:100));
 id_pop(selectedUnits)

%% latency stats on correlation to tgt
% avg-trial latency correlation
latencyOK = (sum(isnan(stratified_latency_pop),1) < 1);
% scatter(expval_tgt_rel(~latencyOK,2),expval_tgt_rel(~latencyOK,1),15,[.5 .5 .5])
scatter(corr_tgt_rel_pop(:,2),corr_tgt_rel_pop(:,1),15,[.5 .5 .5])
hold on
scatter(corr_tgt_rel_pop(latencyOK,2),corr_tgt_rel_pop(latencyOK,1),15,stratified_avgCorr_pop(latencyOK),'filled')
scatter(corr_tgt_rel_pop(selectedIDs,2),corr_tgt_rel_pop(selectedIDs,1),30,[1 1 0]);%stratified_avgCorr_pop(selectedIDs));
%scatter(corr_tgt_rel_pop(selectedIDs2,2),corr_tgt_rel_pop(selectedIDs2,1),30,[1 0 1]);%stratified_avgCorr_pop(selectedIDs));
colormap("winter");colorbar;
xlabel('eye speed'); ylabel('stimulus');
title('expval to tgt, relative to full mdl');
squareplot(gca, [-50 150]);
screen2png(fullfile(saveServer,saveSuffix_p,'corr_tgt_rel_pop_stratified_avgCorr_valid7bins'));

%single-trial latency correlation 2024June
scatter(corr_tgt_rel_pop(:,2),corr_tgt_rel_pop(:,1),15,[.5 .5 .5])
hold on
scatter(corr_tgt_rel_pop(latencyOK,2),corr_tgt_rel_pop(latencyOK,1),15,latency_r_nb_pop(latencyOK),'filled')
scatter(corr_tgt_rel_pop(selectedIDs,2),corr_tgt_rel_pop(selectedIDs,1),30,[1 1 0]);%stratified_avgCorr_pop(selectedIDs));
%scatter(corr_tgt_rel_pop(selectedIDs2,2),corr_tgt_rel_pop(selectedIDs2,1),30,[1 0 1]);%stratified_avgCorr_pop(selectedIDs));
clim([-0 .6]);
colormap("winter");colorbar;
xlabel('eye speed'); ylabel('stimulus');
title(sprintf('corr to tgt, relative to full mdl\nlatency correlation (success pref)'));
squareplot(gca, [-50 150]);
screen2png(fullfile(saveServer,saveSuffix_p,'corr_tgt_rel_pop_latency_r_success_pref'));

%% neuro-bhv correlation v cue-bhv correlation
r_cue = [latency_r_cue_pop.r_success_prev]';
p_cue = [latency_r_cue_pop.p_success_prev]';

pOK = logical((p_cue<0.05) .* (latency_p_nb_pop<0.05));

scatter(latency_r_nb_pop(latencyOK), r_cue(latencyOK));
hold on
scatter(latency_r_nb_pop(pOK), r_cue(pOK), 'filled');

xlabel('neuro-bhv latency correlation');
ylabel('cue-bhv latency correlation');
squareplot(gca, [-0.4 1]);
screen2png(fullfile(saveServer,saveSuffix_p,'neuro-bhv-cue correlation'));


%% latency correlation
scatter(stratified_avgCorr_pop, latency_r_pop, 15); hold on;
scatter(stratified_avgCorr_pop(latencyOK), latency_r_pop(latencyOK),'filled');
xlabel('stratified avgCorr'); ylabel('latency corr (success pref)');
squareplot;
screen2png('latencyCorr_v_stratified_avgCorr_valid7bins');

visionUnits = logical(latencyOK'.*(expval_tgt_rel(:,1)>expval_tgt_rel(:,2))); %TMP
eyespdUnits = logical(latencyOK'.*(expval_tgt_rel(:,1)<expval_tgt_rel(:,2))); %TMP
hh(1)=histogram(latency_r_pop(eyespdUnits),-1:.1:1); hold on;
hh(2)=histogram(latency_r_pop(visionUnits),-1:.1:1); 
%[~,p_vision] = ttest(latency_r_success_pref_pop(visionUnits));
%[~,p_eyespd] = ttest(latency_r_success_pref_pop(eyespdUnits));
[~, p] = ttest2(latency_r_pop(visionUnits), latency_r_pop(eyespdUnits));
sigstar([nanmean(latency_r_pop(visionUnits)) nanmean(latency_r_pop(eyespdUnits))],p)
xlabel('single-trial correlation'); ylabel('#units');
legend(hh(:),'eye speed','stimulus');
screen2png('hist_latencyCorr_valid7bins');


%% effect of cue
respUnits = (p_visRespByCue_pop< 0.05);
scatter(corr_tgt_rel_pop(:,2),corr_tgt_rel_pop(:,1),15,[.5 .5 .5])
hold on
scatter(corr_tgt_rel_pop(respUnits,2),corr_tgt_rel_pop(respUnits,1),15,-log(p_cueModulationByCue_pop(respUnits)),'filled')
clim([0 4]);
colormap("cool");hh=colorbar; hh.Label.String = '-log(p cueModulation)';
xlabel('eye speed'); ylabel('stimulus');
title(sprintf(' wo cue v w cue to tgt at 0deg'));
squareplot(gca, [-50 150]);
screen2png(fullfile(saveServer,saveSuffix_p,'cueEffect_pval'));

respUnits = (p_visRespByCue_pop< 0.05);
scatter(corr_tgt_rel_pop(:,2),corr_tgt_rel_pop(:,1),15,[.5 .5 .5])
hold on
metric = diff(mampRespByCue_pop(:,respUnits))./mean(mampRespByCue_pop(:,respUnits));
scatter(corr_tgt_rel_pop(respUnits,2),corr_tgt_rel_pop(respUnits,1),15, metric,'filled')
clim([-0.5 0.5]);
colormap("jet");hh=colorbar; hh.Label.String = 'woCue - wCue';
xlabel('eye speed'); ylabel('stimulus');
title(sprintf(' wo cue v w cue to tgt at 0deg'));
squareplot(gca, [-50 150]);
screen2png(fullfile(saveServer,saveSuffix_p,'cueEffect_frate_r'));

% %%
 [val, id] = sort(p_cueModulationByCue_pop(respUnits), 'ascend');
 respUnitsID = find(respUnits);
 selectedUnits = respUnitsID(id(1:30));
 scatter(corr_tgt_rel_pop(selectedUnits,2),corr_tgt_rel_pop(selectedUnits,1),15, ...
     -log(p_cueModulationByCue_pop(selectedUnits)),'filled');
hold on;
text(corr_tgt_rel_pop(selectedUnits,2),corr_tgt_rel_pop(selectedUnits,1), id_pop(selectedUnits), 'fontsize',6);
clim([0 5]);
xlabel('eye speed'); ylabel('stimulus');
title(sprintf(' wo cue v w cue (<0.6s delay) v w cue (>0.6s delay) to tgt at 0deg'));
squareplot(gca, [-50 150]);
screen2png(fullfile(saveServer,saveSuffix_p,'cueEffect_pval_top30'));

%% correlation on avg response between kernels
% fig = showScatterTriplets(corr_tgt_pop(2:4,:), ...
%     param.predictorNames, [], selectedIDs);
% screen2png(['corr_tgt_a_' animal]);close;
%
% fig = showScatterTriplets(100*corr_tgt_pop(2:4,:)./corr_tgt_pop(1,:), ...
%     param.predictorNames, [-20 100], selectedIDs);
% screen2png(['corr_tgt_r_' animal]);close;

%% explained variance on avg response between kernels
% fig = showScatterTriplets(expval_avgtgt_pop(2:4,:), ...
%     param.predictorNames, [], selectedIDs);
% screen2png(fullfile(saveServer,saveSuffix_p,['expval_avgtgt_a_' animal]));close;
% 
% fig = showScatterTriplets(100*expval_avgtgt_pop(2:4,:)./expval_avgtgt_pop(1,:), ...
%     param.predictorNames, [-200 100], selectedIDs);
% screen2png(fullfile(saveServer,saveSuffix_p,['expval_avgtgt_r_' animal]));close;

%% correlation on avg response between kernels
fig = showScatterTriplets(corr_avgtgt_pop(2:4,:), ...
    param.predictorNames, [], selectedIDs);
screen2png(fullfile(saveServer,saveSuffix_p,['corr_avgtgt_a_' animal]));close;

fig = showScatterTriplets(100*corr_avgtgt_pop(2:4,:)./corr_avgtgt_pop(1,:), ...
    param.predictorNames, [-20 100], selectedIDs);
screen2png(fullfile(saveServer,saveSuffix_p,['corr_avgtgt_r_' animal]));close;


%% show average kernel before centering
[f, kernel_avg] = showKernel3(kernel_pop, tlags, param.cardinalDir, 0);
savePaperFigure(f,fullfile(saveServer,saveSuffix_p,['avgKernel_' animal]));


%% centerring by preferred direction
tgtRange = [0.05 0.15; 0.03 0.25; -0.1 0.1];
[f, kernel_centered_avg] = showKernel3(kernel_pop, tlags, param.cardinalDir, 1, tgtRange);
savePaperFigure(f,fullfile(saveServer,saveSuffix_p,['avgKernel_centered_' animal]));


%% preferred direction across 3 kernels
 [prefDir, amp] = getKernelPrefDirAmp(kernel_pop, tlags, tgtRange, param.cardinalDir);
tuned = amp>param.ampTh;

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

% %plot3(prefdir{1},prefdir{2},prefdir{3},'.');
% subplot(131);
% %sigUnits = (prefdirPval{1}<0.05) & (prefdirPval{2}<0.05); %EMPTY
% plot(prefdir{1}, prefdir{2},'.');
% %plot(prefdir{1}(sigUnits), prefdir{2}(sigUnits),'b.');
% [rho, pval] = circ_corrcc(prefdir{1}*pi/180,prefdir{2}*pi/180);
% title(['rho:' num2str(rho) ', pval:' num2str(pval)])
% xlabel('vision'); ylabel('eye speed');
% axis equal square;
% subplot(132);
% plot(prefdir{2},prefdir{3},'.');
% [rho, pval] = circ_corrcc(prefdir{2}*pi/180,prefdir{3}*pi/180);
% title(['rho:' num2str(rho) ', pval:' num2str(pval)])
% xlabel('eye speed'); ylabel('eye position');
% axis equal square;
% subplot(133);
% plot(prefdir{3},prefdir{1},'.');
% [rho, pval] = circ_corrcc(prefdir{3}*pi/180,prefdir{1}*pi/180);
% title(['rho:' num2str(rho) ', pval:' num2str(pval)])
% xlabel('eye position'); ylabel('vision');
% axis equal square;
%screen2png(['prefDirCorr_' animal]);
% circular correlation


%% pupil dilation/constriction
% pupilLabels = {'dilation st','constriction st'};
% psthNames_pupil = cat(2,{'psth','alpha','alphaAmp','predicted_all'},param.predictorNames ,{'pdiam'});
% nChannels = size(avgPupilResp_pop,4);
% mPupilResp = squeeze(mean(avgPupilResp_pop,4));
% sePupilResp = 1/sqrt(nChannels)*squeeze(std(avgPupilResp_pop,[],4));
% nvars = size(avgPupilResp_pop,2);
% figure('position',[0 0 1000 1000]);
% ax=[];
% for ivar = 1:nvars
%     ax(ivar)=subplot(5, 2, ivar);
%     %imagesc(winSamps, param.cardinalDir, squeeze(avgDlResp(:,ivar,:)));
%     errorbar(repmat(winSamps_sacc',[1 length(pupilLabels)]), ...
%         squeeze(mPupilResp(:,ivar,:))',squeeze(sePupilResp(:,ivar,:))');
%
%     %set(gca, 'ytick',1:4,'yticklabel',pupilLabels);
%     %xlabel('time from pupil onset [s]');
%     title(psthNames_pupil{ivar});
%     xlim([-0.3 0.3]);
%     vline(0);
% end
% %linksubaxes('y',ax(1:nvars-1));
% legend(pupilLabels);
% marginplots;
% screen2png('pupilOn_pop');

%% saccade resp. recorded vs predicted
% psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);
% avgSacc = permute(squeeze(nanmean(avgSaccResp_pop,4)),[3 1 2]);
% crange = prctile(avgSacc(:),[1 99]);
% for ipred = 1:size(avgSacc,3)
%     subplot(7,1,ipred);
%     imagesc(winSamps_sacc, param.cardinalDir, squeeze(avgSacc(:,:,ipred))');
%     ylabel(psthNames{ipred})
%     %caxis(crange);
%     vline(0);
%     mcolorbar(gca,.5);
% end
% xlabel('time from saccade onset [s]');
% screen2png('Sacc_pop');
% savePaperFigure(gcf,'Sacc_pop');
%
%
% %% saccade resp aligned to preferred saccde direction. recorded vs predicted
% tgtTimes = intersect(find(winSamps>0.03), find(winSamps<0.25));
% psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);
% pavgSaccResp = permute(avgSaccResp_pop, [3 1 4 2]);
% [centeredDir, centeredData_sacc]  = alignMtxDir(pavgSaccResp, tgtTimes, param.cardinalDir);
% %[time x direction x channels x predictors]
% avgCentSacc = squeeze(nanmean(centeredData_sacc,3));
% crange = prctile(avgCentSacc(:),[1 99]);
% for ipred = 1:size(avgCentSacc,3)
%     subplot(6,1,ipred);
%     imagesc(winSamps, centeredDir, squeeze(avgCentSacc(:,:,ipred))');
%     ylabel(psthNames{ipred})
%     caxis(crange);
%     vline(0);
%     mcolorbar(gca,.5);
% end
% xlabel('time from saccade onset [s]');
% screen2png('centeredSacc_pop');
% savePaperFigure(gcf,'centeredSacc_pop');

%
% %% target resp aligned to preferred target direction
% pavgTgtResp = permute(avgTgtResp_pop, [3 1 4 2]);
% [centeredDir, centeredData_tgt]  = alignMtxDir(pavgTgtResp, tgtTimes, param.cardinalDir);
% %[time x direction x channels x predictors]
% avgCentTgt = squeeze(nanmean(centeredData_tgt,3));
% crange = prctile(avgCentTgt(:),[1 99]);
% for ipred = 1:size(avgCentTgt,3)
%     subplot(6,1,ipred);
%     imagesc(winSamps, centeredDir, squeeze(avgCentTgt(:,:,ipred))');
%     ylabel(psthNames{ipred})
%     caxis(crange);
%     vline(0);
%     mcolorbar(gca,.5);
% end
% xlabel('time from target onset [s]');
% screen2png('centeredTgt_pop');
% savePaperFigure(gcf,'centeredTgt_pop');
%
%
% %% kernel aligned to preferred direction
% %centerBin = 4;
% %centeredDir = 180/pi*circ_dist(pi/180*param.cardinalDir, pi/180*param.cardinalDir(centerBin));
% [centeredDir, centeredData_vis]  = alignMtxDir(kernel_pop(:,1:8,:), tgtTimes, param.cardinalDir);
% [~, centeredData_eye]  = alignMtxDir(kernel_pop(:,9:16,:), tgtTimes, param.cardinalDir);
%
%
% figure('position',[680   276   846   702]);
% subplot(221);
% thisImage = mean(centeredData_vis,3)';
% crange = max(abs(thisImage(:)));
% imagesc(kerneltlags, centeredDir, thisImage);
% vline(0);
% caxis([-crange crange]);
% mcolorbar(gca,.5);
% set(gca,'ytick', centeredDir);
% xlabel('time from target onset [s]');
% ylabel('relative target direction [deg]');
%
% subplot(223);
% thisImage = mean(centeredData_eye,3)';
% imagesc(kerneltlags, centeredDir, thisImage);
% crange = max(abs(thisImage(:)));
% caxis([-crange crange]);
% vline(0);
% mcolorbar(gca,.5);
% set(gca,'ytick',centeredDir);
% xlabel('time from eye movement [s]');
% ylabel('relative eye direction [deg]');
%
% subplot(424);
% plot(kerneltlags, squeeze(kernel_pop(:,17,:)), 'color',[.5 .5 .5]);
% hold on;
% plot(kerneltlags, mean(kernel_pop(:,17,:),3), 'linewidth',2);
% xlabel('time from pupil dilation [s]');
% axis tight
% vline(0);
%
% subplot(422);
% plot(kerneltlags, squeeze(kernel_pop(:,18,:)), 'color',[.5 .5 .5]);
% hold on;
% plot(kerneltlags, mean(kernel_pop(:,18,:),3), 'linewidth',2);
% xlabel('time from blink onset [s]');
% axis tight
% vline(0);
%
% subplot(426);
% plot(kerneltlags, squeeze(kernel_pop(:,19,:)), 'color',[.5 .5 .5]);
% hold on;
% plot(kerneltlags, mean(kernel_pop(:,19,:),3), 'linewidth',2);
% xlabel('time from reward on [s]');
% axis tight
% vline(0);
%
% subplot(428);
% plot(kerneltlags, squeeze(kernel_pop(:,20,:)), 'color',[.5 .5 .5]);
% hold on;
% plot(kerneltlags, mean(kernel_pop(:,20,:),3), 'linewidth',2);
% xlabel('time from punish onset [s]');
% axis tight
% vline(0);
%
% screen2png('centeredkernel_pop');
% savePaperFigure(gcf,'centeredkernel_pop');



% %% explained variance, correlation
% subplot(211);
% plot(mFiringRate_pop, expval_pop, '.');
% xlabel('mean firing rate [Hz]');
% ylabel('explained variance [%]');
% axis square;
% marginplot;
%
% subplot(212);
% plot(mFiringRate_pop, corrcoef_pop, '.');
% xlabel('mean firing rate [Hz]');
% ylabel('correlation coef');
% axis square
% marginplot;
% screen2png('mFiringRate_corrcoef_expVar_pop');
% savePaperFigure(gcf,'mFiringRate_corrcoef_expVar_pop');


% %% powerspectrum
% subplot(121);
% semilogy(faxis_common, pspec_parea_pop, 'color',[.5 .5 .5]);
% hold on
% semilogy(faxis_common, mean(pspec_parea_pop,2), 'linewidth',2);
% xlim([0 25]);
% xlabel('Frequency [Hz]'); ylabel('parea PSD');
%
% subplot(122);
% semilogy(faxis_common, pspec_psth_pop, 'color',[.5 .5 .5]);
% hold on
% semilogy(faxis_common, mean(pspec_psth_pop,2), 'linewidth',2);
% xlabel('Frequency [Hz]'); ylabel('psth PSD');
% xlim([0 25]);
%
% screen2png('powerspectrum_pop');


%% direction tuning DOITAGAIN
% subplot(121);
% plot(cmean_pop(1,), cmean_pred_pop, '.');
% squareplot;
% marginplot;
% title('Preferred direction (cirular mean)');
% xlabel('measured [rad]');
% ylabel('fitted [rad]');
%
% subplot(122);
% plot(cvar_pop, cvar_pred_pop, '.');
% squareplot;
% marginplot;
% title('Tuning width (cirular variance)');
% xlabel('measured [rad]');
% ylabel('fitted [rad]');
%
% screen2png('cmean_var_pop');

% %% alpha power
% corrbins = -1:.04:1;
% for ifreq = 1:4
%     subplot(4,2,2*ifreq-1);
%     histogram(corrcoef_parea_spk_pop(ifreq,:), corrbins);
%     [~,Pspk(ifreq)]=ttest(corrcoef_parea_spk_pop(ifreq,:));
%     vline(0);
%     title(['ttest p:' num2str(Pspk(ifreq))]);
%
%     subplot(4,2,2*ifreq);
%     histogram(corrcoef_parea_alpha_pop(ifreq,:), corrbins);
%     [~,Palpha(ifreq)]=ttest(corrcoef_parea_alpha_pop(ifreq,:));
%     title(['ttest p:' num2str(Palpha(ifreq))]);
%     vline(0);
% end
% marginplots;
% savePaperFigure(gcf,'corrcoef_parea-spk-alpha');

