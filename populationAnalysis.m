%16/12/24 created from fitPSTH_pop.m
[saveServer, rootFolder] = getReady();

saveSuffix_p = ['20260514'];%['20250207'];
animals = {'hugo', 'ollie'};

n=load(fullfile(saveServer, ['param' saveSuffix_p '.mat']),'param');
param =n.param;

%% load result of mainscript_assembly.m
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
latency_p_pop = [];
stats_stratified_pop = [];
p_hm_pop = [];
spkNGRate_pop = [];
CueTrRate_pop = [];
nLatencyTrials_pref_success_pop = [];
animalid_pop = [];
prefDir_pop = [];
auc_hm_pop = [];
difflatency_pop = [];
corr_tgt_avg_pop = [];
corr_tgt_avg_rel_pop = [];
saccDirNoTask_spkOkUCue_hist_pop = [];

for aa = 1:numel(animals)
    switch aa
        case 1
            animal = animals{aa};
        case 2
            animal = animals{aa};
    end

    for yy = 1:3
        switch yy
            case 1
                thisyear = '2021';
            case 2
                thisyear = '2022';
            case 3
                thisyear = '2023';
        end
        saveFolder = fullfile(saveServer, thisyear,animal);%17/6/23

        assemblyData = fullfile(saveFolder, ['assembly' saveSuffix_p '.mat']);

        if exist(assemblyData, 'file')
            assembly = load(assemblyData);
        else
            continue;
        end

        entries =  find(1-cellfun(@isempty, assembly.id_pop));

        %% inclusion criteria
        mFiringRate_pop = cat(2,mFiringRate_pop, [assembly.mFiringRate_pop{entries}]);
        ntargetTrials_pop = cat(2,ntargetTrials_pop, [assembly.ntargetTrials_pop{entries}]);
        PtonsetResp_pop = cat(2,PtonsetResp_pop, [assembly.PtonsetResp_pop{entries}]);

        % unit stats
        id_pop = cat(2,id_pop, [assembly.id_pop(entries)]);
        animalid_pop = cat(2, animalid_pop, aa*ones(1,numel(entries)));

        kernel_pop = cat(2,kernel_pop, assembly.kernel_pop(:,entries));
        expval_tgt_pop = cat(2,expval_tgt_pop, assembly.expval_tgt_pop{entries});
        corr_tgt_tmp = [assembly.corr_tgt_pop{entries}];
        corr_tgt_pop = cat(2,corr_tgt_pop, corr_tgt_tmp); %[corr_tgt_pop; corr_tgt_tmp'];
        corr_tgt_rel_tmp = [assembly.corr_tgt_rel_pop{entries}];
        corr_tgt_rel_pop = cat(2, corr_tgt_rel_pop, corr_tgt_rel_tmp);%[corr_tgt_rel_pop; corr_tgt_rel_tmp];
        p_hm_pop = cat(2,p_hm_pop, [assembly.p_hm_pop(entries)]);
        auc_hm_pop = cat(2, auc_hm_pop, [assembly.auc_hm_pop(entries)]);
        corr_tgt_avg_tmp = [assembly.corr_tgt_avg_pop{entries}];
        corr_tgt_avg_pop = cat(2, corr_tgt_avg_pop, corr_tgt_avg_tmp);%[corr_tgt_rel_pop; corr_tgt_rel_tmp];
        corr_tgt_avg_rel_tmp = [assembly.corr_tgt_avg_rel_pop{entries}];
        corr_tgt_avg_rel_pop = cat(2, corr_tgt_avg_rel_pop, corr_tgt_avg_rel_tmp(2:4,:));%[corr_tgt_rel_pop; corr_tgt_rel_tmp];
        
        latency_r_pop = cat(2,latency_r_pop, [assembly.latency_r_pop(entries)]);
        latency_p_pop = cat(2,latency_p_pop, [assembly.latency_p_pop(entries)]);
        difflatency_pop = cat(2,difflatency_pop, [assembly.difflatency_pop(entries)]);
        spkNGRate_pop = cat(2,spkNGRate_pop, [assembly.spkNGRate_pop{entries}]);
        CueTrRate_pop = cat(2, CueTrRate_pop, [assembly.CueTrRate_pop{entries}]);
        nLatencyTrials_pref_success_pop = cat(2, nLatencyTrials_pref_success_pop, assembly.nLatencyTrials_pref_success_pop{entries});
        prefDir_pop = cat(2,prefDir_pop, [assembly.prefDir_pop{entries}]);

        binEdges = [param.cardinalDir 360] - 0.5*mean(diff(param.cardinalDir)); %cf.getSaccDir
        hist_tmp = nan(numel(binEdges)-1, numel(entries));
        for ee = 1:numel(entries)
            data = assembly.saccDirNoTask_spkOkUCue_pop(entries(ee));
            hist_tmp(:,ee) = histcounts(data{1}, 'BinEdges', binEdges);
        end
        saccDirNoTask_spkOkUCue_hist_pop = cat(2, saccDirNoTask_spkOkUCue_hist_pop, hist_tmp);
        
        %load just once
        tlags = assembly.tlags_pop{entries(1)};
    end
end

display(['fraction of trials whose spike rates are too high: ' num2str(mean(spkNGRate_pop))]);
display(['fraction of trials where cue was presented and excluded from analysis: ' num2str(mean(CueTrRate_pop))]);

%% apply inclusion critetia
param.corr_tgtTh = 0.5; %0.1
[~, corr_tgtOK, ntargetTrOK, ptonsetRespOK] ...
    = inclusionCriteria(corr_tgt_avg_pop(1,:), ntargetTrials_pop, PtonsetResp_pop, param);
okunits = find(ntargetTrOK .* ptonsetRespOK); % for Revision 1

%% unit rejection in both animals
nTotalUnits = numel(ntargetTrOK);
disp(['criteria 1: ' num2str(nTotalUnits - sum(ntargetTrOK))]);
disp(['criteria 2:' num2str(nTotalUnits - sum(ptonsetRespOK))]);

%% unit rejection in each animal
for aa = 1:2
    disp(animals{aa})
    nTotalUnits = numel(intersect(1:numel(ntargetTrOK), find(animalid_pop==aa)));
    disp(['amongst ' num2str(nTotalUnits) ' that passed criteria 2, ']);
    disp(['rejected units by criteria 1: ' num2str(nTotalUnits - sum(ntargetTrOK.*(animalid_pop==aa)))]);
    disp(['rejected units by criteria 2:' num2str(nTotalUnits - sum(ptonsetRespOK.*(animalid_pop==aa)))]);
    rejectRate = 100*(1 - numel(intersect(okunits, find(animalid_pop==aa)))/nTotalUnits);
    disp([num2str(rejectRate) '% units rejected in total']);
end

disp(['mean correlation that passed all the criteria:' num2str(mean(corr_tgt_avg_pop(1,okunits)))]);
disp(['std correlation that passed all the criteria:' num2str(std(corr_tgt_avg_pop(1,okunits)))]);



% units_ng4 = find(ntargetTrOK .* ptonsetRespOK .* ~corr_tgtOK);
% id_pop(units_ng4)'

% omit data with the following criteria??
% no saccade response
% low spontaneous firing
% low number of successful trials


%% retain included units
% units x 1
kernel_pop = kernel_pop(:,okunits);
corr_tgt_pop = corr_tgt_pop(:,okunits);
corr_tgt_rel_pop = corr_tgt_rel_pop(:,okunits);
id_pop = id_pop(okunits);
mFiringRate_pop = mFiringRate_pop(okunits);
latency_r_pop = latency_r_pop(okunits);
latency_p_pop = latency_p_pop(okunits);
p_hm_pop = p_hm_pop(okunits);
auc_hm_pop = auc_hm_pop(okunits);
nLatencyTrials_pref_success_pop = nLatencyTrials_pref_success_pop(okunits);
animalid_pop = animalid_pop(okunits);
prefDir_pop = prefDir_pop(okunits);
difflatency_pop = cellfun(@(a)a(1), difflatency_pop(okunits));
corr_tgt_avg_pop = corr_tgt_avg_pop(:,okunits);
corr_tgt_avg_rel_pop = corr_tgt_avg_rel_pop(:,okunits);
saccDirNoTask_spkOkUCue_hist_pop = saccDirNoTask_spkOkUCue_hist_pop(:,okunits);

disp(['total ' num2str(numel(animalid_pop)) ' units, M1: ' num2str(sum(animalid_pop==1)) ', M2:' num2str(sum(animalid_pop==2))]);
save('ppcPaper_id','id_pop'); % To Jo 5/11/2025; 17/5/2026

%% convert from cell to matrix
latency_r_nb_pop = cellfun(@(a)a(1), latency_r_pop);
latency_p_nb_pop = cellfun(@(a)a(1), latency_p_pop);
p_hm_pop = [cellfun(@(a)a(1), p_hm_pop); cellfun(@(a)a(2), p_hm_pop)];
auc_hm_pop = [cellfun(@(a)a(1), auc_hm_pop); cellfun(@(a)a(2), auc_hm_pop)];

%% selected units
theseIDs = {'hugo/2021/08August/25/27',... %vision
      'hugo/2022/07July/26/19',... %eye speed
    'hugo/2022/08August/15/4'}; %integrator new 2025 
    [~, selectedIDs] = intersect(id_pop, theseIDs);

% high latency correlation, also eye-position driven
theseIDs_lat = {'hugo/2022/07July/08/1', ... % visual, latency independent
    'hugo/2021/03March/19/8'}; %dependent
[~, selectedIDs_lat] = intersect(id_pop, theseIDs_lat);


theseIDs_hm = {  'hugo/2021/03March/23/29' ...
    'hugo/2021/09September/07/26'}; %no reduction 
[~, selectedIDs_hm] = intersect(id_pop, theseIDs_hm);

selectedIDs_lc = find(corr_tgt_avg_pop(1,:)<0.5);

%% preferred direction (FIG1)
disp('fig 1');
scriptFileName = 'script_fig1.m';
lines = readlines(scriptFileName);
for i = 1:numel(lines)
    eval(lines(i));
end
savePaperFigure(f, fullfile(saveServer, saveSuffix_p, ['hist_prefDir_' animals{:}]));
save('data_fig1','prefDir_pop','animalid_pop','lines');
close(f);

%preferred direction across 3 kernels
% tgtRange = [0.05 0.15; 0.03 0.25; -0.1 0.1];
% [f, pval] = showKernelPrefDirScatter(kernel_pop, tlags, tgtRange, param, animalid_pop);
% savePaperFigure(gcf,fullfile(saveServer,saveSuffix_p,['tuning_pop_' animals{:}]));
% close(f);

%% spontaneous saccade direction (Fig1 SUPPLEMENT)
disp('fig supplement 1');
scriptFileName = 'script_fig_supplement1.m';
lines = readlines(scriptFileName);
for i = 1:numel(lines)
    eval(lines(i));
end
savePaperFigure(f, fullfile(saveServer, saveSuffix_p, ['hist_saccDir_' animals{:}]));
save('data_fig_supplement1','saccDirNoTask_spkOkUCue_hist_pop','param','binEdges');
close(f);

%% cell type distibution (FIG2)
disp('fig 2JKL');
scriptFileName = 'script_fig2JKL.m';
lines = readlines(scriptFileName);
eval(lines); %fig_avg, fig3D_avg
savePaperFigure(fig_avg, fullfile(saveServer,saveSuffix_p,['corr_tgt_avg_' animals{:}]), 'dum');close(fig_avg);
save('data_fig2JKL','corr_tgt_avg_pop','animalid_pop','selectedIDs','param','lines');


%% histogram of correlation of the full model (FIG SUPPLEMENT 2-1)
disp('fig supplement 2-1');
figs2_1 = function_fig_supplement2_1(corr_tgt_avg_pop, animalid_pop, param);
savePaperFigure(figs2_1,fullfile(saveServer,saveSuffix_p,['hist_corr_tgt_avg_pop_' animals{:}]));
save('data_fig_supplement2_1','corr_tgt_avg_pop','animalid_pop','param');
close(figs2_1);


%% recording depth (FIG SUPPLEMENT 2-3)
disp('fig supplement 2-3');
load('ppcDepths.mat','ppcdays','ppcdepths');
[f, depth_pop] = function_fig_supplement2_3(corr_tgt_avg_pop, animalid_pop, ...
    animals, ppcdays, ppcdepths, id_pop, param);
savePaperFigure(f, fullfile(saveServer,saveSuffix_p,['p_depth_' animals{:} ]), 'dum');close;
save('data_fig_supplement2_3','corr_tgt_avg_pop','animalid_pop','animals',...
    'param','ppcdays','ppcdepths','id_pop');


%% hit v miss (FIG 3)
disp('fig 3');
scriptFileName = 'script_fig3.m';
lines = readlines(scriptFileName);
for i = 1:numel(lines)
    eval(lines(i));
end
savePaperFigure(f, fullfile(saveServer,saveSuffix_p,['p_hm_' param.predictorNames{tgtModalities(1)} ...
    '_'  param.predictorNames{tgtModalities(2)} '_' animals{:} ]), 'dum');close;
save('data_fig3','corr_tgt_avg_pop', 'auc_hm_pop', 'selectedIDs_hm', 'animalid_pop', 'tgtModalities', 'param','lines');
% disp(['The number of significant units before regression: '  num2str(sum(p_hm_pop_before<0.05))]);
% disp(['The number of significant units after regression: '  num2str(sum(p_hm_pop_after<0.05))]);
disp(['AUC of ' num2str(selectedIDs_hm(1)) ': ' num2str(auc_hm_pop(1, selectedIDs_hm(1)))]);
disp(['AUC of ' num2str(selectedIDs_hm(2)) ': ' num2str(auc_hm_pop(1, selectedIDs_hm(2)))]);


%% hit v miss (FIG SUPPLEMENT 3)
disp('fig supplement 3_2');
scriptFileName = 'script_fig_supplement3_2.m';
lines = readlines(scriptFileName);
for i = 1:numel(lines)
    eval(lines(i));
end
savePaperFigure(f, fullfile(saveServer,saveSuffix_p,['p_hm_' param.predictorNames{tgtModalities(1)} ...
    '_'  param.predictorNames{tgtModalities(2)} '_' animals{:} ]), 'dum');close;
save('data_fig_supplement3_2','corr_tgt_avg_pop', 'auc_hm_pop', 'selectedIDs_hm', 'animalid_pop', ...
    'tgtModalities', 'param','lines');
% disp(['The number of significant units before regression: '  num2str(sum(p_hm_pop_before<0.05))]);
% disp(['The number of significant units after regression: '  num2str(sum(p_hm_pop_after<0.05))]);
disp(['AUC of ' num2str(selectedIDs_hm(1)) ' afer regression: ' num2str(auc_hm_pop(2, selectedIDs_hm(1)))]);
disp(['AUC of ' num2str(selectedIDs_hm(2)) ' afer regression: ' num2str(auc_hm_pop(2, selectedIDs_hm(2)))])



%% latency stats on correlation to tgt (FIG4)
disp('fig 4');
scriptFileName = 'script_fig4CD.m';
lines = readlines(scriptFileName);
for i = 1:numel(lines)
    eval(lines(i));
end
savePaperFigure(f, fullfile(saveServer,saveSuffix_p,['p_latency_r_pref_success_' param.predictorNames{tgtModalities(1)} ...
    '_'  param.predictorNames{tgtModalities(2)} '_' animals{:}]),'dum'); close(f);
save('data_fig4CD','corr_tgt_avg_pop', 'latency_r_nb_pop', 'latency_p_nb_pop',  ...
    'nLatencyTrials_pref_success_pop', 'param', 'selectedIDs_lat', 'animalid_pop','tgtModalities','lines');
disp(['unit ID: ' num2str(selectedIDs_lat(1)) ', r: ' num2str(latency_r_nb_pop(selectedIDs_lat(1))) ', p: ' num2str(latency_p_nb_pop(selectedIDs_lat(1)))])
disp(['unit ID: ' num2str(selectedIDs_lat(2)) ', r: ' num2str(latency_r_nb_pop(selectedIDs_lat(2))) ', p: ' num2str(latency_p_nb_pop(selectedIDs_lat(2)))])

% histogram of neural - behavioural latencies
scriptFileName = 'script_fig4E.m';
lines = readlines(scriptFileName);
for i = 1:numel(lines)
    eval(lines(i));
end
savePaperFigure(f, fullfile(saveServer,saveSuffix_p,['hist_difflatency_' animals{:}]),'dum'); close(f);
save('data_fig4E','difflatency_pop','nLatencyTrials_pref_success_pop', 'latency_r_nb_pop','latency_p_nb_pop', 'param','tgtUnits','lines');


%% latency stats on correlation to tgt (FIG SUPPLEMENT 4)
disp('fig supplement 4');
scriptFileName = 'script_fig_supplement4.m';
lines = readlines(scriptFileName);
for i = 1:numel(lines)
    eval(lines(i));
end
savePaperFigure(f, fullfile(saveServer,saveSuffix_p,['p_latency_r_pref_success_' param.predictorNames{tgtModalities(1)} ...
    '_'  param.predictorNames{tgtModalities(2)} '_' animals{:}]),'dum'); close(f);
save('data_fig_supplement4','corr_tgt_avg_pop', 'latency_r_nb_pop', 'latency_p_nb_pop',  'nLatencyTrials_pref_success_pop', 'param', 'selectedIDs_lat', 'animalid_pop', 'tgtModalities','lines');

