%16/12/24 created from fitPSTH_pop.m
[saveServer, rootFolder] = getReady();

saveSuffix_p = ['20241212'];
animals = {'hugo', 'ollie'};

n=load(fullfile(saveServer,'param20241223.mat'),'param');
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
nLatencyTrials_pref_success_pop = [];
animalid_pop = [];

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
        latency_r_pop = cat(2,latency_r_pop, [assembly.latency_r_pop(entries)]);
        spkNGRate_pop = cat(2,spkNGRate_pop, [assembly.spkNGRate_pop{entries}]);
        CueTrRate_pop = cat(2, CueTrRate_pop, [assembly.CueTrRate_pop{entries}]);
        nLatencyTrials_pref_success_pop = cat(2, nLatencyTrials_pref_success_pop, assembly.nLatencyTrials_pref_success_pop{entries});

        %load just once
        tlags = assembly.tlags_pop{entries(1)};
    end
end


%% apply inclusion critetia
% param.corr_tgtTh = 0.05;
param.ptonsetRespTh=.2; %.1;
% param.nLatencyTr_pref_th = 5;
[okunits, corr_tgtOK, ntargetTrOK, ptonsetRespOK] ...
    = inclusionCriteria(corr_tgt_pop(1,:), ntargetTrials_pop, PtonsetResp_pop, param);

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
p_hm_pop = p_hm_pop(okunits);
nLatencyTrials_pref_success_pop = nLatencyTrials_pref_success_pop(okunits);
animalid_pop = animalid_pop(okunits);

%% convert from cell to matrix
latency_r_nb_pop = cellfun(@(a)a(1), latency_r_pop);
p_hm_pop_before = cellfun(@(a)a(1), p_hm_pop);
p_hm_pop_after = cellfun(@(a)a(2), p_hm_pop);
p_hm_pop = [p_hm_pop_before; p_hm_pop_after];


%% selected units
theseIDs = {'hugo/2021/08August/25/27',... %vision
    'hugo/2022/07July/26/19',... %eye speed
    'hugo/2022/08August/15/4'}; %integrator new 2025 
[~, selectedIDs] = intersect(id_pop, theseIDs);
%integrator candidates
% 'hugo/2022/08August/05/2' 2023 JNS
% 'hugo/2022/08August/23/7'
% 'hugo/2022/08August/9/3';
% 'hugo/2022/08August/12/1';
% 'hugo/2022/08August/15/5';
% 'hugo/2022/08August/05/2'; %integrater used in 2023 JNS
% 'hugo/2021/12December/14/13'
% 'hugo/2021/12December/07/6'

% high latency correlation, also eye-position driven
theseIDs_lat = {'hugo/2022/07July/08/1', ... % visual, latency independent
    'hugo/2022/09September/15/3',... %dependent
    };
[~, selectedIDs_lat] = intersect(id_pop, theseIDs_lat);
% latency candidates
% 'hugo/2021/08August/25/27',... %vision example unit, latency independent
% 'hugo/2022/09September/8/11' %excluded by inclusionCriteria
%    'hugo/2022/08August/05/2', ... %integrator, latency independent
% 'hugo/2022/08August/23/2',... %latency dependent
%    'hugo/2022/08August/15/3',... %eye position, mildly latency dependent
%    'hugo/2021/09September/07/23'}; %dependent
% 'hugo/2023/02February/02/11',... %independent ... excluded by inclusion criteria
%    'hugo/2021/12December/10/26' %independent, eye driven


theseIDs_hm = {    'hugo/2021/09September/07/26',... %no reduction
    'hugo/2022/03March/01/29'}; %massive reduction
[~, selectedIDs_hm] = intersect(id_pop, theseIDs_hm);
% H&M candidates
% 'hugo/2022/03March/01/29' %massive reduction in d' after regression
% 'hugo/2022/03March/21/13' % incomplete reduction in d'
% 'hugo/2022/07July/07/17' %massive reduction in d'
% 'hugo/2021/03March/25/24'; %no reduction in d'
% 'hugo/2022/03March/01/29' %massive reduction in d'
% 'hugo/2022/03March/15/29' %excluded by inclusionCriteria
%     'hugo/2022/03March/21/13',... %massive reduction

%% kernel avg across units (FIG2)
%average kernel before centering
[f, kernel_avg] = showKernel3(kernel_pop, tlags, param.cardinalDir, 0);
savePaperFigure(f,fullfile(saveServer,saveSuffix_p,['avgKernel_' animals{:}]));
close(f);

%centerring by preferred direction
tgtRange = [0.05 0.15; 0.03 0.25; -0.1 0.1];
[f, kernel_centered_avg] = showKernel3(kernel_pop, tlags, param.cardinalDir, 1, tgtRange);
savePaperFigure(f, fullfile(saveServer,saveSuffix_p,['avgKernel_centered_' animals{:}]));
close(f);

%preferred direction across 3 kernels
[f, pval] = showKernelPrefDirScatter(kernel_pop, tlags, tgtRange, param, animalid_pop);
savePaperFigure(gcf,fullfile(saveServer,saveSuffix_p,['tuning_pop_' animals{:}]));
close(f);


%% cell type distibution (FIG3)

% correlation during target presentation
fig_abs = showScatterTriplets(corr_tgt_pop, ...
    param.predictorNames, [-.2 1], selectedIDs,'linear',animalid_pop);
squareplots(fig_abs, [-.2 1]);
savePaperFigure(fig_abs, fullfile(saveServer,saveSuffix_p,['corr_tgt_' animals{:}]));close(fig_abs);

fig_rel = showScatterTriplets(corr_tgt_rel_pop, ...
    param.predictorNames, [-50 150], selectedIDs,'linear',animalid_pop);
squareplots(fig_rel, [-50 150]);
savePaperFigure(fig_rel, fullfile(saveServer,saveSuffix_p,['corr_tgt_r_' animals{:}]));close(fig_rel);



%% latency stats on correlation to tgt (FIG4)
param.nLatencyTr_pref_th = 7;
f = showLatencyScatter(corr_tgt_rel_pop, latency_r_nb_pop, nLatencyTrials_pref_success_pop, ...
    param, selectedIDs_lat, animalid_pop);
screen2png(fullfile(saveServer,saveSuffix_p,['corr_tgt_rel_pop_latency_r_pref_success_' animals{:}]));close(f);



%% hit v miss (FIG 5)
f = showHMScatter(corr_tgt_rel_pop, p_hm_pop, selectedIDs_hm, animalid_pop);
screen2png(fullfile(saveServer,saveSuffix_p,['p_hm' animals{:}]));close;
