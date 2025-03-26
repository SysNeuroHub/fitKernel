%16/12/24 created from fitPSTH_pop.m
[saveServer, rootFolder] = getReady();

saveSuffix_p = ['20250207'];
animals = {'hugo', 'ollie'};

n=load(fullfile(saveServer, ['param' saveSuffix_p '.mat']),'param');
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
latency_p_pop = [];
stats_stratified_pop = [];
p_hm_pop = [];
spkNGRate_pop = [];
CueTrRate_pop = [];
nLatencyTrials_pref_success_pop = [];
animalid_pop = [];
prefDir_pop = [];
ranksumval_hm_pop = [];
ranksumz_hm_pop = [];
difflatency_pop = [];

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
        ranksumval_hm_pop = cat(2, ranksumval_hm_pop, [assembly.ranksumval_hm_pop(entries)]);
        ranksumz_hm_pop = cat(2, ranksumz_hm_pop, [assembly.ranksumz_hm_pop(entries)]);
        
        latency_r_pop = cat(2,latency_r_pop, [assembly.latency_r_pop(entries)]);
        latency_p_pop = cat(2,latency_p_pop, [assembly.latency_r_pop(entries)]);
        difflatency_pop = cat(2,difflatency_pop, [assembly.difflatency_pop(entries)]);
        spkNGRate_pop = cat(2,spkNGRate_pop, [assembly.spkNGRate_pop{entries}]);
        CueTrRate_pop = cat(2, CueTrRate_pop, [assembly.CueTrRate_pop{entries}]);
        nLatencyTrials_pref_success_pop = cat(2, nLatencyTrials_pref_success_pop, assembly.nLatencyTrials_pref_success_pop{entries});
        prefDir_pop = cat(2,prefDir_pop, [assembly.prefDir_pop{entries}]);

        %load just once
        tlags = assembly.tlags_pop{entries(1)};
    end
end

display(['fraction of trials whose spike rates are too high: ' num2str(mean(spkNGRate_pop))]);
display(['fraction of trials where cue was presented and excluded from analysis: ' num2str(mean(CueTrRate_pop))]);

%% apply inclusion critetia
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
latency_p_pop = latency_p_pop(okunits);
p_hm_pop = p_hm_pop(okunits);
ranksumval_hm_pop = ranksumval_hm_pop(okunits);
ranksumz_hm_pop = ranksumz_hm_pop(okunits);
nLatencyTrials_pref_success_pop = nLatencyTrials_pref_success_pop(okunits);
animalid_pop = animalid_pop(okunits);
prefDir_pop = prefDir_pop(okunits);
difflatency_pop = cellfun(@(a)a(1), difflatency_pop(okunits));

%% convert from cell to matrix
latency_r_nb_pop = cellfun(@(a)a(1), latency_r_pop);
latency_p_nb_pop = cellfun(@(a)a(1), latency_p_pop);
p_hm_pop = [cellfun(@(a)a(1), p_hm_pop); cellfun(@(a)a(2), p_hm_pop)];
ranksumval_hm_pop = [cellfun(@(a)a(1), ranksumval_hm_pop); cellfun(@(a)a(2), ranksumval_hm_pop)];
ranksumz_hm_pop = [cellfun(@(a)a(1), ranksumz_hm_pop); cellfun(@(a)a(2), ranksumz_hm_pop)];

%% selected units
theseIDs = {'hugo/2021/08August/25/27',... %vision
      'hugo/2022/07July/26/19',... %eye speed
    'hugo/2022/08August/15/4'}; %integrator new 2025 
%'hugo/2022/08August/17/3',... %eye speed+position
    [~, selectedIDs] = intersect(id_pop, theseIDs);

%    'hugo/2021/03March/22/8',... %eye position
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
    'hugo/2021/03March/19/8'}; %dependent
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
% 'hugo/2022/09September/15/3',... %dependent

theseIDs_hm = {  'hugo/2021/03March/23/29' ...
    'hugo/2021/09September/07/26'}; %no reduction 
[~, selectedIDs_hm] = intersect(id_pop, theseIDs_hm);
% H&M candidates
%'hugo/2021/03March/25/24',... %no reduction but NG
% 'hugo/2022/07July/27/26' ... NG
% 'hugo/2021/09September/07/26' ... NG
% 'hugo/2022/03March/01/29' %massive reduction but NG
% 'hugo/2022/03March/21/13' % incomplete reduction in d'
% 'hugo/2022/07July/07/17' %massive reduction in d'
% 'hugo/2021/03March/25/24'; %no reduction in d'
% 'hugo/2022/03March/15/29' %excluded by inclusionCriteria
%     'hugo/2022/03March/21/13',... %massive reduction
% candidates for massive reduction
    % {'hugo/2021/03March/18/13'    }
    % {'hugo/2021/03March/19/4'     }
    % {'hugo/2021/03March/19/7'     } good
    % {'hugo/2021/03March/22/20'    }
    % {'hugo/2021/03March/23/29'    } good
    % {'hugo/2021/03March/29/20'    }
    % {'hugo/2021/03March/29/7'     }
    % {'hugo/2021/03March/29/8'     }
    % {'hugo/2022/07July/29/28'     } good
    % {'hugo/2022/07July/29/30'     }
    % {'hugo/2022/08August/17/1'    } good
    % {'hugo/2022/08August/18/1'    } good
    % {'hugo/2022/08August/18/10'   }
    % {'hugo/2022/08August/18/2'    }
    % {'hugo/2022/09September/08/11'}
    % {'hugo/2022/09September/15/3' }

%% preferred direction (FIG1)
f = figure('position', [0 0 168 168]);
histogram(prefDir_pop(animalid_pop == 1), 0:5:359,'FaceColor',[.5 .5 .5]); hold on;
histogram(prefDir_pop(animalid_pop == 2), 0:5:359,'FaceColor',[.1 .1 .1]);
axis tight square; box off;
ylim([0 35])
set(gca,'tickdir','out','xtick',[0 90 180 270 360]);
title(['n=' num2str(numel(animalid_pop))]);
xlabel('Target stimulus direction [deg]'); ylabel('# Units');
legend('M1','M2');
savePaperFigure(f, fullfile(saveServer, saveSuffix_p, ['hist_prefDir_' animals{:}]));
close(f);

%preferred direction across 3 kernels
tgtRange = [0.05 0.15; 0.03 0.25; -0.1 0.1];
[f, pval] = showKernelPrefDirScatter(kernel_pop, tlags, tgtRange, param, animalid_pop);
savePaperFigure(gcf,fullfile(saveServer,saveSuffix_p,['tuning_pop_' animals{:}]));
close(f);

%% kernel avg across units (FIG2)
% %average kernel before centering
% [f, kernel_avg] = showKernel_pop(kernel_pop, tlags, param.cardinalDir, 0);
% savePaperFigure(f,fullfile(saveServer,saveSuffix_p,['avgKernel_' animals{:}]));
% close(f);
% 
% %centerring by preferred direction
% tgtRange = [0.05 0.15; 0.03 0.25; -0.1 0.1];
% [f, kernel_centered_avg] = showKernel_pop(kernel_pop, tlags, param.cardinalDir, 1, tgtRange);
% savePaperFigure(f, fullfile(saveServer,saveSuffix_p,['avgKernel_centered_' animals{:}]));
% close(f);

% histogram of correlation of the full model
f = figure('position', [0 0 200 200]);
histogram(corr_tgt_pop(1,animalid_pop==1),.05:.02:.8,'FaceColor',[.5 .5 .5]); hold on;
histogram(corr_tgt_pop(1,animalid_pop==2),.05:.02:.8,'FaceColor',[.1 .1 .1])
axis tight square 
box off
set(gca,'tickdir','out')
xlabel('correlation of the full model'); ylabel('# units');
legend('M1','M2');
savePaperFigure(gcf,fullfile(saveServer,saveSuffix_p,['hist_corr_tgt_pop_' animals{:}]));
close(f);

%% cell type distibution (FIG2)
% correlation during target presentation
% fig_abs = showScatterTriplets(corr_tgt_pop, ...
%     param.predictorNames, [-.2 1], selectedIDs,'linear',animalid_pop);
% squareplots(fig_abs, [-.2 1]);
% savePaperFigure(fig_abs, fullfile(saveServer,saveSuffix_p,['corr_tgt_' animals{:}]));close(fig_abs);

fig_rel = showScatterTriplets(corr_tgt_rel_pop, ...
    param.predictorNames, [-50 120], selectedIDs,'linear',animalid_pop);
% squareplots(fig_rel, [-50 120]);
savePaperFigure(fig_rel, fullfile(saveServer,saveSuffix_p,['corr_tgt_r_' animals{:}]), 'dum');close(fig_rel);


%% hit v miss (FIG 3)
param.ranksumz_th = 2;
tgtModalities = [1 2];
f = showHMScatter(corr_tgt_rel_pop, ranksumz_hm_pop, selectedIDs_hm, animalid_pop, tgtModalities, param);
savePaperFigure(f, fullfile(saveServer,saveSuffix_p,['p_hm_' param.predictorNames{tgtModalities(1)} ...
    '_'  param.predictorNames{tgtModalities(2)} '_' animals{:} ]), 'dum');close;
% disp(['The number of significant units before regression: '  num2str(sum(p_hm_pop_before<0.05))]);
% disp(['The number of significant units after regression: '  num2str(sum(p_hm_pop_after<0.05))]);

tgtModalities = [1 3];
f = showHMScatter(corr_tgt_rel_pop, ranksumz_hm_pop, [], animalid_pop, tgtModalities, param);
savePaperFigure(f, fullfile(saveServer,saveSuffix_p,['p_hm_' param.predictorNames{tgtModalities(1)} ...
    '_'  param.predictorNames{tgtModalities(2)} '_' animals{:} ]), 'dum');close;
% disp(['The number of significant units before regression: '  num2str(sum(p_hm_pop_before<0.05))]);
% disp(['The number of significant units after regression: '  num2str(sum(p_hm_pop_after<0.05))]);



%% latency stats on correlation to tgt (FIG4)
param.nLatencyTr_pref_th = 7;
param.r_latency_th = 0.5; 
tgtModalities = [1 2];
f = showLatencyScatter(corr_tgt_rel_pop, latency_r_nb_pop,   ...
    nLatencyTrials_pref_success_pop, param, selectedIDs_lat, animalid_pop, tgtModalities);
savePaperFigure(f, fullfile(saveServer,saveSuffix_p,['p_latency_r_pref_success_' param.predictorNames{tgtModalities(1)} ...
    '_'  param.predictorNames{tgtModalities(2)} '_' animals{:}]),'dum'); close(f);

f = figure('position',[0 0 200 200]);
tgtUnits = find(latency_r_nb_pop > param.r_latency_th &...
    nLatencyTrials_pref_success_pop > param.nLatencyTr_pref_th);
histogram(1e3*difflatency_pop(tgtUnits), 10, 'facecolor', "#7E2F8E");
 axis square; box off;
 set(gca,'tickdir','out');
 xlabel('Neural - behavioural latency [ms]');ylabel('#Units');title(['n=' num2str(numel(tgtUnits))]);
 savePaperFigure(f, fullfile(saveServer,saveSuffix_p,['hist_difflatency_' animals{:}]),'dum'); close(f);


tgtModalities = [1 3];
f = showLatencyScatter(corr_tgt_rel_pop, latency_r_nb_pop,   ...
    nLatencyTrials_pref_success_pop, param, selectedIDs_lat, animalid_pop, tgtModalities);
savePaperFigure(f, fullfile(saveServer,saveSuffix_p,['p_latency_r_pref_success_' param.predictorNames{tgtModalities(1)} ...
    '_'  param.predictorNames{tgtModalities(2)} '_' animals{:}]),'dum'); close(f);


