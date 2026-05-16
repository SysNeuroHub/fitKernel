%equivalent of figure 2
f2 = showScatterTriplets(corr_tgt_avg_pop(2:4,selectedIDs_lc), param.predictorNames, [-0.4 1], [],'linear',animalid_pop(selectedIDs_lc));
savePaperFigure(f2, fullfile(saveServer,saveSuffix_p,['corr_tgt_avg_' animals{:} '_lc']), 'dum');close(fig_avg);

%equivalent of figure 3
param.auc_th = .8;
tgtModalities = [1 2];
f3 = showHMScatter(corr_tgt_avg_pop(2:4,selectedIDs_lc), auc_hm_pop(:,selectedIDs_lc), [], animalid_pop(selectedIDs_lc), tgtModalities, param);
ylim(gca,[0 12]);
savePaperFigure(f3, fullfile(saveServer,saveSuffix_p,['p_hm_' param.predictorNames{tgtModalities(1)} ...
    '_'  param.predictorNames{tgtModalities(2)} '_' animals{:} '_lc']), 'dum');close;

%equivalent of figure4
param.nLatencyTr_pref_th = 3;
param.r_latency_th = 0.5; 
param.p_latency_th = 0.05; 
tgtModalities = [1 2];
f4 = showLatencyScatter(corr_tgt_avg_pop(2:4,selectedIDs_lc), latency_r_nb_pop(selectedIDs_lc),   latency_p_nb_pop(selectedIDs_lc),  nLatencyTrials_pref_success_pop(selectedIDs_lc), param, [], animalid_pop(selectedIDs_lc), tgtModalities);
ylim(gca,[0 12]);
savePaperFigure(f4, fullfile(saveServer,saveSuffix_p,['p_latency_r_pref_success_' param.predictorNames{tgtModalities(1)} ...
    '_'  param.predictorNames{tgtModalities(2)} '_' animals{:} '_lc']),'dum'); close(f4);
