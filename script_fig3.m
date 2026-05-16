%load('data_fig3.mat');

param.auc_th = .8;
tgtModalities = [1 2];
f = showHMScatter(corr_tgt_avg_pop(2:4,:), auc_hm_pop, selectedIDs_hm, animalid_pop, tgtModalities, param);
