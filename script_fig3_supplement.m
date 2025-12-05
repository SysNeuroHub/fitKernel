param.auc_th = .8;
tgtModalities = [1 3];
f = showHMScatter(corr_tgt_avg_pop(2:4,:), auc_hm_pop, [], animalid_pop, tgtModalities, param);
