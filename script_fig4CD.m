%load('data_fig4CD.mat')
param.nLatencyTr_pref_th = 3;
param.r_latency_th = 0.5; 
param.p_latency_th = 0.05; 
tgtModalities = [1 2];
f = showLatencyScatter(corr_tgt_avg_pop(2:4,:), latency_r_nb_pop,   latency_p_nb_pop,  nLatencyTrials_pref_success_pop, param, selectedIDs_lat, animalid_pop, tgtModalities);
