param.nLatencyTr_pref_th = 3;
param.r_latency_th = 0.5; 
param.p_latency_th = 0.05; 
f = figure('position',[0 0 200 200]);
tgtUnits = find(nLatencyTrials_pref_success_pop >= param.nLatencyTr_pref_th   & latency_r_nb_pop > param.r_latency_th & latency_p_nb_pop < param.p_latency_th);
histogram(1e3*difflatency_pop(tgtUnits), 10, 'facecolor', "#7E2F8E");
 axis square; box off;
 set(gca,'tickdir','out');
 xlabel('Neural - behavioural latency [ms]');ylabel('#Units');title(['n=' num2str(numel(tgtUnits)) ', median=' num2str(median(1e3*difflatency_pop(tgtUnits)))]);
