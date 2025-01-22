function f = showLatencyScatter(corr_tgt_rel_pop, latency_r_nb_pop, nLatencyTrials_pref_pop, param, selectedIDs, animalid_pop)
%f = showLatencyScatter(corr_tgt_rel_pop, latency_r_nb_pop, param, selectedIDs, animalid_pop)
%single-trial latency correlation 2024June
%
%  latency_r_nb_pop: neuro-bhv latency correlation
% nLatencyTrials_pref_pop: number of trials in which latency was successfully detected

nUnits = size(corr_tgt_rel_pop,2);
msize = 30;

if nargin < 6
    animalid_pop = ones(1, nUnits);
end
if nargin < 5
    selectedIDs = [];
end

f = figure('position',[0 0 1000 1000]);


%latencyOK = (sum(isnan(stratified_latency_pop),1) < 1);
units_latency = nLatencyTrials_pref_pop >= param.nLatencyTr_pref_th;
units_selectedIDs = zeros(1, nUnits);
units_selectedIDs(selectedIDs) = 1;

for aa = 1:numel(unique(animalid_pop))
    switch aa
        case 1
            asymbol = 'o'; acolor = [.5 .5 .5];%'b';
        case 2
            asymbol = 'diamond'; acolor = [.5 .5 .5];%'r';
    end

    % all units
    s0 = scatter(corr_tgt_rel_pop(2, animalid_pop==aa), corr_tgt_rel_pop(1, animalid_pop==aa), ...
        msize, [.5 .5 .5], asymbol);
    s0.MarkerEdgeColor = acolor;
    hold on;

    % units with sufficient number of detected latency trials
    s1 = scatter(corr_tgt_rel_pop(2,units_latency.*animalid_pop==aa), corr_tgt_rel_pop(1, units_latency.*animalid_pop==aa),...
        msize, latency_r_nb_pop(units_latency.*animalid_pop==aa), 'filled', asymbol);
    s1.MarkerEdgeColor = acolor;
    s1.MarkerFaceAlpha = .5;

    % units used as examples
    s2 = scatter(corr_tgt_rel_pop(2, units_selectedIDs.*animalid_pop==aa), corr_tgt_rel_pop(1, units_selectedIDs.*animalid_pop==aa), ...
        msize, [0 0 0], asymbol, 'linewidth', 2);%stratified_avgCorr_pop(selectedIDs));
    end
clim([-0 .6]);
colormap("cool");colorbar;
% colormap("winter");colorbar;
%colormap(flipud(summer));colorbar;
xlabel('eye speed'); ylabel('stimulus');
title(sprintf('corr to tgt, relative to full mdl\nlatency correlation (success pref)'));
squareplot(gca, [-50 120]);
set(gca,'tickdir','out');
