function f = showLatencyScatter(corr_tgt_rel_pop, latency_r_nb_pop,  ...
    nLatencyTrials_pref_pop, param, selectedIDs, animalid_pop, tgtModalities)
%f = showLatencyScatter(corr_tgt_rel_pop, latency_r_nb_pop, param, selectedIDs, animalid_pop)
%single-trial latency correlation 2024June
%
%  latency_r_nb_pop: neuro-bhv latency correlation
% latency_p_nb_pop: p-value for neuro-bhv latency correlation
% nLatencyTrials_pref_pop: number of trials in which latency was successfully detected
% param.nLatencyTr_pref_th

nUnits = size(corr_tgt_rel_pop,2);
msize = 7;%10;%30;
edgeColor = [.75 .75 .75];%'none';
figPosition = [0 0 400 200];
alpha = 1;
scatterLimit = [0 1];

if nargin < 7
    tgtModalities = [1 2];
end

if nargin < 6
    animalid_pop = ones(1, nUnits);
end

if nargin < 5
    selectedIDs = [];
end

f = figure('position', figPosition); %[0 0 500 1000]


%latencyOK = (sum(isnan(stratified_latency_pop),1) < 1);
units_latency = nLatencyTrials_pref_pop >= param.nLatencyTr_pref_th;
units_selectedIDs = zeros(1, nUnits);
units_selectedIDs(selectedIDs) = 1;

for aa = 1:numel(unique(animalid_pop))
    switch aa
        case 1
            asymbol = 'o';% acolor = [.5 .5 .5];%'b';
        case 2
            asymbol = 'diamond'; %acolor = [.5 .5 .5];%'r';
    end

    ax = subplot(121);
    % all units
    % s0 = scatter(corr_tgt_rel_pop(tgtModalities(2), animalid_pop==aa), ...
    %     corr_tgt_rel_pop(tgtModalities(1), animalid_pop==aa), ...
    %     msize, [.5 .5 .5], asymbol);
    % s0.MarkerEdgeColor = acolor;
    % hold on
    
    % units with sufficient number of detected latency trials
    s1 = scatter(corr_tgt_rel_pop(tgtModalities(2), units_latency.*animalid_pop==aa), ...
        corr_tgt_rel_pop(tgtModalities(1), units_latency.*animalid_pop==aa),...
        msize, latency_r_nb_pop(units_latency.*animalid_pop==aa), 'filled', asymbol, 'linewidth', .25);
    s1.MarkerEdgeColor = edgeColor;
    s1.MarkerFaceAlpha = alpha;
    hold on;

    % units used as examples
    s2 = scatter(corr_tgt_rel_pop(tgtModalities(2), units_selectedIDs.*animalid_pop==aa), ...
        corr_tgt_rel_pop(tgtModalities(1), units_selectedIDs.*animalid_pop==aa), ...
        msize, [0 0 0], asymbol, 'linewidth', .5);%stratified_avgCorr_pop(selectedIDs));


    if aa == 2
        clim(scatterLimit);
        colormap(flipud(slanCM('plasma'))); 
        xlabel(param.predictorNames(tgtModalities(2)));
        ylabel(param.predictorNames(tgtModalities(1)));
        %title(sprintf('corr to tgt, relative to full mdl\nlatency correlation (success pref)'));
        squareplot(gca, [-50 120]);
        line([120 -50], [-50 120],'color', 'r')
        [h,g] = mcolorbar(gca, 0.5);
        g.Label.String = 'r';
        title(['n=' num2str(sum(units_latency))]);

        set(gca,'tickdir','out');
    end

      switch tgtModalities(1)
            case 1 
                ax.YColor = 'r';
            case 2
                ax.YColor = 'b';
            case 3
                ax.YColor = 'g';
        end
        switch tgtModalities(2)
            case 1 
                ax.XColor = 'r';
            case 2
                ax.XColor = 'b';
            case 3
                ax.XColor = 'g';
        end


    if aa==2 %use both animals for stats
        thisAxis = corr_tgt_rel_pop(tgtModalities(2), units_latency) ...
            - corr_tgt_rel_pop(tgtModalities(1), units_latency);
        significant = latency_r_nb_pop(units_latency)  > param.r_latency_th;
        
        ax(1)= subplot(122);
        histogram(thisAxis(significant == 0), -170:20:170, 'facecolor', "#FFFF00");%yellow
        set(gca,'tickdir','out');
        hold on;
        histogram(thisAxis(significant == 1), -170:20:170, 'facecolor', "#7E2F8E"); %purple
        axis square; box off;
        legend(['r<' num2str(param.r_latency_th)], ['r>' num2str(param.r_latency_th)],...
            'Location','northwest');

        %p_skew_all =  signrank(thisAxis);
        %p_skew_selected =  signrank(thisAxis(significant==1));
        %p = ranksum(thisAxis(significant==0), thisAxis(significant==1));
        p = mediantest(thisAxis(significant==0), thisAxis(significant==1));
        %[~, p] = adtest2(thisAxis(significant==0)', thisAxis(significant==1)');
        % [p, tbl] = anova1(thisAxis, significant)
        % p = kruskalwallis(thisAxis', significant');


        title(ax(1),['p: ' num2str(p)]);
        xlabel('x - y');
        set(gca,'tickdir','out', 'xcolor','r');
        linkaxes(ax);
        vline(0,gca);
    end

end


