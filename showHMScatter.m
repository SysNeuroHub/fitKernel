function f = showHMScatter(corr_tgt_rel_pop, ranksumz_hm_pop, selectedIDs, animalid_pop, tgtModalities, param)
%f = showHMScatter(corr_tgt_rel_pop, p_hm_pop, animalid_pop)

if nargin < 5
    tgtModalities = [1 2];
    %1: vision, 2: eye spd, 3: eye position
end

units_hm = ~isnan(ranksumz_hm_pop(1,:));
msize = 7;%10;
figPosition = [0 0 400 400];
edgeColor = [.75 .75 .75];%'none';
alpha = 1;
scatterLimit = [0 4];

nUnits = size(corr_tgt_rel_pop,2);

if nargin < 3
    animalid_pop = ones(1, nUnits);
end

units_selectedIDs = zeros(1, nUnits);
units_selectedIDs(selectedIDs) = 1;

f = figure('position', figPosition);

for aa = 1:numel(unique(animalid_pop))
    switch aa
        case 1
            asymbol = 'o'; acolor = edgeColor;%[.5 .5 .5];%'b';
        case 2
            asymbol = 'diamond'; acolor = edgeColor;% [.5 .5 .5]; %'r';
    end

    for iregress = 1:2
        ax(2*iregress-1) = subplot(2,2,2*iregress-1);
        s0.MarkerEdgeColor = acolor;
        hold on

        %% alll units
        s1 = scatter(corr_tgt_rel_pop(tgtModalities(2), units_hm.*animalid_pop==aa), ...
            corr_tgt_rel_pop(tgtModalities(1), units_hm.*animalid_pop==aa), ...
            msize, abs(ranksumz_hm_pop(iregress, units_hm.*animalid_pop==aa)), 'filled', asymbol);
        s1.MarkerEdgeColor = acolor;
        s1.MarkerFaceAlpha = alpha;
        s1.LineWidth = .25;

        %% selected units
        s2 = scatter(corr_tgt_rel_pop(tgtModalities(2), animalid_pop==aa.*units_selectedIDs), ...
            corr_tgt_rel_pop(tgtModalities(1), animalid_pop==aa.*units_selectedIDs), msize, [0 0 0], asymbol);%,'LineWidth',2);
        s2.LineWidth = .5;

        clim(scatterLimit);
        %colormap("cool");
        %colormap('winter')
        colormap(flipud(slanCM('plasma')))

        xlabel(param.predictorNames(tgtModalities(2)));
        ylabel(param.predictorNames(tgtModalities(1)));

        switch tgtModalities(1)
            case 1 
                ax(2*iregress-1).YColor = 'r';
            case 2
                ax(2*iregress-1).YColor = 'b';
            case 3
                ax(2*iregress-1).YColor = 'g';
        end
        switch tgtModalities(2)
            case 1 
                ax(2*iregress-1).XColor = 'r';
            case 2
                ax(2*iregress-1).XColor = 'b';
            case 3
                ax(2*iregress-1).XColor = 'g';
        end
        tname = sprintf('Hit v Miss at preferred direction');
        if iregress == 2
            tname = [tname ' after regression'];
        end
        title(tname);
        if aa==2
            squareplots(gca, [-50 120]);
            line([120 -50], [-50 120],'color', 'r')
            [hh,gg] = mcolorbar(gca, 0.5);
            gg.Label.String = '|z|';
        end

        if aa==2 %use both animals for stats
            %% histogram
            ax(2*iregress) = subplot(2,2,2*iregress);

            thisAxis = corr_tgt_rel_pop(tgtModalities(2), units_hm) ...
                - corr_tgt_rel_pop(tgtModalities(1), units_hm);
            significant = (abs(ranksumz_hm_pop(iregress, units_hm)) > param.ranksumz_th);

            histogram(thisAxis(significant == 0), -170:20:170,'facecolor',"#FFFF00"); %yellow
            hold on
            histogram(thisAxis(significant == 1), -170:20:170,'facecolor',"#7E2F8E"); %purple
            axis square; box off;
            
            %p_skew(iregress) =  signrank(thisAxis(significant==1)); %whether the highly dissociable units have non-zero value
            p_skew(iregress) =  mediantest(thisAxis(significant==1), thisAxis(significant==0)); % whetehr the value changes 
            title(['p ' num2str(p_skew(iregress))]);
            xlabel('x - y');
            legend(['|z|<' num2str(param.ranksumz_th)], ['|z|>' num2str(param.ranksumz_th)],...
                'location','northwest');
            set(ax(2*iregress),'tickdir','out', 'xcolor','r');
        end
    end
end
linkaxes(ax([2 4]));
vline(0,ax([2 4]));

