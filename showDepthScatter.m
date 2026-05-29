function f = showDepthScatter(corr_tgt_rel_pop,  depth_pop,  selectedIDs,  animalid_pop, param)
%f = showHMScatter(corr_tgt_rel_pop, p_hm_pop, animalid_pop)


showRange = [-0.4 1];%[0 1]; %[-50 120]
units_depthDetected = logical(~isnan(depth_pop).*(depth_pop>0));
msize = 7;%10;
figPosition = [0 0 600 400];
edgeColor = [.75 .75 .75];%'none';
alpha = 1;
scatterLimit = [0 3900];

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

    for imodality = 1:2
        switch imodality
            case 1
                tgtModalities = [1 2];
                %1: vision, 2: eye spd, 3: eye position
            case 2
                tgtModalities = [2 3];
        end

        ax(2*imodality-1) = subplot(2,2,2*imodality-1);
        s0.MarkerEdgeColor = acolor;
        hold on

        if aa==1
            line(fliplr(showRange), showRange,'linestyle',':','color','k','linewidth',0.25)
        end
        %% alll units
        s1 = scatter(corr_tgt_rel_pop(tgtModalities(2), units_depthDetected.*animalid_pop==aa), ...
            corr_tgt_rel_pop(tgtModalities(1), units_depthDetected.*animalid_pop==aa), ...
            msize, depth_pop(units_depthDetected.*animalid_pop==aa), 'filled', asymbol);
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
        %colormap(RedWhiteBlue);

        xlabel(param.predictorNames(tgtModalities(2)));
        ylabel(param.predictorNames(tgtModalities(1)));
        title(['n=' num2str(sum(units_depthDetected))]);

        switch tgtModalities(1)
            case 1 
                ax(2*imodality-1).YColor = 'r';
            case 2
                ax(2*imodality-1).YColor = 'g';
            case 3
                ax(2*imodality-1).YColor = 'b';
        end
        switch tgtModalities(2)
            case 1 
                ax(2*imodality-1).XColor = 'r';
            case 2
                ax(2*imodality-1).XColor = 'g';
            case 3
                ax(2*imodality-1).XColor = 'b';
        end
      
        if aa==2
            squareplots(gca, showRange);
            [hh,gg] = mcolorbar(gca, 0.5);
            gg.Label.String = 'Depth [um]';
        end

        if aa==2 %use both animals for stats
            %% 
            ax(2*imodality) = subplot(2,2,2*imodality);

            thisAxis = corr_tgt_rel_pop(tgtModalities(2), units_depthDetected) ...
                - corr_tgt_rel_pop(tgtModalities(1), units_depthDetected);

            plot(thisAxis, depth_pop(units_depthDetected), '.');
            hold on
            [rho, p_corr] = corr(thisAxis', depth_pop(units_depthDetected)','type','Spearman');
            title(sprintf('rho=%.1e, p=%.1e', rho, p_corr));

            xlim([-1 1]);
            ylim(scatterLimit);
            xlabel([param.predictorNames{tgtModalities(2)} ' - ' param.predictorNames{tgtModalities(1)}]);
            
            %% stats for each animal
            for ianimal = 1:2
                thisAxis = corr_tgt_rel_pop(tgtModalities(2), animalid_pop==ianimal.*units_depthDetected) ...
                    - corr_tgt_rel_pop(tgtModalities(1), animalid_pop==ianimal.*units_depthDetected);
                %significant = (abs(auc_hm_pop(iregress, animalid_pop==ianimal.*units_depthDetected)) > param.auc_th);
                %p_skew_animal(iregress,ianimal) =  ranksum(thisAxis(significant==1), thisAxis(significant==0)); % whetehr the value changes
                [rho, p_corr] = corr(thisAxis', depth_pop(animalid_pop==ianimal.*units_depthDetected)','type','Spearman');
                disp(['tgt modality: ' num2str(imodality) ', animal: ' num2str(ianimal) ',rho: ' num2str(rho) 'p: ' num2str(p_corr)]);
            end
        end
    end
end
%linkaxes(ax([2 4]));
%vline(0,ax([2 4]),[],[],1);

