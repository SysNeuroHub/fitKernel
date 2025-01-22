function f = showHMScatter(corr_tgt_rel_pop, p_hm_pop, selectedIDs, animalid_pop)
%f = showHMScatter(corr_tgt_rel_pop, p_hm_pop, animalid_pop)


units_hm = ~isnan(p_hm_pop(1,:));
msize = 30;

nUnits = size(corr_tgt_rel_pop,2);

if nargin < 3
    animalid_pop = ones(1, nUnits);
end

units_selectedIDs = zeros(1, nUnits);
units_selectedIDs(selectedIDs) = 1;

f = figure('position',[0 0 1000 2000]);

for aa = 1:numel(unique(animalid_pop))
    switch aa
        case 1
            asymbol = 'o'; acolor = [.5 .5 .5];%'b';
        case 2
            asymbol = 'diamond'; acolor = [.5 .5 .5]; %'r';
    end

    for iregress = 1:2
        subplot(2,1,iregress);
        s0.MarkerEdgeColor = acolor;
        hold on
        
        s1 = scatter(corr_tgt_rel_pop(2, units_hm.*animalid_pop==aa),corr_tgt_rel_pop(1, units_hm.*animalid_pop==aa), ...
            msize, -log(p_hm_pop(iregress, units_hm.*animalid_pop==aa)), 'filled', asymbol);
        s1.MarkerEdgeColor = acolor;
         s1.MarkerFaceAlpha = .5;

        s2 = scatter(corr_tgt_rel_pop(2, animalid_pop==aa.*units_selectedIDs), ...
            corr_tgt_rel_pop(1, animalid_pop==aa.*units_selectedIDs), msize, [0 0 0], asymbol,'LineWidth',2);
        s2.LineWidth = 2;

        clim([0 4]);
        colormap("cool");
        %colormap(flipud(summer));
        hh=colorbar; hh.Label.String = '-log(p hit v miss)';
        xlabel('eye speed'); ylabel('stimulus');
        tname = sprintf('Hit v Miss at preferred direction');
        if iregress == 2
            tname = [tname ' after regression'];
        end
        title(tname);
    end
end
squareplots(f, [-50 120]);
