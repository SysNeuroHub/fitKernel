function [f, depth_pop] = function_fig_supplement2_3(corr_tgt_avg_pop, ...
    animalid_pop, animals, ppcdays, ppcdepths, id_pop, param)

depth_pop = nan(size(id_pop));
for ianimal = 1:numel(animals)
    thesedays = ppcdays.(animals{ianimal});
    thesedepths = ppcdepths.(animals{ianimal});
    theseids_depth = cell(size(thesedepths));
    selectedUnit = zeros(size(thesedepths,1)*size(thesedepths,2),1);
    for isession = 1:numel(thesedays)

        % Convert to datetime
        dt = datetime(thesedays{isession}, 'InputFormat', 'yyyy/MM/dd/');

        % Create reformatted string: '2021/03March/16/'
        newDateStr = sprintf('%04d/%02d%s/%02d/', ...
            year(dt), ...
            month(dt), ...
            char(month(dt,'name')), ...
            day(dt));
        for iunit = 1:size(thesedepths,1)
            theseids_depth{iunit, isession} = sprintf('%s/%s%d',animals{ianimal}, newDateStr, iunit);
        end
    end
    theseids_depth_1d = theseids_depth(:);
    thesedepths_1d = thesedepths(:);

    [~,idx_pop, idx_depth]=intersect(id_pop, theseids_depth_1d);
    depth_pop(idx_pop) = thesedepths_1d(idx_depth);
    selectedUnit(idx_depth) = 1;
    selectedUnit = reshape(selectedUnit, size(thesedepths,1), size(thesedepths,2));
end

f = showDepthScatter(corr_tgt_avg_pop(2:4,:), depth_pop, [], animalid_pop, param);
