function [trIdx_r, trIdx] = retrieveTrIdx_r(t_cat, t_cat_r, t_tr)

catidx = 1;
trIdx = [];
trIdx_r = [];
for itr = 1:numel(t_tr)
    trIdx{itr} = catidx:catidx+numel(t_tr{itr})-1;
    catidx = catidx + numel(trIdx{itr});
    
    startTime = min(t_cat(trIdx{itr}));
    endTime = max(t_cat(trIdx{itr}));
    
    trIdx_r{itr} = find(t_cat_r>startTime & t_cat_r<=endTime);
end