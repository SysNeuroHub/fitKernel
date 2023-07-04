function [trIdx_r, trIdx] = retrieveTrIdx_r(t_cat, t_cat_r, t_tr)
%[trIdx_r, trIdx] = retrieveTrIdx_r(t_cat, t_cat_r, t_tr)
%
% INPUT:
%    t_cat: consecutive time (product of concatenate_spk.m)
%    t_cat_r: consecutive time, resampled ( = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))')
%    t_tr: trial time (product of processEyeData = {eyeData.t})
% OUTPUT:
%    trIdx_r: trial time, resampled
%    trIdx: trial index
%
% 4/7/2023 overhaul so that trIdx_r have no gap between trials

catidx = 1;
trIdx = [];
trIdx_re = [];
for itr = 1:numel(t_tr)
    trIdx{itr} = catidx:catidx+numel(t_tr{itr})-1;
    catidx = catidx + numel(trIdx{itr});
    
    endTime = max(t_cat(trIdx{itr}));
    
    trIdx_re(itr) = find(t_cat_r<=endTime, 1, 'last' );
end

trIdx_r{1} = (1:trIdx_re(1))';
for itr = 2:numel(t_tr)
    trIdx_r{itr} = (trIdx_re(itr-1)+1:trIdx_re(itr))';
end

% %% sanity check
% tidx_r = [];
% for itr = 1:numel(trIdx_r)
%     tidx_r = [tidx_r; trIdx_r{itr}];
% end

