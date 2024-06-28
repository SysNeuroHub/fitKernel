function [okunits, mfiringRateOK, expvalOK, ntargetTrOK, ptonsetRespOK] ...
    = inclusionCriteria(mFiringRate_pop, expval_pop, ntargetTrials_pop, PtonsetResp_pop, param)
%[okunits, mfiringRateOK, expvalOK, ntargetTrOK, ptonsetRespOK] ...
% = inclusionCriteria(mFiringRate_pop, expval_pop, ntargetTrials_pop, PtonsetResp_pop, param)
%
% INPUT:
% mFiringRate_pop: 
% expval_pop: 
% ntargetTrials_pop:
% PtonsetResp_pop:
%
% param requires following fields:
%   mfiringRateTh
%   expvalTh
%   ntargetTrTh
%   ptonsetRespTh

if size(mFiringRate_pop(1)) < size(mFiringRate_pop(2))
    mFiringRate_pop = mFiringRate_pop';
end
if size(expval_pop(1)) < size(expval_pop(2))
    expval_pop = expval_pop';
end
if size(ntargetTrials_pop(1)) < size(ntargetTrials_pop(2))
    ntargetTrials_pop = ntargetTrials_pop';
end
if size(PtonsetResp_pop(1)) < size(PtonsetResp_pop(2))
    PtonsetResp_pop = PtonsetResp_pop';
end

mfiringRateOK = (mFiringRate_pop > param.mfiringRateTh);
expvalOK = (expval_pop>param.expvalTh);
ntargetTrOK = (ntargetTrials_pop>param.ntargetTrTh );
ptonsetRespOK = (PtonsetResp_pop < param.ptonsetRespTh);

okunits = find(mfiringRateOK.*expvalOK.*ntargetTrOK.*ptonsetRespOK);
