function [okunits,corr_tgtOK, ntargetTrOK, ptonsetRespOK] ...
    = inclusionCriteria(corr_tgt_pop, ntargetTrials_pop, PtonsetResp_pop, param)
%[okunits, mfiringRateOK, corr_tgtOK, ntargetTrOK, ptonsetRespOK] ...
% = inclusionCriteria(mFiringRate_pop, corr_tgt_pop, ntargetTrials_pop, PtonsetResp_pop, param)
%
% INPUT:
% mFiringRate_pop: 
% corr_tgt_pop: 
% ntargetTrials_pop:
% PtonsetResp_pop:
%
% param requires following fields:
%   mfiringRateTh
%   corr_tgtTh
%   ntargetTrTh
%   ptonsetRespTh


if size(corr_tgt_pop(1)) < size(corr_tgt_pop(2))
    corr_tgt_pop = corr_tgt_pop';
end
if size(ntargetTrials_pop(1)) < size(ntargetTrials_pop(2))
    ntargetTrials_pop = ntargetTrials_pop';
end
if size(PtonsetResp_pop(1)) < size(PtonsetResp_pop(2))
    PtonsetResp_pop = PtonsetResp_pop';
end

corr_tgtOK = (corr_tgt_pop>param.corr_tgtTh);
ntargetTrOK = (ntargetTrials_pop>param.ntargetTrTh );
ptonsetRespOK = (PtonsetResp_pop < param.ptonsetRespTh);

okunits = find(corr_tgtOK.*ntargetTrOK.*ptonsetRespOK);
