function [okunits, mfiringRateOK, expvalOK, ntargetTrOK, ptonsetRespOK] ...
    = inclusionCriteria(mFiringRate_pop, expval_pop, ntargetTrials_pop, PtonsetResp_pop, param)
mfiringRateOK = (mFiringRate_pop > param.mfiringRateTh);
expvalOK = (expval_pop>param.expvalTh);
ntargetTrOK = (ntargetTrials_pop>param.ntargetTrTh );
ptonsetRespOK = (PtonsetResp_pop < param.ptonsetRespTh);

okunits = find(mfiringRateOK+expvalOK+ntargetTrOK+ptonsetRespOK ==4);
