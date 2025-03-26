
for aa =1:2
    [okunits, corr_tgtOK, ntargetTrOK, ptonsetRespOK] ...
        = inclusionCriteria(corr_tgt_pop(1,animalid_pop==aa), ntargetTrials_pop(animalid_pop==aa),...
        PtonsetResp_pop(animalid_pop==aa), param);
    
    switch aa
        case 1
            animal =  'hugo';%'ollie';% % %'andy';%
        case 2
            animal =  'ollie';% % %'andy';%
    end

    disp(['In animal ' animal])
    disp(['Criteria2 (firing rate) ' num2str(sum(animalid_pop==aa))]);
    disp('amongst neurons that passed this criterion, ')
    disp(['Criteria 1 (#trials) ' num2str(sum(ntargetTrOK))])
    disp(['Criteria 3 (direction sensitivity) ' num2str(sum(ptonsetRespOK))])
    disp(['Criteria 4 (fit to model) ' num2str(sum(corr_tgtOK))])
    end
