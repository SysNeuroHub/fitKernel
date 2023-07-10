tgtKernel = 'eyespeed';
tgtidx = find(strcmp(param.predictorNames, tgtKernel))+1;

c1=expval_ind_pop(1,:)>3;
others = setxor(2:4, tgtidx);
c2=expval_tgt_pop(tgtidx,:)./expval_tgt_pop(others(1),:)>2;
c3=expval_tgt_pop(tgtidx,:)./expval_tgt_pop(others(2),:)>2;

candidates=find(c1.*c2.*c3==1);
[sorted, sortedIDc]=sort(expval_tgt_pop(3,candidates)./expval_tgt_pop(1,candidates),'descend');
sortedID = candidates(sortedIDc);

thisIdx = sortedID(11:20);
thisID = id_pop(thisIdx)'

expval_ind_pop(:,thisIdx)'
expval_tgt_pop(:,thisIdx)'


%new inclusion criteria?:
% expval_ind_pop(1,:)>3;
