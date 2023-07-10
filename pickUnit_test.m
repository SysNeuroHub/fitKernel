data = 100*expval_tgt_pop(2:4,:)./expval_tgt_pop(1,:);
plot(data(1,:),data(3,:), '.');
xlim([-50 200]);ylim([-50 200]);
margin = 2;
% xtgt = -13;
% ytgt = 118;

gg=ginput(1);
xtgt = gg(1);ytgt=gg(2);

thisIdx=intersect(find(abs(data(1,:)-xtgt)<margin),find(abs(data(3,:)-ytgt)<margin))
thisID = id_pop(thisIdx)

expval_ind_pop(:,thisIdx)'
expval_tgt_pop(:,thisIdx)'
