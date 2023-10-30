Rsqadj_pop_a = Rsqadj_pop([2 3 4],:);
Rsqadj_pop_r = Rsqadj_pop([2 3 4],:)./Rsqadj_pop([1],:);

idx_v = find(Rsqadj_pop_r(1,:) < 0.7);
idx_es = find(Rsqadj_pop_r(2,:) < 0.5);
idx_ep = find(Rsqadj_pop_r(3,:) < 0.5);


id_pop([idx_v idx_es idx_ep])'

