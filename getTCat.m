function [t_cat] = getTCat(t_tr)
%  created from concatenate_spk(spk_tr, t_tr)

%t_tr={dd.eye.t};
dt = median(diff(t_tr{1}));
t_cat = [];
for itr = 1:length(t_tr)
    if isempty(t_cat)
        t0 = t_tr{1}(1);
    else
        t0 =  max(t_cat)-t_tr{itr}(1)+dt;
    end
    t_cat = cat(1, t_cat, t_tr{itr}+t0);
end