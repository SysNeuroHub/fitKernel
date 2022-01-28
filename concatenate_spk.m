function [spk_cat, t_cat] = concatenate_spk(spk_tr, t_tr)
%  [spk_cat, t_cat] = concatenate_spk(spk_tr, t_tr)

%t_tr={dd.eye.t};
dt = median(diff(t_tr{1}));
spk_cat = [];
t_cat = [];
for itr = 1:length(spk_tr)
    if isempty(t_cat)
        t0 = t_tr{1}(1);
    else
        t0 =  max(t_cat)-t_tr{itr}(1)+dt;
    end
    t_cat = cat(1, t_cat, t_tr{itr}+t0);
    
    spk_cat = cat(1, spk_cat, spk_tr{itr}+t0);
end