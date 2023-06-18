function [spk_cat, t_cat] = concatenate_spk(spk_tr, eyeData)
%  [spk_cat, t_cat] = concatenate_spk(spk_tr, eyeData)
% 18/6/23 refactored so t_cat compuation is identical to concatenate_eye.m

%t_tr={dd.eye.t};
%dt = median(diff(t_tr{1}));
spk_cat = [];
t_cat = [];
for itr = 1:length(spk_tr)
    if isempty(t_cat)
       t0 = eyeData(itr).t(1);
 elseif ~isempty(eyeData(itr).t)
        t0 =  max(t_cat) -eyeData(itr).t(1)+eyeData(itr).dt;
    end

    if ~isempty(eyeData(itr).t)
        t_cat = cat(1, t_cat, eyeData(itr).t+t0);
        spk_cat = cat(1, spk_cat, spk_tr{itr}+t0);
    end
end
t_cat = t_cat(~isnan(t_cat));
[t_cat, ix] = unique(t_cat);
