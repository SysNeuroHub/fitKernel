% https://towardsdatascience.com/detecting-stationarity-in-time-series-data-d29e0a21e638
% https://en.wikipedia.org/wiki/Structural_break
%
%         load(loadNames{idata},'ephysdata','dd');
%         spk_all = ephysdata.spikes.spk;
%         [spk_all_cat, t_cat] = concatenate_spk(spk_all, {dd.eye.t});
%         
%         dt_r = median(diff(t_r));
%         
%         PSTH_r = getPSTH(spk_all_cat, t_r);
%         PSTH_f = filtPSTH(PSTH_r, dt_r, param.psth_sigma, 2);
% 
%         %Cusum test for structural change
%         h=cusumtest((1:numel(PSTH_f))', PSTH_f, plot='on',display='summary')
%         [,p]=pptest(PSTH_f)
%         [,p]=vratiotest(PSTH_f)
%
% other candidates:
% chowtest: need to supply a break point
% https://github.com/zhaokg/Rbeast
