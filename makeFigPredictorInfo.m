animal = 'hugo'; %'m1899' 'andy' 'ollie' 
year = '2021';

saveServer = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\Daisuke\cuesaccade_data';%'E:/tmp/cuesaccade_data';

load(fullfile(saveServer,year,animal,'predictorInfo_hugo12December_14.mat'));
load(fullfile(saveServer,'param20230405.mat'),'param');
load('\\storage.erc.monash.edu.au\shares\R-MNHS-Physio\SysNeuroData\Monash Data\Joanita\2021\cuesaccade_data\12December\14\saved_oephysdata\hugo_oephysdata_ch1.mat','dd')
load(fullfile(saveServer,year,animal,'eyeCat_hugo12December_14.mat'))

tlarge = [560 590];
fig = showPredictorInfo(eyeData_rmotl_cat,predictorInfo, param, catEvTimes, dd,tlarge);
xlim(tlarge);

saveas(gcf, 'FigPredictorInfo.fig');
savePaperFigure(gcf, 'FigPredictorInfo');

load('Z:\Shared\Daisuke\cuesaccade_data\2021\hugo\hugo12December_14_1.mat');
plot(t_r, PSTH_f, 'k');hold on; plot(t_r, predicted_all, 'color',[1 .8 0], 'linewidth',2);
xlim(tlarge);
savePaperFigure(gcf, 'FigPredictorInfo_spk');

% tsmall=[571.1 573.7];%[566.5 569.5];
% fig = showPredictorInfo(eyeData_rmotl_cat,predictorInfo, param, catEvTimes, dd,tsmall);
% xlim(tsmall);
% savePaperFigure(gcf, 'FigPredictorInfo_enlarge');



