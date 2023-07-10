animal = 'hugo'; %'m1899' 'andy' 'ollie' 
% biased eye position
% year = '2022';
% months = '07July';
% dates = '29';
% channels = 28;%30;
%
% too low firing rate
% year = '2021';
% months = '11November';
% dates = '22';
% channels = 28;
%
% firing rate is tad too low. behavior is good
% year = '2022';
% months = '09September';
% dates = '19';
% channels = 10;
%
% good behavior. fitting is not so great
% year = '2022';
% months = '08August';
% dates = '10';
% channels = 4;
%
%
year = '2022';
months = '08August';
dates = '17';
channels = 3;
 

rootFolder = '//storage.erc.monash.edu.au/shares/R-MNHS-Physio/SysNeuroData/Monash Data/Joanita/';
saveServer = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\Daisuke\cuesaccade_data';%'E:/tmp/cuesaccade_data';
thisDate = [months '_' dates];
datech = [months '/' dates '/' num2str(channels)];
load(fullfile(saveServer,'param20230405.mat'),'param');

saveFolder = fullfile(saveServer, year,animal);%17/6/23
load(fullfile(saveFolder,['predictorInfo_' animal thisDate '.mat']), ...
                    'predictorInfo');
%load('\\storage.erc.monash.edu.au\shares\R-MNHS-Physio\SysNeuroData\Monash Data\Joanita\2021\cuesaccade_data\12December\14\saved_oephysdata\hugo_oephysdata_ch1.mat','dd')
[loadNames] = getMonthDateCh(animal, year, rootFolder);
thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
    regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data',months,dates,['*_ch' num2str(channels) '*'])))));
load(loadNames{thisdata},'dd');

eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
load(eyeName);
saveSuffix = [animal replace(datech,'/','_') ];%'_cue'];

saveName = fullfile(saveFolder, [saveSuffix '.mat']);
load(saveName);

tlarge = [840 937];%[860 930];[584 646];%[528.2 530.3];
fig = showPredictorInfo(eyeData_rmotl_cat,predictorInfo, param, catEvTimes, ...
    dd, PSTH_f, predicted_all,tlarge);

savePaperFigure(gcf, 'FigPredictorInfo');

xlim([900.8 902.6]);
savePaperFigure(gcf, 'FigPredictorInfo_e');

xlim([907.9 909.8])
savePaperFigure(gcf, 'FigPredictorInfo_e2');



