%% make eyedat concatenated across trials from scratch
% register FAIL trials with getChoice.m called in concatenate_eye.m
rootFolder = '//storage.erc.monash.edu.au/shares/R-MNHS-Physio/SysNeuroData/Monash Data/Joanita/2021/cuesaccade_data/';
saveFolder = 'E:/tmp/cuesaccade_data';

animal = 'hugo';

[rawDataNames, months, dates, channels] = getMonthDateCh(animal, rootFolder);

for idata = 1:length(months)
    monthDate{idata} = [months{idata} '_' dates{idata}];
end
[monthDate, idx] = unique(monthDate);

load('E:\tmp\cuesaccade_data\param20220626','param');

tic
for ii = 84%:numel(monthDate)
    
    disp(monthDate{ii});
    
    saveName = fullfile(saveFolder,['eyeCat_' animal monthDate{ii} '.mat']);
    load(rawDataNames{idx(ii)},'dd');
    
    [eyeData_rmotl_cat, catEvTimes, t_tr, onsets_cat, meta_cat, blinks, outliers] ...
        = processEyeData(dd.eye, dd, param);
    close all;
    save(saveName,'eyeData_rmotl_cat','catEvTimes',...
        'onsets_cat','meta_cat','blinks','outliers','t_tr');
    
    
    %% prepare predictor variables
    if max(diff(eyeData_rmotl_cat.t)) > 0.02
        %< jump in eyeData_rmotl_cat.t
        %'03March_24' '06June_09' '06June_11'
        disp(['NG ' monthDate{ii}]);
        continue;
    end
    t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
    predictorInfo = preparePredictors(dd, eyeData_rmotl_cat, t_r, param, catEvTimes);
    save(fullfile(saveFolder,['predictorInfo_' animal monthDate{ii} '.mat']), ...
        'predictorInfo','t_r');

end
t= toc