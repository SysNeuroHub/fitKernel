

%% recorded data
animal = 'hugo';
rootFolder = '\\storage.erc.monash.edu.au\shares\R-MNHS-Physio\SysNeuroData\Monash Data\Joanita\2021/cuesaccade_data/';
saveFolder = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\Daisuke\cuesaccade_data';

cutoffFreqs = [1e-3 1e-2 1e-1 1e0];

[loadNames, months, dates, channels] = getMonthDateCh(animal, rootFolder);

%idata = find(1-cellfun(@isempty, regexp(loadNames, regexptranslate('wildcard','06June\18\*_ch17.mat'))));

for idata = 1:length(channels) %1061;%865;%

datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
disp(datech);

saveSuffix = [animal replace(datech,'/','_')];

thisDate = [months{idata} '_' dates{idata}];

    
    saveName = fullfile(saveFolder, [saveSuffix '.mat']);
    
    if ~exist(saveName,'file')
        continue;
    end
    
    load(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');

    load(saveName, 'mFiringRate',...
        'param','eyeData_rmotl_cat');

    load(fullfile(saveFolder,['eyeCat_' thisDate '.mat']), 't_tr');
    load(loadNames{idata},'ephysdata');
   
    
    %% concatenate across trials
    spk_all = ephysdata.spikes.spk;
    clear ephysdata
    [spk_all_cat, t_cat] = concatenate_spk(spk_all, t_tr);
    clear spk_all
    
    
    Rms = corrPSTH_eye_ms(spk_all_cat, eyeData_rmotl_cat, ...
        predictorInfo.predictors_r(17,:)', param.dt_r, cutoffFreqs);
    screen2png(['corrPSTH_eye_ms_' saveSuffix]);
    close;
    
    save(saveName, 'Rms','cutoffFreqs','-append');
end


