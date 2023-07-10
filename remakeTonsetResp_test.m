rootFolder = '//storage.erc.monash.edu.au/shares/R-MNHS-Physio/SysNeuroData/Monash Data/Joanita/';
%saveServer = 'E:/tmp/cuesaccade_data';
saveServer = 'Z:\Shared\Daisuke\cuesaccade_data';

%load('fitPSTH_pop20230704hugo.mat')

theseIDs = {'hugo/2021/09September/01/25',...
    'hugo/2022/03March/10/20',...
    'hugo/2022/07July/29/19'};


for ii = 1:numel(theseIDs)
    
    aaa = split(theseIDs{ii},'/');
    animal = aaa{1};
    year = aaa{2};
    months = aaa{3};
    dates = aaa{4};
    channels = aaa{5};

    saveFigFolder = fullfile(saveServer, '20230704',year,animal);
    mkdir(saveFigFolder);
    
     datech = [months filesep dates filesep channels];
     saveSuffix = [animal replace(datech,filesep,'_') ];%'_cue'];

    saveFolder = fullfile(saveServer, year,animal);%17/6/23
    saveName = fullfile(saveFolder, [saveSuffix '.mat']);
    
    load(saveName, 'PSTH_f','predicted_all', 'predicted','kernelInfo'...
        ,'t_r','param','t_cat');%,'dds');
    psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);
    
    eyeName = fullfile(saveFolder,['eyeCat_' animal months '_' dates '.mat']);
     load(eyeName,'catEvTimes','startSaccNoTask','saccDirNoTask');

     loadName = getCuesaccadeName(rootFolder, theseIDs{ii});
     load(loadName,'dd');
     
    y_r = cat(2,PSTH_f,predicted_all, predicted);
    
    if min(PSTH_f)<0
        k = -min(PSTH_f)+0.5;
    else
        k = 0.5;
    end
    
    [f] = showTonsetResp(t_r, (y_r), catEvTimes, dd, psthNames, ...
        startSaccNoTask, saccDirNoTask, param, [-0.5 0.5], 1);
    screen2png(fullfile(saveFigFolder,['cellclassFig_' saveSuffix '_allTr']), f);
    close(f);

    [f] = showTonsetResp(t_r, y_r, catEvTimes, dd, psthNames, ...
        startSaccNoTask, saccDirNoTask, param, [-0.5 0.5], 0);
    screen2png(fullfile(saveFigFolder,['cellclassFig_' saveSuffix]), f);
    close(f);
    
    %% traces in log space
    % problem of log(0) was dealt with pseudolog transformation
    [f] = showTonsetResp(t_r, log(y_r+k), catEvTimes, dd, psthNames, ...
        startSaccNoTask, saccDirNoTask, param, [-0.5 0.5], 0);
    screen2png(fullfile(saveFigFolder,['cellclassFig_' saveSuffix '_log']), f);
    close(f);
    
    [f] = showTonsetResp(t_r, log(y_r+k), catEvTimes, dd, psthNames, ...
        startSaccNoTask, saccDirNoTask, param, [-0.5 0.5], 1);
    screen2png(fullfile(saveFigFolder,['cellclassFig_' saveSuffix '_allTr_log']), f);
    close(f);

end
