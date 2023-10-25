[saveServer, rootFolder] = getReady();
saveServer = 'Z:\Shared\Daisuke\cuesaccade_data';
animal = 'hugo';% 'andy' 'ollie'
year = '2022';
figTWin = [-0.5 0.5];
onlySuccess = 0;
respWin = [0.05 0.35]; %[s]

load(fullfile(saveServer,'param20230405.mat'),'param');


saveFigFolder = fullfile(saveServer, '20230713',year,animal);
mkdir(saveFigFolder);

[loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);

 thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
        regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','08August','05','*_ch2*')))));

ng =[];
for idata = 1:length(channels)
    try        
        
        datech = [months{idata} filesep dates{idata} filesep num2str(channels{idata})];
        disp(datech);
        
        saveSuffix = [animal replace(datech,filesep,'_') '_linear_rReg'];%'_cue'];
        
        thisDate = [months{idata} '_' dates{idata}];
        saveFolder = fullfile(saveServer, year,animal);%17/6/23
        
        %eyeName = 'Z:\Shared\Daisuke\cuesaccade_data\2021\hugo\eyeCat_hugo08August_13.mat';
        eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
        load(eyeName,'catEvTimes');
        %saveName = 'Z:\Shared\Daisuke\cuesaccade_data\2021\hugo\hugo08August_13_21_linear_rReg.mat';
        saveName = fullfile(saveFolder, [saveSuffix '.mat']);
        load(saveName,'PSTH_f','predicted_all','t_r','dd');
        
        y_r = cat(2,PSTH_f,predicted_all);
        
        %% obtain gain
        gainInfo = getGainInfo(t_r, y_r(:,1:2), param.cardinalDir, catEvTimes, ...
            dd, figTWin, onlySuccess, respWin);
        f=showGainInfo(gainInfo);
        savefigname = fullfile(saveFigFolder,[saveSuffix '_gainInfo']);
        screen2png(savefigname);
        close(f);
        save(saveName, 'gainInfo','-append');
        clear gainInfo y_r catEvTimes PSTH_f
        
    catch err
        disp(err);
        ng = [ng idata];
        close all;
        
    end
end
    
