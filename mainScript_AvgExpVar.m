%% TO BE DELETED

set(0,'DefaultFigureVisible','off');
[saveServer, rootFolder] = getReady();
%saveServer = 'Z:\Shared\Daisuke\cuesaccade_data';
warning('off','all');

%% recorded data

for aaa = 1%:2
    switch aaa
        case 1
            animal = 'hugo'; %'ollie';%%'andy';% 'andy' '
            yidx=1%:3;
        case 2
            animal = 'ollie';
            yidx=3;
    end

    for yyy = yidx
        switch yyy
            case 1
                year = '2021';
            case 2
                year = '2022';
            case 3
                year = '2023';
        end
        logName = fullfile(saveServer,'20250207',year,animal,'log_mainScript_latency');

        saveFigFolder = fullfile(saveServer, '20250207',year,animal);
        mkdir(saveFigFolder);


        [loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);

        % thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
        %     regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','07July','26','*_ch19'))))); %2022
        % thisdata = [thisdata find(1-cellfun(@isempty, regexp(loadNames, ...
        %     regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','08August','15','*_ch4')))))]; %2022
        thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
            regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','09September','17','*_ch32'))))); %2021
        % thisdata = 1:length(channels);


        % parameters
        n=load(fullfile(saveServer,'param20250207.mat'),'param');
        param =n.param;
        n=[];
        ncDirs = length(param.cardinalDir);

        psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);

        ng = [];
        previousDate = [];
        for idata = thisdata
           try
                % datech = [years{idata} filesep months{idata} filesep dates{idata} filesep num2str(channels{idata})];
                datech = [months{idata} filesep dates{idata} filesep num2str(channels{idata})];
                disp([num2str(idata) '/' num2str(numel(thisdata)) ', ' datech ]);

                saveSuffix = [animal replace(datech,filesep,'_') '_linear_rReg'];%'_cue'];

                thisDate = [months{idata} '_' dates{idata}];

                saveFolder = fullfile(saveServer, year,animal);%17/6/23
                saveName = fullfile(saveFolder, [saveSuffix '.mat']);

                load(saveName,'mFiringRate');
                if ~exist("mFiringRate",'var') || mFiringRate < 5
                    clear mFiringRate
                    continue;
                end
                load(saveName, 't_r','t_cat','spkOkUCueTrials','PSTH_f','predicted_all','predicted','spkOkTrials');

                EE = load(loadNames{idata},'dd');
                dd = EE.dd;
                clear EE

                y_r = cat(2,PSTH_f,predicted_all, predicted);

                eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
                load(eyeName,'catEvTimes');

                [corr_tgt_avg, corr_tgt_avg_rel] = getCorr_tgt_avg(t_r, y_r, catEvTimes, dd, param, spkOkUCueTrials);

                %% event triggered traces
                load(eyeName,'eyeData_rmotl_cat', 'blinks', 'outliers', 't_tr', 'onsets_cat');
                [startSaccNoTask_spkOkUCue, endSaccNoTask_spkOkUCue, saccDirNoTask_spkOkUCue, dirIndexNoTask_spkOkUCue] = ...
                    getSaccNoTask(t_cat, catEvTimes, eyeData_rmotl_cat, blinks, outliers, t_tr, onsets_cat, spkOkTrials, param);

                [f, cellclassInfo] = showTonsetResp(t_r, y_r, catEvTimes, dd, psthNames, ...
                    startSaccNoTask_spkOkUCue, saccDirNoTask_spkOkUCue, param, [], spkOkUCueTrials);
                cellclassInfo.datech = datech;
                savePaperFigure(f, fullfile(saveFigFolder,['cellclassFig_' saveSuffix]));
                close all

                %% save results
                save(saveName,'corr_tgt_avg','corr_tgt_avg_rel', '-append');

                clear dd mFiringRate t_r f
            catch err
                disp(err);
                ng = [ng idata];
                clear kernelInfo kernelInfo_norm mFiringRate t_r 
                close all;
            end
            save(logName, "ng",'idata');
        end
    end
end
