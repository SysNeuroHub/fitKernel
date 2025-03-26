%% TOBE DELETED

set(0,'DefaultFigureVisible','off');
[saveServer, rootFolder] = getReady();
%saveServer = 'Z:\Shared\Daisuke\cuesaccade_data';
warning('off','all');

%% recorded data

for aaa = 1:2
    switch aaa
        case 1
            animal = 'hugo'; %'ollie';%%'andy';% 'andy' '
        case 2
            animal = 'ollie';
    end

    for yyy = 1:3
        switch yyy
            case 1
                year = '2021';
            case 2
                year = '2022';
            case 3
                year = '2023';
        end
        logName = fullfile(saveServer,'20241212',year,animal,'log_mainScript_latency');

        saveFigFolder = fullfile(saveServer, '20241212',year,animal);
        mkdir(saveFigFolder);


        [loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);

        thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
              regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','03March','18','*_ch13'))))); %2023
        % thisdata = 1:length(channels);


        % parameters
        n=load(fullfile(saveServer,'param20241223.mat'),'param');
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
                load(saveName,'PSTH_f','t_r','spkOkUCueTrials');

                load(loadNames{idata}, 'dd');

                %% prepare behavioral data (common across channels per day)
                eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);

                load(eyeName,'catEvTimes');%,...
      

                %% from showTonsetResp.m
                % figTWin = param.figTWin; %figure temporal window
                % compTWin = [figTWin(1) - 0.2 figTWin(2)]; %triggered response is computed
                %
                % [choiceOutcome] = getChoiceOutcome(dd);
                %
                % successEvents = intersect(find((choiceOutcome==1)), spkOkUCueTrials);
                % successEvents = intersect(successEvents, find(catEvTimes.tOnset + compTWin(2) < max(t_r)));
                %
                % onsetTimes = catEvTimes.tOnset(successEvents);
                % tgtDir = getTgtDir(dd.targetloc(successEvents), param.cardinalDir);
                %
                % [~,dirIdx] = intersect(param.cardinalDir, unique(tgtDir));
                %
                % [~, winSamps, singleOnsetResp, ...
                %     sortedOnsetLabels, uniqueOnsetLabels] ...
                %     = eventLockedAvg(y_r', t_r, onsetTimes, tgtDir, compTWin);
                % tonsetRespAmp = characteriseResp(singleOnsetResp, ...
                %     winSamps, param.tOnRespWin, param.baseWin, 'mean');

                % prefDir = getPrefDir(y_r(:,1), t_r, onsetTimes, tgtDir, param);
                param.cardinalDir = 0:359;
                prefDir = getPrefDir_wrapper(PSTH_f, t_r, dd, catEvTimes, param, spkOkUCueTrials);


                %% save results
                save(saveName,'prefDir', '-append');

                clear cellclassInfo mFiringRate PSTH_f t_r latency_bhv latency_neuro latency_r  latencyStats nLatencyTrials_pref_success
            catch err
                disp(err);
                ng = [ng idata];
                clear cellclassInfo mFiringRate PSTH_f t_r latency_bhv latency_neuro latency_r  mdiffCueFOnset stddiffCueFOnset
                close all;
            end
            save(logName, "ng",'idata');
        end

    end
end
