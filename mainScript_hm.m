%% TOBE DELETED

set(0,'DefaultFigureVisible','off');
[saveServer, rootFolder] = getReady();
%saveServer = 'Z:\Shared\Daisuke\cuesaccade_data';
warning('off','all');

%% recorded data

for aaa = 1%:2
    switch aaa
        case 1
            animal = 'hugo'; %'ollie';%%'andy';% 'andy' '
        case 2
            animal = 'ollie';
    end

    for yyy = 1
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
   

        [loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);

        thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
            regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','03March','23','*_ch29'))))); %2021
        thisdata = [thisdata find(1-cellfun(@isempty, regexp(loadNames, ...
              regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','09September','07','*_ch26')))))]; %2021
        % thisdata = 1:length(channels);


        % parameters
        n=load(fullfile(saveServer,'param20250207.mat'),'param');
         %n=load(fullfile(saveServer,'param20241223.mat'),'param');
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
                load(saveName,'PSTH_f','t_r','spkOkUCueTrials', 'predicted','predicted_all');

                load(loadNames{idata}, 'dd');

                %% prepare behavioral data (common across channels per day)
                eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);

                load(eyeName,'catEvTimes');%,...
      
                y_r = cat(2, PSTH_f, predicted_all, predicted);

                targetTrials_c = find(~isnan(catEvTimes.tOnset) .* ~isnan(catEvTimes.fOnset)); %with or without cue
                winSamps = param.latencyTWin(1):median(diff(t_r)):param.latencyTWin(2);
                periEventTimes = bsxfun(@plus, catEvTimes.tOnset, winSamps); % rows of absolute time points around each event
                okEvents = intersect(find(periEventTimes(:,end)<=max(t_r)), find(periEventTimes(:,1)>=min(t_r)));
                targetTrials = intersect(intersect(targetTrials_c, okEvents), spkOkUCueTrials);


                [fig, avgAmp_hm, p_hm, ranksumval_hm, ranksumz_hm] = showTonsetResp_hm(y_r, ...
                    t_r, catEvTimes, param.figTWin,  param, dd, targetTrials);
                  savePaperFigure(fig, fullfile(saveFigFolder, ['tOnsetResp_hm_' saveSuffix]));
          
                close(fig);
                close all
            
                %% save results
                save(saveName,'ranksumval_hm', 'ranksumz_hm', '-append');

                clear cellclassInfo mFiringRate PSTH_f t_r latency_bhv latency_neuro latency_r  latencyStats nLatencyTrials_pref_success ranksumval_hm ranksumz_hm
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
