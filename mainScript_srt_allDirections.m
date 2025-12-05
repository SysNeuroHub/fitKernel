%% created from mainScript_latency.m
% TOBE DELETED

set(0,'DefaultFigureVisible','off');
[saveServer, rootFolder] = getReady();
%saveServer = 'Z:\Shared\Daisuke\cuesaccade_data';
warning('off','all');

for aa = 1:2
    switch aa
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

        saveFigFolder = fullfile(saveServer, '20250207',year,animal);
        mkdir(saveFigFolder);


        [loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);

        thisdata = 1:length(channels);


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
                load(saveName,'t_r','spkOkUCueTrials');

                load(loadNames{idata}, 'dd');

                %% prepare behavioral data (common across channels per day)
                eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);

                load(eyeName,'catEvTimes','onsets_cat');%,...
                % 'eyeData_rmotl_cat','startSaccNoTask', 'saccDirNoTask');%result of processEyeData

                targetTrials_c = find(~isnan(catEvTimes.tOnset) .* ~isnan(catEvTimes.fOnset)); %with or without cue
                winSamps = param.latencyTWin(1):median(diff(t_r)):param.latencyTWin(2);
                periEventTimes = bsxfun(@plus, catEvTimes.tOnset, winSamps); % rows of absolute time points around each event
                okEvents = intersect(find(periEventTimes(:,end)<=max(t_r)), find(periEventTimes(:,1)>=min(t_r)));
                targetTrials = intersect(intersect(targetTrials_c, okEvents), spkOkUCueTrials);

                % [latency_neuro, thresh_neuro, tgtDir, fig_latency, nLatencyTrials_pref_success, latencyStats] = ...
                %     getTgtNeuroLatency(PSTH_f, t_r, onsets_cat, catEvTimes, param.latencyTWin, ...
                %     param.threshParam, param, dd, targetTrials);
                % savePaperFigure(fig_latency, fullfile(saveFigFolder, ['latencyCorr_' saveSuffix]));close(fig_latency);
                %delete(fullfile(saveFigFolder, ['latencySingle_' saveSuffix]));


                srt_allDirections.latency_bhv_srt = getTgtBhvLatency(onsets_cat, dd, targetTrials, 1);
                srt_allDirections.tgtDir = getTgtDir(dd.targetloc(targetTrials), param.cardinalDir);


                %% save results
                save(saveName,'srt_allDirections', '-append');

                clear srt_allDirections
            catch err
                disp(err);
                ng = [ng idata];
                clear mFiringRate PSTH_f t_r srt_allDirections
                close all;
            end
        end

    end
end
