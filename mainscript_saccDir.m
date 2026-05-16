
[saveServer, rootFolder] = getReady();
saveSuffix_p = ['20250207'];
n=load(fullfile(saveServer, ['param' saveSuffix_p '.mat']),'param');
param =n.param;

for aa = 1:2
    switch aa
        case 1
            animal = 'hugo'; yidx=1:3;
        case 2
            animal =  'ollie'; yidx = 3;
    end
    for yy = yidx
        switch yy
            case 1
                year = '2021';
            case 2
                year = '2022';
            case 3
                year = '2023';
        end
        [loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);
        nData = numel(loadNames);
        thisdata = 1:nData;


        saveFolder = fullfile(saveServer, year,animal);%17/6/23
     

        for idata = thisdata

            datech = [months{idata} filesep dates{idata} filesep num2str(channels{idata})];
            thisid = [animal '/' year '/' datech];
            disp(thisid);

            saveSuffix = [animal replace(datech,filesep,'_') '_linear_rReg'];

            saveName = fullfile(saveFolder, [saveSuffix '.mat']);
             if ~exist(saveName,'file'); continue; end;
          
            load(saveName, 't_r', 'spkOkTrials','t_cat');
            if ~exist('t_r','var'); continue; end;

            %% prepare behavioral data (common across channels per day)
            thisDate = [months{idata} '_' dates{idata}];
            
            eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
            load(eyeName, 'eyeData_rmotl_cat','catEvTimes','blinks','outliers',...
                'onsets_cat','t_tr');


            %%  detect saccades outside the task
            [startSaccNoTask_spkOkUCue, endSaccNoTask_spkOkUCue, saccDirNoTask_spkOkUCue, dirIndexNoTask_spkOkUCue] = ...
                getSaccNoTask(t_cat, catEvTimes, eyeData_rmotl_cat, blinks, outliers, ...
                t_tr, onsets_cat, spkOkTrials, param);

            close all

            save(saveName,'saccDirNoTask_spkOkUCue','-append');

            % %% Figure for saccade direction outside task
            binEdges = [param.cardinalDir 360] - 0.5*mean(diff(param.cardinalDir)); %cf.getSaccDir
            fhist=figure('position',[0 0 200 200]);
            histogram(saccDirNoTask_spkOkUCue, 'BinEdges', binEdges);
            % set(gca,'tickdir','out','XTick',param.cardinalDir);
            % ylabel('# saccades');xlabel('Direction [deg]');
            % savePaperFigure(fhist, fullfile(saveFigFolder,['sptSaccDir_' saveSuffix '_pupil_blink']));
            % close(fhist);

            clear 'eyeData_rmotl_cat' 'catEvTimes' 'blinks' 'outliers' 'onsets_cat' 't_tr' 't_r' 'spkOkTrials' 't_cat'
        end
    end
end