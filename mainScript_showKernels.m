set(0,'DefaultFigureVisible','off');
[saveServer, rootFolder] = getReady();
%saveServer = 'Z:\Shared\Daisuke\cuesaccade_data';
warning('off','all');

%% recorded data
animal = 'hugo'; %'ollie';%%'andy';% 'andy' '
tWin_t = [-0.5 0.5];


for yyy = 1
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
           regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','12December','16','*_ch26'))))); %2021
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
            load(saveName, 't_r','t_cat','spkOkUCueTrials','kernelInfo','PSTH_f','predicted_all','predicted');

            y_r = cat(2,PSTH_f,predicted_all, predicted);

            eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
            load(eyeName,'t_tr');

            predictorInfoName = fullfile(saveFolder,['predictorInfo_' animal thisDate '.mat']);
            load(predictorInfoName, 'predictorInfo');

            [trIdx_r] = retrieveTrIdx_r(t_cat, t_r, t_tr);

            t_r_idx = [];
            for itr = spkOkUCueTrials
                t_r_idx = [t_r_idx; trIdx_r{itr}];
            end

            predictorInfo.t_r = predictorInfo.t_r(t_r_idx);
            predictorInfo.predictors_r = predictorInfo.predictors_r(:,t_r_idx);

            kernelInfo_norm = getKernelInfo_norm(kernelInfo, predictorInfo);

            f = showKernel( t_r, y_r, kernelInfo_norm, param.cardinalDir);

            screen2png(fullfile(saveFigFolder,['normkernels_exp' saveSuffix]), f);
            close(f);


            %% save results
            save(saveName,'kernelInfo_norm', '-append');

            clear kernelInfo kernelInfo_norm mFiringRate t_r t_tr predictorInfo
        catch err
            disp(err);
            ng = [ng idata];
            clear kernelInfo kernelInfo_norm mFiringRate t_r t_tr predictorInfo
            close all;
        end
        save(logName, "ng",'idata');
    end

end
