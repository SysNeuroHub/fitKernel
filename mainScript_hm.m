set(0,'DefaultFigureVisible','off');
[saveServer, rootFolder] = getReady();
%saveServer = 'Z:\Shared\Daisuke\cuesaccade_data';
warning('off','all');

%% recorded data
animal = 'hugo';%'ollie';%'andy';% 'andy' '
tWin_t = [-0.8 0.5]; %17/7/24


for yyy = 1:3
    switch yyy
        case 1
            year = '2021';
        case 2
            year = '2022';
        case 3
            year = '2023';
    end
    logName = fullfile(saveServer,'20240715',year,animal,'log_mainScript_latency');

    saveFigFolder = fullfile(saveServer, '20240715',year,animal);
    mkdir(saveFigFolder);


    [loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);

    % to obtain index of specified month&date&channel
    % thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
    %     regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','09September','16','*_ch9'))))); %2022 positive corr
    % thisdata = [thisdata find(1-cellfun(@isempty, regexp(loadNames, ...
    %     regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','09September','19','*_ch20')))))]; %2022 positive corr
    % thisdata = [thisdata find(1-cellfun(@isempty, regexp(loadNames, ...
    %     regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','09September','19','*_ch4')))))]; %2022 no corr
    % thisdata = [thisdata find(1-cellfun(@isempty, regexp(loadNames, ...
    %     regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','09September','19','*_ch15')))))]; %2022 no corr
    % thisdata = [thisdata find(1-cellfun(@isempty, regexp(loadNames, ...
    %     regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','02February','25','*_ch24')))))]; %2022 vision
    % thisdata = [find(1-cellfun(@isempty, regexp(loadNames, ...
    %    regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','03March','16','*_ch29')))))]; %2022 vision

     thisdata = 1:length(channels);

    % parameters
    %n=load(fullfile(saveServer,'param20230405_copy.mat'),'param');
    n=load(fullfile(saveServer,'param20240625.mat'),'param');
    param =n.param;
    n=[];
    ncDirs = length(param.cardinalDir);
    %param.lagRange(2,:)=[-1 0.5];

    psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);

    ng = [];
    previousDate = [];
    for idata = thisdata

        % n=load(fullfile(saveServer,'param20230405_copy.mat'),'param');
        % param =n.param;
        % n=[];

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
            load(saveName,'PSTH_f','t_r', 'predicted_all');
            y_r = cat(2, PSTH_f, predicted_all);

            load(loadNames{idata}, 'dd');

            %% prepare behavioral data (common across channels per day)
            eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
            
            load(eyeName,'eyeData_rmotl_cat','catEvTimes','onsets_cat');%result of processEyeData
            
            %% only use events whose time window is within the recording (from eventLockedAvg)
            validEvents_c = find(~isnan(catEvTimes.tOnset) .* ~isnan(catEvTimes.fOnset)); %with or without cue
            winSamps = tWin_t(1):median(diff(t_r)):tWin_t(2);
            periEventTimes = bsxfun(@plus, catEvTimes.tOnset, winSamps); % rows of absolute time points around each event
            okEvents = intersect(find(periEventTimes(:,end)<=max(t_r)), find(periEventTimes(:,1)>=min(t_r)));
            validEvents = intersect(validEvents_c, okEvents);

         
            %% target response hit v miss 
             [fig, avgAmp_hm, p_hm] = showTonsetResp_hm(y_r, t_r, onsets_cat, catEvTimes, tWin_t,  ...
                param, dd, validEvents);
            screen2png(fullfile(saveFigFolder, ['tOnsetResp_hm_' saveSuffix]), fig);
            close(fig);
            save(saveName,'avgAmp_hm', 'p_hm','-append');

             close all
              clear mFiringRate y_r avgAmp_hm p_hm
        catch err
            disp(err);
            ng = [ng idata];
            clear y_r mFiringRate PSTH_f t_r clear avgAmp_hm p_hm

            close all;
        end
        save(logName, "ng",'idata');
    end

end
