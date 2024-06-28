set(0,'DefaultFigureVisible','off');
[saveServer, rootFolder] = getReady();
%saveServer = 'Z:\Shared\Daisuke\cuesaccade_data';
warning('off','all');

%% recorded data
animal = 'hugo';%'ollie';%'andy';% 'andy' '
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
    logName = fullfile(saveServer,'20240625',year,animal,'log_mainScript_latency');

    saveFigFolder = fullfile(saveServer, '20240625',year,animal);
    mkdir(saveFigFolder);


    [loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);

    % to obtain index of specified month&date&channel
    % thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
    %        regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','11November','25','*_ch24'))))); %2021 no corr
    % thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
    %        regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','09September','16','*_ch9'))))); %2022 positive corr
        % thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
        %    regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','09September','19','*_ch20'))))); %2022 positive corr
     % thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
     %       regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','09September','19','*_ch4'))))); %2022 no corr
     % thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
     %       regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','09September','19','*_ch15'))))); %2022 no corr
     % thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
     %       regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','02February','25','*_ch24'))))); %2022 vision
     % thisdata = [thisdata find(1-cellfun(@isempty, regexp(loadNames, ...
     %       regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','08August','09','*_ch3')))))]; %2022 vision
     % thisdata =  find(1-cellfun(@isempty, regexp(loadNames, ...
     %       regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','08August','25','*_ch27'))))); %2022 vision
     % thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
     %       regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','07July','08','*_ch17'))))); %2022
     % thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
     %     regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','03March','22','*_ch21')))));
     % thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
     %     regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','08August','24','*_ch27'))))); %2021 non bhv
        thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
         regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','12December','13','*_ch*')))));

     %thisdata = 1:length(channels);
  

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
            load(saveName,'PSTH_f','t_r');

            load(loadNames{idata}, 'dd');

            %% prepare behavioral data (common across channels per day)
            eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
            
            load(eyeName,'eyeData_rmotl_cat','catEvTimes','onsets_cat',...
                'startSaccNoTask', 'saccDirNoTask');%result of processEyeData
            
            %% only use events whose time window is within the recording (from eventLockedAvg)
            validEvents_c = find(~isnan(catEvTimes.tOnset) .* ~isnan(catEvTimes.fOnset)); %with or without cue
            winSamps = tWin_t(1):median(diff(t_r)):tWin_t(2);
            periEventTimes = bsxfun(@plus, catEvTimes.tOnset, winSamps); % rows of absolute time points around each event
            okEvents = intersect(find(periEventTimes(:,end)<=max(t_r)), find(periEventTimes(:,1)>=min(t_r)));
            validEvents = intersect(validEvents_c, okEvents);

            %% single-trial latency
            [latency_neuro, latency_bhv, latency_r, fig_latency, fig_neurolatency] = ...
                getTgtLatencyCorr(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, param, ...
                dd, validEvents);
            screen2png(fullfile(saveFigFolder, ['latencyCorr_' saveSuffix]), fig_latency);
            screen2png(fullfile(saveFigFolder, ['latencySingle_' saveSuffix]), fig_neurolatency);
            close(fig_latency);
            close(fig_neurolatency);
            save(saveName,'latency_bhv','latency_neuro','latency_r', '-append');

            %% avg-trial latency after stratification
            [fig, stats_stratified] = showTonsetResp_stratified(PSTH_f, t_r, onsets_cat, catEvTimes, ...
                tWin_t,  param, dd, validEvents, param.nDiv);
            screen2png(fullfile(saveFigFolder, ['tOnset_stratified_' saveSuffix]), fig)
            close(fig);
           
            %% ?
             [diffCueFOnset, mdiffCueFOnset,stddiffCueFOnset] = ...
                 getDiffCueTgtOnset(onsets_cat, catEvTimes); %3/6/24


             %% tgt response to obtain cellClass info
             load(saveName,'predicted_all','predicted');
             y_r = cat(2,PSTH_f, predicted_all, predicted);
             [f, cellclassInfo] = showTonsetResp(t_r, y_r, catEvTimes, dd, psthNames, ...
                 startSaccNoTask, saccDirNoTask, param);
             cellclassInfo.datech = datech;
                screen2png(fullfile(saveFigFolder,['cellclassFig_' saveSuffix '_allTr']), f);
                close(f);

             %% save results
             save(saveName,'stats_stratified','mdiffCueFOnset','stddiffCueFOnset','cellclassInfo','-append');

              clear cellclassInfo mFiringRate PSTH_f t_r latency_bhv latency_neuro latency_r  mdiffCueFOnset stddiffCueFOnset
        catch err
            disp(err);
            ng = [ng idata];
            clear cellclassInfo mFiringRate PSTH_f t_r latency_bhv latency_neuro latency_r  mdiffCueFOnset stddiffCueFOnset
            close all;
        end
        save(logName, "ng",'idata');
    end

end
