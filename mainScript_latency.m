set(0,'DefaultFigureVisible','off');
[saveServer, rootFolder] = getReady();
%saveServer = 'Z:\Shared\Daisuke\cuesaccade_data';


%% recorded data
animal = 'hugo';%'ollie';%'andy';% 'andy' '
tWin_t = [-0.5 0.5];
Thresh = 2;%4; %threshold to define latency

for yyy = 2
    switch yyy
        case 1
            year = '2021';
        case 2
            year = '2022';
        case 3
            year = '2023';
    end

    saveFigFolder = fullfile(saveServer, '20240619',year,animal);
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
     % thisdata = [thisdata find(1-cellfun(@isempty, regexp(loadNames, ...
     %       regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','08August','10','*_ch24')))))]; %2022 vision
     thisdata = 915:1097;

    if isempty(thisdata)
        thisdata = 1:length(channels);
    end


    % parameters
    n=load(fullfile(saveServer,'param20230405_copy.mat'),'param');
    param =n.param;
    n=[];
    ncDirs = length(param.cardinalDir);
    %param.lagRange(2,:)=[-1 0.5];

    psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);

    ng = [];
    previousDate = [];
    for idata = thisdata

        n=load(fullfile(saveServer,'param20230405_copy.mat'),'param');
        param =n.param;
        n=[];

        try
            % datech = [years{idata} filesep months{idata} filesep dates{idata} filesep num2str(channels{idata})];
            datech = [months{idata} filesep dates{idata} filesep num2str(channels{idata})];
            disp([num2str(idata) '/' num2str(numel(thisdata)) ', ' datech ]);

            saveSuffix = [animal replace(datech,filesep,'_') '_linear_rReg'];%'_cue'];

            thisDate = [months{idata} '_' dates{idata}];
        
            saveFolder = fullfile(saveServer, year,animal);%17/6/23
            saveName = fullfile(saveFolder, [saveSuffix '.mat']);

            load(saveName,'mFiringRate');
            if mFiringRate < 5
                clear mFiringRate
                continue;
            end
            load(saveName,'PSTH_f','t_r');

            load(loadNames{idata}, 'dd');

            %% prepare behavioral data (common across channels per day)
            eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
            
            load(eyeName,'eyeData_rmotl_cat','catEvTimes','onsets_cat');%result of processEyeData
            
            [latency_neuro, latency_bhv, latency_r, fig_latency, fig_neurolatency] = ...
                getTgtLatencyCorr(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, Thresh, param, dd);
            screen2png([saveName(1:end-4) '_latency.png'], fig_latency);
            screen2png([saveName(1:end-4) '_neurolatency.png'], fig_neurolatency);
            close(fig_latency);
            close(fig_neurolatency);
            save(saveName,'latency_bhv','latency_neuro','latency_r', '-append');

            
            validEvents = find(~isnan(catEvTimes.tOnset) .* ~isnan(catEvTimes.fOnset)); %with or without cue
            [fig, stats_stratified] = showTonsetResp_stratified(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t,  param, dd, validEvents);
            screen2png([saveName(1:end-4) '_tOnset_stratified.png'], fig);
            close(fig);
           
             [diffCueFOnset, mdiffCueFOnset,stddiffCueFOnset] = getDiffCueTgtOnset(onsets_cat, catEvTimes); %3/6/24
             save(saveName,'stats_stratified','mdiffCueFOnset','stddiffCueFOnset','-append');

            clear PSTH_f t_r latency_bhv latency_neuro latency_r  mdiffCueFOnset stddiffCueFOnset
        catch err
            disp(err);
            ng = [ng idata];
            close all;
        end
    end

end
