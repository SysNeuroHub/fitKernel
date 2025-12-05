rootFolder =  '/media/daisuke/cuesaccade_data';

%% recorded data
animal =  'hugo';%'ollie';%
useGPU = 1; %13/12/24
dataType = 0;%0: each channel, 1: all channels per day


saveFigFolder = fullfile(rootFolder, animal);
mkdir(saveFigFolder);

for idata = 1:3

    n=load(fullfile(rootFolder,'param20250207.mat'),'param');
    param =n.param;
    n=[];
    psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);

    switch idata
        case 1
            year = '2021';
            month = '08August';
            date = '25';
            channel = 27;
        case 2
            year = '2022';
            month = '07July';
            date = '26';
            channel = 19;
        case 3
            year = '2022';
            month = '08August';
            date = '15';
            channel = 4;
    end

    loadName = fullfile(rootFolder,year,'cuesaccade_data',month,date,'saved_oephysdata',[animal '_oephysdata_ch' num2str(channel) '.mat']);
    
    datech = [month filesep date filesep num2str(channel)];
    thisid = [animal '/' year '/' datech];
    disp(thisid);

    saveSuffix = [animal replace(datech,filesep,'_') '_linear_rReg'];

    thisDate = [month '_' date];
    saveFolder = fullfile(rootFolder, year,animal);%17/6/23
    if ~exist(saveFolder, 'dir')
        mkdir(saveFolder);
    end

    EE = load(loadName,'ephysdata','dd');
    dd = EE.dd;

    %% prepare predictor variables after downsampling
    spk_all = EE.ephysdata.spikes.spk;
    EE = [];
    if ~isempty(spk_all)
        %% concatenate across trials
        [spk_all_cat, t_cat] = concatenate_spk(spk_all,  dd.eye);
        mFiringRate = length(spk_all_cat)/(t_cat(end)-t_cat(1)); %spks/s
    else
        mFiringRate = 0;
    end
    clear spk_all;

    %% prepare behavioral data (common across channels per day)
    eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);

    [eyeData_rmotl_cat, catEvTimes, t_tr, onsets_cat,meta_cat,blinks,outliers] ...
        = processEyeData(dd.eye, dd, param);

    %% prepare predictor variables after downsampling
    t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
    predictorInfoName = fullfile(saveFolder,['predictorInfo_' animal thisDate '.mat']);

    predictorInfo = preparePredictors(dd, eyeData_rmotl_cat, t_r, param, catEvTimes);
    save(predictorInfoName, 'predictorInfo');

    %% nomalize predictor variables so the resulting kernels are comparable 28/1/25
    predictorInfo.predictors_r(9:end,:) = normalize(predictorInfo.predictors_r(9:end,:),2,'zscore');

    [trIdx_r] = retrieveTrIdx_r(t_cat, t_r, t_tr);

    %% detect trials where firing rate is extremely large
    [spkOkTrials, spkOk_th] = getSpkOKtrials(spk_all_cat, t_r, trIdx_r, param);

    [spkOkUCueTrace, spkOkUCueTrials] = getIncludeTrace(t_cat, t_r, t_tr, onsets_cat, spkOkTrials);
    spkNGRate = (numel(t_tr)-numel(spkOkTrials))/numel(t_tr)*100;
    CueTrRate = (numel(spkOkTrials)-numel(spkOkUCueTrials))/numel(spkOkTrials)*100;

    ntargetTrials = numel(intersect(find(~isnan(catEvTimes.tOnset)), spkOkUCueTrials));


    %%  detect saccades outside the task
    [startSaccNoTask_spkOkUCue, endSaccNoTask_spkOkUCue, saccDirNoTask_spkOkUCue, dirIndexNoTask_spkOkUCue] = ...
        getSaccNoTask(t_cat, catEvTimes, eyeData_rmotl_cat, blinks, outliers, t_tr, onsets_cat, spkOkTrials, param);


    %% obtain kernels!
    disp('fit kernels');
    [predicted_all, predicted, PSTH_f, kernelInfo] = fitPSTH_cv(spk_all_cat, ...
        predictorInfo.t_r, param.predictorNames,   predictorInfo.predictors_r, ...
        predictorInfo.npredVars,param.psth_sigma, param.kernelInterval, ...
        param.lagRange, param.ridgeParams, trIdx_r(spkOkUCueTrials),param.fitoption,useGPU, ...
        param.kfolds);

    y_r = cat(2,PSTH_f,predicted_all, predicted);


    %% preferred direction
    param_tmp = param;
    param_tmp.cardinalDir = 0:359;
    prefDir = getPrefDir_wrapper(PSTH_f, t_r, dd, catEvTimes, param_tmp, spkOkUCueTrials);


    %% explained variance for target response
    nPredictorNames = numel(param.predictorNames);

    expval_tgt = zeros(nPredictorNames, 1);
    corr_tgt = zeros(nPredictorNames, 1);
    [expval_tgt(1,1), corr_tgt(1,1)] = ...
        getExpVal_tgt(PSTH_f, predicted_all, catEvTimes, t_r, param.tOnRespWin, spkOkUCueTrials);
    [expval_tgt(2:nPredictorNames+1,1), corr_tgt(2:nPredictorNames+1,1)] = ...
        getExpVal_tgt(PSTH_f, predicted, catEvTimes, t_r, param.tOnRespWin, spkOkUCueTrials);
    corr_tgt_rel = 100*corr_tgt(2:4)./corr_tgt(1);


    %% Figure for target onset response (only to preferred direction)
    [f, cellclassInfo] = showTonsetResp(t_r, y_r, catEvTimes, dd, psthNames, ...
        startSaccNoTask_spkOkUCue, saccDirNoTask_spkOkUCue, param, [], spkOkUCueTrials);
    cellclassInfo.datech = datech;
    screen2png(fullfile(saveFigFolder,['cellclassFig_' saveSuffix]), f);
    close(f);

end
