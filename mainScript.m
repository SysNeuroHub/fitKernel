set(0,'DefaultFigureVisible','off');
[saveServer, rootFolder] = getReady();
%saveServer = 'Z:\Shared\Daisuke\cuesaccade_data';


%% recorded data
animal =  'hugo';%'ollie';% % %'andy';%
useGPU = 1; %13/12/24
dataType = 0;%0: each channel, 1: all channels per day
for yyy = 2
    switch yyy
        case 1
            year = '2021'; %hugo
        case 2
            year = '2022'; %hugo
        case 3
            year = '2023'; %hugo, ollie
    end

    saveFigFolder = fullfile(saveServer, '20250207',year,animal);
    if ~exist(saveFigFolder, 'dir')
        mkdir(saveFigFolder);
    end

    [loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);

    % to obtain index of specified month&date&channel
    thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
    regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','07July','26','*_ch19')))));
 
    thisdata = thisdata:numel(loadNames);
    nData = numel(thisdata);

    id_pop = cell(nData,1);

    expval_tgt_pop = cell(nData,1);
    corr_tgt_pop = cell(nData,1);
    corr_tgt_rel_pop = cell(nData,1);

    latency_r_pop = cell(nData,1);
    avgAmp_hm_pop = cell(nData,1);
    p_hm_pop = cell(nData,1);
    spkOk_th_pop = cell(nData,1);
    spkOkTrials_pop = cell(nData,1);
    spkOkUCueTrials_pop = cell(nData,1);
    mFiringRate_pop = cell(nData,1);
    PtonsetResp_pop = cell(nData,1);
    ntargetTrials_pop = cell(nData, 1);
    errorIDs = cell(nData,1);

    ng = [];
    previousDate = [];
    for idata = thisdata


        n=load(fullfile(saveServer,'param20250207.mat'),'param');
        param =n.param;
        n=[];
        psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);

        try
            % datech = [years{idata} filesep months{idata} filesep dates{idata} filesep num2str(channels{idata})];
            datech = [months{idata} filesep dates{idata} filesep num2str(channels{idata})];
            thisid = [animal '/' year '/' datech];
            disp(thisid);

            saveSuffix = [animal replace(datech,filesep,'_') '_linear_rReg'];

            thisDate = [months{idata} '_' dates{idata}];
            % if sum(strcmp(thisDate, {'06June_06','06June_11','06June_09'}))>0
            %     %june11  Sample points must be unique.
            %     %june09
            %     %june06 weird blank period in time around 500-600s
            %     continue;
            % end
            saveFolder = fullfile(saveServer, year,animal);%17/6/23
            if ~exist(saveFolder, 'dir')
                mkdir(saveFolder);
            end
            saveName = fullfile(saveFolder, [saveSuffix '.mat']);


            EE = load(loadNames{idata},'ephysdata','dd');
            dd = EE.dd;

            %% prepare predictor variables after downsampling
            %predictorInfoName = fullfile(saveFolder,['predictorInfo_' animal thisDate '.mat']);
            % predictorInfo = preparePredictors(dd, eyeData_rmotl_cat, t_r, param, catEvTimes);
            % save(predictorInfoName, 'predictorInfo');
            % m=matfile(predictorInfoName,'writable',true);
            % m.predictorInfo=predictorInfo;


            spk_all = EE.ephysdata.spikes.spk;
            EE = [];
            if ~isempty(spk_all)
                %% concatenate across trials
                [spk_all_cat, t_cat] = concatenate_spk(spk_all,  dd.eye);
                %clear spk_all
                mFiringRate = length(spk_all_cat)/(t_cat(end)-t_cat(1)); %spks/s
            else
                mFiringRate = 0;
            end
            clear spk_all;
            %clear ephysdata

            if mFiringRate < param.mfr_th
                disp(['skipped as mFiringRate < ' num2str(param.mfr_th)]);
                %save(saveName,'mFiringRate');
                m=matfile(saveName,'writable',true);
                m.FiringRate = mFiringRate;
                clear mFiringRate
                continue;
            end


            %% prepare behavioral data (common across channels per day)
            eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
            if  ~exist(eyeName, 'file') %~strcmp(thisDate, previousDate)

                [eyeData_rmotl_cat, catEvTimes, t_tr, onsets_cat,meta_cat,blinks,outliers] ...
                    = processEyeData(dd.eye, dd, param);

                m=matfile(eyeName,'writable',true);
                m.eyeData_rmotl_cat = eyeData_rmotl_cat;
                m.catEvTimes = catEvTimes;
                m.onsets_cat = onsets_cat;
                m.meta_cat = meta_cat;
                m.blinks = blinks;
                m.outliers=outliers;
                m.t_tr=t_tr;
                close all
            else
                %                 if exist(saveName,'file')
                %                     continue;
                %                 end

                disp('loading eye/predictor data');
                n=load(eyeName,'eyeData_rmotl_cat','catEvTimes',...
                    'onsets_cat','meta_cat','blinks','outliers','t_tr');
                eyeData_rmotl_cat = n.eyeData_rmotl_cat;
                catEvTimes = n.catEvTimes;
                onsets_cat=n.onsets_cat;
                meta_cat=n.meta_cat;
                blinks=n.blinks;
                outliers=n.outliers;
                t_tr=n.t_tr;

                n=[];

                % t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
                % %predictorInfo = preparePredictors(dd, eyeData_rmotl_cat, t_r, param, catEvTimes);
                % n=load(predictorInfoName, 'predictorInfo');
                % predictorInfo = n.predictorInfo;
                % n=[];
                %load(fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']));
            end



            %% prepare predictor variables after downsampling
            t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
            predictorInfoName = fullfile(saveFolder,['predictorInfo_' animal thisDate '.mat']);
            % if exist(predictorInfoName, 'file')
            %     n=load(predictorInfoName, 'predictorInfo');
            %     predictorInfo = n.predictorInfo;
            %     n=[];
            % else
                predictorInfo = preparePredictors(dd, eyeData_rmotl_cat, t_r, param, catEvTimes);
                save(predictorInfoName, 'predictorInfo');
            % end

            %% nomalize predictor variables so the resulting kernels are comparable 28/1/25
            predictorInfo.predictors_r(9:end,:) = normalize(predictorInfo.predictors_r(9:end,:),2,'zscore');
          
            %% remove trials with too short duration 28/10/23
            %[t_tr, catEvTimes, validTrials] = trimInvalids(t_tr,
            %catEvTimes); %commented out 17/12/24

            [trIdx_r] = retrieveTrIdx_r(t_cat, t_r, t_tr);

            %% detect trials where firing rate is extremely large
            [spkOkTrials, spkOk_th] = getSpkOKtrials(spk_all_cat, t_r, trIdx_r, param);

            [spkOkUCueTrace, spkOkUCueTrials] = getIncludeTrace(t_cat, t_r, t_tr, onsets_cat, spkOkTrials);
            spkNGRate = (numel(t_tr)-numel(spkOkTrials))/numel(t_tr)*100;
            CueTrRate = (numel(spkOkTrials)-numel(spkOkUCueTrials))/numel(spkOkTrials)*100;

            ntargetTrials = numel(intersect(find(~isnan(catEvTimes.tOnset)), spkOkUCueTrials));


            %%  detect saccades outside the task
            tOnset = catEvTimes.tOnset;
            cOnset = catEvTimes.cOnset; %choice onset not cue %% WHY THIS CONDITION??
            validEvents = intersect(find(~isnan(tOnset)), find(~isnan(cOnset)));
            tOnset = tOnset(validEvents);
            cOnset = cOnset(validEvents);

            tcOnset_trace = event2Trace(t_cat, [tOnset; cOnset], 2*0.5);
            excEventT_cat = (tcOnset_trace + blinks + outliers > 0); %28/1/22
            [startSaccNoTask, endSaccNoTask] = selectSaccades(catEvTimes.saccadeStartTimes, ...
                catEvTimes.saccadeEndTimes, t_cat, excEventT_cat);%param.minSaccInterval);
            saccNoTaskTrace = event2Trace(t_cat, [startSaccNoTask endSaccNoTask]);
            spkOkUCueTrace_tmp = getIncludeTrace(t_cat, t_cat, t_tr, onsets_cat, spkOkTrials);
            saccNoTask_spkOkUCue_Trace = saccNoTaskTrace.*spkOkUCueTrace_tmp;

            [~, startSaccNoTask_spkOkUCue, endSaccNoTask_spkOkUCue] = trace2Event(saccNoTask_spkOkUCue_Trace, t_cat);
            [saccDirNoTask_spkOkUCue, dirIndexNoTask_spkOkUCue] = getSaccDir(startSaccNoTask_spkOkUCue, endSaccNoTask_spkOkUCue, ...
                eyeData_rmotl_cat, param.cardinalDir);

            %% obtain kernels!
            % if exist(saveName,'file')
            %      disp('load kernel fit results'); `
            %      load (saveName, 'predicted_all', 'predicted', 'PSTH_f', 'kernelInfo','mFiringRate');
            %  else
            disp('fit kernels');
            [predicted_all, predicted, PSTH_f, kernelInfo] = fitPSTH_cv(spk_all_cat, ...
                predictorInfo.t_r, param.predictorNames,   predictorInfo.predictors_r, ...
                predictorInfo.npredVars,param.psth_sigma, param.kernelInterval, ...
                param.lagRange, param.ridgeParams, trIdx_r(spkOkUCueTrials),param.fitoption,useGPU, ...
                param.kfolds);
            % end
            y_r = cat(2,PSTH_f,predicted_all, predicted);


            %% explained variance for target response
            nPredictorNames = numel(param.predictorNames);

            expval_tgt = zeros(nPredictorNames, 1);
            corr_tgt = zeros(nPredictorNames, 1);
            [expval_tgt(1,1), corr_tgt(1,1)] = ...
                getExpVal_tgt(PSTH_f, predicted_all, catEvTimes, t_r, param.tOnRespWin, spkOkUCueTrials);
            [expval_tgt(2:nPredictorNames+1,1), corr_tgt(2:nPredictorNames+1,1)] = ...
                getExpVal_tgt(PSTH_f, predicted, catEvTimes, t_r, param.tOnRespWin, spkOkUCueTrials);
            corr_tgt_rel = 100*corr_tgt(2:3)./corr_tgt(1);


            %% figure for kernel fitting
            f = showKernel( t_r, y_r, kernelInfo, param.cardinalDir);
            screen2png(fullfile(saveFigFolder,['kernels_exp' saveSuffix]), f);
            close(f);

            %% Figure for target onset response (only to preferred direction)
            [f, cellclassInfo] = showTonsetResp(t_r, y_r, catEvTimes, dd, psthNames, ...
                startSaccNoTask_spkOkUCue, saccDirNoTask_spkOkUCue, param, [], spkOkUCueTrials);
            cellclassInfo.datech = datech;
            %savePaperFigure(f, fullfile(saveFigFolder,['cellclassFig_' saveSuffix]));
            screen2png(fullfile(saveFigFolder,['cellclassFig_' saveSuffix '_allTr']), f);
            close(f);


            %% single-trial latency
            %temporal window was [-0.5 0.5] in early 2024
            % only use events whose time window is within the recording (from eventLockedAvg)
            targetTrials_c = find(~isnan(catEvTimes.tOnset) .* ~isnan(catEvTimes.fOnset)); %with or without cue
            winSamps = param.latencyTWin(1):median(diff(t_r)):param.latencyTWin(2);
            periEventTimes = bsxfun(@plus, catEvTimes.tOnset, winSamps); % rows of absolute time points around each event
            okEvents = intersect(find(periEventTimes(:,end)<=max(t_r)), find(periEventTimes(:,1)>=min(t_r)));
            targetTrials = intersect(intersect(targetTrials_c, okEvents), spkOkUCueTrials);

            [latency_neuro, thresh_neuro, tgtDir, fig_latency, nLatencyTrials_pref_success, latencyStats] = ...
                getTgtNeuroLatency(PSTH_f, t_r, onsets_cat, catEvTimes, param.latencyTWin, ...
                param.threshParam, param, dd, targetTrials);
            screen2png(fullfile(saveFigFolder, ['latencyCorr_' saveSuffix]), fig_latency);close(fig_latency);


            %% target response hit v miss
            [fig, avgAmp_hm, p_hm] = showTonsetResp_hm(y_r, t_r, onsets_cat, catEvTimes, param.figTWin,  ...
                param, dd, targetTrials);

            screen2png(fullfile(saveFigFolder, ['tOnsetResp_hm_' saveSuffix]), fig);
            close(fig);

            %% for population analysis
            id_pop{idata} = thisid;

            expval_tgt_pop{idata} = expval_tgt;
            corr_tgt_pop{idata} = corr_tgt;
            corr_tgt_rel_pop{idata} = corr_tgt_rel;
            latency_r_pop{idata} = latencyStats.latency_r;
            avgAmp_hm_pop{idata} = avgAmp_hm;
            p_hm_pop{idata} = p_hm;
            spkOk_th_pop{idata} = spkOk_th;
            spkOkTrials_pop{idata} = spkOkTrials;
            spkOkUCueTrials_pop{idata} = spkOkUCueTrials;
            mFiringRate_pop{idata} = mFiringRate;
            PtonsetResp_pop{idata} = cellclassInfo.PtonsetResp;
            ntargetTrials_pop{idata} = ntargetTrials;
            errorIDs{idata} = 0;
            spkNGRate_pop{idata} = spkNGRate;
            CueTrRate_pop{idata} = CueTrRate;
            nLatencyTrials_pref_success_pop{idata} = nLatencyTrials_pref_success;


            %% save results
            mm=matfile(saveName,'writable',true);
            mm.PSTH_f = PSTH_f;
            mm.predicted_all = predicted_all;
            mm.predicted = predicted;
            mm.kernelInfo = kernelInfo;
            mm.t_r = t_r;
            mm.param=param;
            mm.mFiringRate=mFiringRate;
            mm.t_cat=t_cat;
            mm.dd=[];%dd; 31/1/25 dd is too large can be >7GB
            mm.latencyStats = latencyStats;
            mm.avgAmp_hm = avgAmp_hm;
            mm.p_hm = p_hm;
            mm.expval_tgt = expval_tgt;
            mm.corr_tgt = corr_tgt;
            mm.corr_tgt_rel = corr_tgt_rel;
            mm.avgAmp_hm = avgAmp_hm;
            mm.p_hm = p_hm;
            mm.spkOk_th = spkOk_th;
            mm.spkOkTrials = spkOkTrials;
            mm.spkOkUCueTrials = spkOkUCueTrials;
            mm.cellclassInfo = cellclassInfo;
            mm.ntargetTrials = ntargetTrials;
            mm.spkNGRate = spkNGRate;
            mm.CueTrRate = CueTrRate;
            mm.nLatencyTrials_pref_success = nLatencyTrials_pref_success;
            clear mm mFiringRate;
            close all
        catch err
            clear mFiringRate

            disp(err);
            errorIDs{idata} = 1;
            ng = [ng idata];
            close all;
        end
    end

    % % save(fullfile(saveFolder, 'assembly20241212.mat'),'param',...
    % %     'id_pop','expval_tgt_pop','corr_tgt_pop','corr_tgt_rel_pop',...
    % %     'latency_r_pop','avgAmp_hm_pop','p_hm_pop','spkOk_th_pop',...
    % %     'spkOkTrials_pop','spkOkUCueTrials_pop','mFiringRate_pop',...
    % %     'PtonsetResp_pop','errorIDs');
    % % assembly = [];
end
