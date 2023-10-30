set(0,'DefaultFigureVisible','off');
[saveServer, rootFolder] = getReady();
%saveServer = 'Z:\Shared\Daisuke\cuesaccade_data';


%% recorded data
animal = 'hugo';%'ollie';%'andy';% 'andy' '
% fitoption = 1; %'linear'
fitoption = 5; %linear_rReg', as of 13/7/2023
useGPU = 1;
dataType = 0;%0: each channel, 1: all channels per day
fitIt = 1;
splitPredictor = 1; %whether to split predictors by cue. 28/10/2023


for yyy = 2%1:3
    switch yyy
        case 1
            year = '2021'; 
        case 2 
            year = '2022';
        case 3
            year = '2023';
    end 
    
    saveFigFolder = fullfile(saveServer, '20230713',year,animal);
    mkdir(saveFigFolder);
    
    
    [loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);
    
    % to obtain index of specified month&date&channel
    % thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
    %     regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','02February','08','*_ch18')))));
    thisdata = [];
    if isempty(thisdata)
        thisdata = 1:length(channels);
    end
    
    %% omit data
    % no saccade response
    % low spontaneous firing
    % low number of successful trials
    
    % parameters
    n=load(fullfile(saveServer,'param20230405.mat'),'param');
    param =n.param;
    n=[];
    ncDirs = length(param.cardinalDir);
    %param.lagRange(2,:)=[-1 0.5];
    
    psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);
    
    ng = [];
    previousDate = [];
    for idata = thisdata(49:end)

        n=load(fullfile(saveServer,'param20230405.mat'),'param');
        param =n.param;
        n=[];

        try
            % datech = [years{idata} filesep months{idata} filesep dates{idata} filesep num2str(channels{idata})];
            datech = [months{idata} filesep dates{idata} filesep num2str(channels{idata})];
            disp([num2str(idata) '/' num2str(numel(thisdata)) ', ' datech ]);
            
            saveSuffix = [animal replace(datech,filesep,'_') '_linear_rReg'];%'_cue'];
            
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
            delete(saveName); %TEMP
            predictorInfoName = fullfile(saveFolder,['predictorInfo_' animal thisDate '.mat']);
           
            EE = load(loadNames{idata},'ephysdata','dd');
            dd = EE.dd;
            
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
            
            if mFiringRate < 5
                disp(['skipped as mFiringRate<5']);
                %save(saveName,'mFiringRate');
                m=matfile(saveName,'writable',true);
                m.FiringRate = mFiringRate;
                continue;
            end
            
            
            %% prepare behavioral data (common across channels per day)
            eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
            if  ~exist(eyeName, 'file') %~strcmp(thisDate, previousDate)
                
                [eyeData_rmotl_cat, catEvTimes, t_tr, onsets_cat,meta_cat,blinks,outliers] ...
                    = processEyeData(dd.eye, dd, param);
                %             [pspec_parea,faxis_parea] =
                %             pmtm(eyeData_rmotl_cat.parea, 10, ...
                %                 length(eyeData_rmotl_cat.parea), fs_eye);%slow
                
                tOnset = catEvTimes.tOnset;
                cOnset = catEvTimes.cOnset; %choice onset not cue
                validEvents = intersect(find(~isnan(tOnset)), find(~isnan(cOnset)));
                tOnset = tOnset(validEvents);
                cOnset = cOnset(validEvents);
                
                tcOnset_trace = event2Trace(t_cat, [tOnset; cOnset], 2*0.5);
                excEventT_cat = (tcOnset_trace + blinks + outliers > 0); %28/1/22
                
                [startSaccNoTask, endSaccNoTask] = selectSaccades(catEvTimes.saccadeStartTimes, ...
                    catEvTimes.saccadeEndTimes, t_cat, excEventT_cat);%param.minSaccInterval);
                %<slow
                
                [saccDirNoTask, dirIndexNoTask] = getSaccDir(startSaccNoTask, endSaccNoTask, ...
                    eyeData_rmotl_cat, param.cardinalDir);
                %<slow
                
                %             save(fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']), 'startSaccNoTask', 'endSaccNoTask', ...
                %                 'saccDirNoTask', 'dirIndexNoTask','-append');
                

                % save(eyeName,'eyeData_rmotl_cat','catEvTimes',...
                %     'onsets_cat','meta_cat','blinks','outliers','t_tr',...
                %     'startSaccNoTask', 'endSaccNoTask', ...
                %     'saccDirNoTask', 'dirIndexNoTask');
                m=matfile(eyeName,'writable',true);
                m.eyeData_rmotl_cat = eyeData_rmotl_cat;
                m.catEvTimes = catEvTimes;
                m.onsets_cat = onsets_cat;
                m.meta_cat = meta_cat;
                m.blinks = blinks;
                m.outliers=outliers;
                m.t_tr=t_tr;
                m.startSaccNoTask=startSaccNoTask;
                m.endSaccNoTask=endSaccNoTask;
                m.saccDirNoTask=saccDirNoTask;
                m.dirIndexNoTask=dirIndexNoTask;
                close all
                
                %% prepare predictor variables after downsampling
                t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
                predictorInfo = preparePredictors(dd, eyeData_rmotl_cat, t_r, param, catEvTimes);
                %save(predictorInfoName, 'predictorInfo');
                m=matfile(predictorInfoName,'writable',true);
                m.predictorInfo=predictorInfo;
            else
%                 if exist(saveName,'file')
%                     continue;
%                 end
                
                disp('loading eye/predictor data');
                %load(fullfile(saveFolder,['predictorInfo_' animal thisDate '.mat']), 'predictorInfo');
                n=load(eyeName,'eyeData_rmotl_cat','catEvTimes',...
                    'onsets_cat','meta_cat','blinks','outliers','t_tr',...
                    'startSaccNoTask', 'endSaccNoTask', ...
                    'saccDirNoTask', 'dirIndexNoTask');
                eyeData_rmotl_cat = n.eyeData_rmotl_cat;
                catEvTimes = n.catEvTimes;
                onsets_cat=n.onsets_cat;
                meta_cat=n.meta_cat;
                blinks=n.blinks;
                outliers=n.outliers;
                t_tr=n.t_tr;
                startSaccNoTask=n.startSaccNoTask;
                endSaccNoTask=n.endSaccNoTask;
                saccDirNoTask=n.saccDirNoTask;
                dirIndexNoTask=n.dirIndexNoTask;
                n=[];

                t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
                %predictorInfo = preparePredictors(dd, eyeData_rmotl_cat, t_r, param, catEvTimes);
                n=load(predictorInfoName, 'predictorInfo');
                predictorInfo = n.predictorInfo;
                n=[];
                %load(fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']));
            end
            

            if splitPredictor
                [predictorInfo, param] = splitPredictorByCue(predictorInfo, dd, onsets_cat, param);

                disp('fit kernels')
                [trIdx_r] = retrieveTrIdx_r(t_cat, t_r, t_tr);
                [predicted_all, predicted, PSTH_f, kernelInfo] = fitPSTH_cv(spk_all_cat, ...
                    predictorInfo.t_r, param.predictorNames,   predictorInfo.predictors_r, ...
                    predictorInfo.npredVars,param.psth_sigma, param.kernelInterval, ...
                    param.lagRange, param.ridgeParams, trIdx_r,fitoption,useGPU);

                y_r = cat(2,PSTH_f,predicted_all, predicted);

                %% figure for kernel fitting
                f = showKernel(t_r, y_r, kernelInfo, param.cardinalDir);
                screen2png(fullfile(saveFigFolder,['kernels_exp' saveSuffix '_wocue']), f);
                close(f);

                y_r_wcue = cat(2,PSTH_f,predicted_all, predicted(:,6:10));
                kernelInfo_wcue = kernelInfo;
                kernelInfo_wcue.kernel = kernelInfo.kernel(6:10);
                kernelInfo_wcue.tlags = kernelInfo.tlags(6:10);
                f = showKernel(t_r, y_r_wcue, kernelInfo_wcue, param.cardinalDir);
                screen2png(fullfile(saveFigFolder,['kernels_exp' saveSuffix '_wcue']), f);
                close(f);

                %% save results
                mm_s=matfile([saveName '_splitPredictor'],'writable',true);
                mm_s.PSTH_f = PSTH_f;
                mm_s.predicted_all = predicted_all;
                mm_s.predicted = predicted;
                mm_s.kernelInfo = kernelInfo;
                mm_s.t_r = t_r;
                mm_s.param=param;
                mm_s.predictorInfo = predictorInfo;

                clear mm_s predictorInfo param kernelInfo predicted predicted_all PSTH_f;
            end

            
            %% remove trials with too short duration 28/10/23
            [t_tr, catEvTimes, validTrials] = trimInvalids(t_tr, catEvTimes);


            if fitIt && ~splitPredictor
                %% obtain kernels!
                disp('fit kernels')
                [trIdx_r] = retrieveTrIdx_r(t_cat, t_r, t_tr);
                [predicted_all, predicted, PSTH_f, kernelInfo] = fitPSTH_cv(spk_all_cat, ...
                    predictorInfo.t_r, param.predictorNames,   predictorInfo.predictors_r, ...
                    predictorInfo.npredVars,param.psth_sigma, param.kernelInterval, ...
                    param.lagRange, param.ridgeParams, trIdx_r,fitoption,useGPU);
                
                % load(saveName, 'predicted_all','predicted','PSTH_f','kernelInfo');
                
                
                % %% extract time course of cue [0 1] for each target direction
                %   param.predictorNames = {'cue'};
                %   predictorInfo_cue = preparePredictors(dd, eyeData_rmotl_cat, t_r, param, catEvTimes);
                %
                % %% fit with ridgeX
                % [predicted_cue, gainInfo] = fitMultiplicative(PSTH_f, predicted_all, t_r, ...
                %     predictorInfo.predictors_r);
                
                y_r = cat(2,PSTH_f,predicted_all, predicted);
                
                %% figure for kernel fitting
                    f = showKernel( t_r, y_r, kernelInfo, param.cardinalDir);
                    screen2png(fullfile(saveFigFolder,['kernels_exp' saveSuffix]), f);
                    close(f);
                 
                %% Figure for target onset response (only to preferred direction)
                [f, cellclassInfo] = showTonsetResp(t_r, y_r, catEvTimes, dd, psthNames, ...
                    startSaccNoTask, saccDirNoTask, param, [-0.5 0.5]);
                cellclassInfo.datech = datech;
                %savePaperFigure(f, fullfile(saveFigFolder,['cellclassFig_' saveSuffix]));
                screen2png(fullfile(saveFigFolder,['cellclassFig_' saveSuffix '_allTr']), f);
                close(f);
                
                
                %% Figure for target onset, w/wo cue
                [f, avgTOnsetByCue, winSamps_tonsetByCue] = showTonsetByCue(t_r, ...
                    y_r, param.cardinalDir, catEvTimes, dd, psthNames, [-0.5 0.5], 1);
                screen2png(fullfile(saveFigFolder,['tonsetByCue_' saveSuffix '_onlySuccess']), f);
                close(f);
                
                [f, avgTOnsetByCue_parea, winSamps_tonsetByCue_parea] = showTonsetByCue(t_r, ...
                    predictorInfo.predictors_r(17,:)', param.cardinalDir, catEvTimes, dd, ...
                    {'parea'}, [-0.5 0.5]);
                set(f,'position',[0 0 1920 300]);
                screen2png(fullfile(saveFigFolder,['tonsetByCue_' saveSuffix '_parea']), f);
                close(f);
                
                %% response to saccades outside the task
                [f,avgSaccResp, winSamps_sacc, singleSaccResp, sortedSaccLabels] = ...
                    showSaccOnsetResp(t_r, y_r, param.cardinalDir, psthNames, ...
                    startSaccNoTask, saccDirNoTask, [-0.5 0.5]);
                screen2png(fullfile(saveFigFolder,['saccOn_' saveSuffix]));
                close(f);
                
                %% response to fixation and cue onsets
                [f, avgfOnsetResp, avgCueResp, winSamps_fc] = showFixCueOnsetResp(t_r, ...
                    y_r, catEvTimes, dd, psthNames, [-0.5 1]);
                screen2png(fullfile(saveFigFolder,['fixCueOnsetResp_' saveSuffix]),f);
                close(f);
                
                close all;
                
                %%response of pdiam to fixation and cue onsets
                [f, avgfOnsetResp_pdiam, avgCueResp_pdiam, winSamps_fc_pdiam] = showFixCueOnsetResp(t_r, ...
                    predictorInfo.predictors_r(17,:)', catEvTimes, dd, {'parea'}, [-0.5 1]);
                screen2png(fullfile(saveFigFolder,['fixCueOnsetResp_' saveSuffix '_pdiam']),f);
                close(f);
                

                %% obtain gain
                onlySuccess = 0;
                respWin = [0.05 0.35]; %[s] %time after stim onset to obtain preferred direction
                gainInfo = getGainInfo(t_r, y_r(:,1:2), param.cardinalDir, catEvTimes, ...
                    dd, [-0.5 1], onlySuccess, respWin);
                f=showGainInfo(gainInfo);
                savefigname = fullfile(saveFigFolder,[saveSuffix '_gainInfo']);
                screen2png(savefigname);
                close(f);
                


                %% save results
                
                % save(saveName, 'PSTH_f','predicted_all', 'predicted','kernelInfo'...
                %     ,'t_r','cellclassInfo','param','mFiringRate','t_cat','dd');
                % 'avgfOnsetResp', 'avgCueResp', 'winSamps_fc', ...
                % 'avgTOnsetByCue','winSamps_sacc', 'singleSaccResp', 'sortedSaccLabels',...
                %);
                %'pspec_psth','pspec_parea','faxis_psth','faxis_parea');
                %         clear spk_all dd kernel kernel_x kernel_y psth_all mDir seDir mDir_pred seDir_pred


                mm=matfile(saveName,'writable',true);
                mm.PSTH_f = PSTH_f;
                mm.predicted_all = predicted_all;
                mm.predicted = predicted;
                mm.kernelInfo = kernelInfo;
                mm.t_r = t_r;
                mm.cellclassInfo = cellclassInfo;
                mm.param=param;
                mm.mFiringRate=mFiringRate;
                mm.t_cat=t_cat;
                mm.dd=dd;
                mm.gainInfo = gainInfo; %26/10/2023
                
                clear mm;

                %previousDate = thisDate;

            end
            
        catch err
            disp(err);
            ng = [ng idata];
            close all;
        end
    end
end
