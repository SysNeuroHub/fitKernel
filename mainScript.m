
if ispc
    addpath(genpath('C:/Users/dshi0006/git'))
    saveFolder = 'E:/tmp/cuesaccade_data';
    saveFigFolder = [saveFolder, '/20220718'];
    mkdir(saveFigFolder);
    %saveFolder = '//storage.erc.monash.edu/shares/R-MNHS-Syncitium/Shared/Daisuke/cuesaccade_data';
    rootFolder = '//storage.erc.monash.edu.au/shares/R-MNHS-Physio/SysNeuroData/Monash Data/Joanita/2021/cuesaccade_data/';
elseif isunix
    addpath(genpath('/home/localadmin/Documents/MATLAB'));
    saveFolder = '/mnt/syncitium/Daisuke/cuesaccade_data';
    rootFolder = '/mnt/physio/Monash Data/Joanita/2021/cuesaccade_data/';
end

%% recorded data
animal = 'hugo';
dataType = 0;%0: each channel, 1: all channels per day

[loadNames, months, dates, channels] = getMonthDateCh(animal, rootFolder);

% to obtain index of specified month&date&channel
thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
    regexptranslate('wildcard','03March\16\*_ch29.mat'))));

%% omit data
% no saccade response
% low spontaneous firing
% low number of successful trials

% parameters
load('E:\tmp\cuesaccade_data\param20220719','param');
ncDirs = length(param.cardinalDir);

psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);

%ng = [];
previousDate = [];
for idata = 1:6%length(channels) %1:795
    
    try
        datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
        disp(datech);
        
        saveSuffix = [animal replace(datech,'/','_')];
        
        thisDate = [months{idata} '_' dates{idata}];
        if sum(strcmp(thisDate, {'06June_06','06June_11','06June_09'}))>0
            %june11  Sample points must be unique.
            %june09
            %june06 weird blank period in time around 500-600s
            continue;
        end
        
        
        load(loadNames{idata});
        
        if dataType == 0
            spk_all = ephysdata.spikes.spk;
            chName = ['_ch' num2str(ephysdata.spikes.chanIds)];
            clear ephysdata
        end
        
        
        if ~isempty(spk_all)
            
            
            nTrials = length(dd.eye);
            fs_eye = median([dd.eye.fs]);
            eyeData = dd.eye;
            
            %% concatenate across trials
            [spk_all_cat, t_cat] = concatenate_spk(spk_all, {dd.eye.t});
            clear spk_all
            mFiringRate = length(spk_all_cat)/(t_cat(end)-t_cat(1)); %spks/s
            if mFiringRate < 5
                disp([chName 'skipped as mFiringRate<5']);
                continue;
            end
            
            
            %% prepare behavioral data (common across channels per day)
            eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
            if  ~exist(eyeName, 'file')%~strcmp(thisDate, previousDate)
                
                [eyeData_rmotl_cat, catEvTimes, t_tr, onsets_cat,meta_cat,blinks,outliers] ...
                    = processEyeData(eyeData, dd, param);
                %             [pspec_parea,faxis_parea] = pmtm(eyeData_rmotl_cat.parea, 10, ...
                %                 length(eyeData_rmotl_cat.parea), fs_eye);%slow
                
                tOnset = catEvTimes.tOnset;
                cOnset = catEvTimes.cOnset;
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
                
                save(eyeName,'eyeData_rmotl_cat','catEvTimes',...
                    'onsets_cat','meta_cat','blinks','outliers','t_tr',...
                    'startSaccNoTask', 'endSaccNoTask', ...
                    'saccDirNoTask', 'dirIndexNoTask');
                
                %% prepare predictor variables
                t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
                predictorInfo = preparePredictors(dd, eyeData_rmotl_cat, t_r, param, catEvTimes);
                save(fullfile(saveFolder,['predictorInfo_' animal thisDate '.mat']), 'predictorInfo');
            else
                disp('loading eye/predictor data');
                load(fullfile(saveFolder,['predictorInfo_' animal thisDate '.mat']), 'predictorInfo');
                load(fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']));
                t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
            end
            
            
            %% obtain kernels!
            disp('fit kernels')
            [trIdx_r] = retrieveTrIdx_r(t_cat, t_r, t_tr);
            [predicted_all, predicted, PSTH_f, kernelInfo] = fitPSTH_cv(spk_all_cat, ...
                predictorInfo.t_r, param.predictorNames,   predictorInfo.predictors_r, predictorInfo.npredVars,...
                param.psth_sigma, param.kernelInterval, param.lagRange, param.ridgeParams, trIdx_r);
            
            y_r = cat(2,PSTH_f,predicted_all, predicted);
            
            %% figure for kernel fitting
            f = showKernel( t_r, y_r(:,1:2), kernelInfo, param.cardinalDir);
            screen2png(fullfile(saveFigFolder,['kernels_exp' saveSuffix]), f);
            close(f);
            
            %% Figure for target onset response
            [f, cellclassInfo] = showTonsetResp(t_r, y_r, catEvTimes, dd, psthNames, ...
                startSaccNoTask, saccDirNoTask, param, [-0.5 0.5]);
            cellclassInfo.datech = datech;
            screen2png(fullfile(saveFigFolder,['cellclassFig_' saveSuffix '_allTr']), f);
            close(f);
            
            %% Figure for target onset, w/wo cue
            [f, avgTOnsetByCue, winSamps_tonsetByCue] = showTonsetByCue(t_r, ...
                y_r, param.cardinalDir, catEvTimes, dd, psthNames, [-0.5 0.5]);
            screen2png(fullfile(saveFigFolder,['tonsetByCue_' saveSuffix]), f);
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
            
            %% save results
            saveName = fullfile(saveFolder, [saveSuffix '.mat']);
            
            save(saveName, 'PSTH_f','predicted_all', 'predicted','kernelInfo',...
                't_r','cellclassInfo','avgfOnsetResp', 'avgCueResp', 'winSamps_fc', ...
                'avgTOnsetByCue','param','winSamps_sacc', 'singleSaccResp', 'sortedSaccLabels','mFiringRate');
            %'pspec_psth','pspec_parea','faxis_psth','faxis_parea');
            %         clear spk_all dd kernel kernel_x kernel_y psth_all mDir seDir mDir_pred seDir_pred
            
            previousDate = thisDate;
            
        end
    catch err
        %ng = [ng idata];
        continue;
    end
end
