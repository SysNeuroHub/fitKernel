switch getenv('COMPUTERNAME')

    case 'MU00175834'
        addpath(genpath('C:/Users/dshi0006/git'))
        saveServer = 'E:/tmp/cuesaccade_data';
        %saveFigFolder = [saveServer, '/20220722'];
        %mkdir(saveFigFolder);
        rootFolder = '//storage.erc.monash.edu.au/shares/R-MNHS-Physio/SysNeuroData/Monash Data/Joanita/';

    case 'MU00011697'
        saveServer = '~/Documents/cuesaccade_data';
        rootFolder = '/mnt/MBI/Monash Data/Joanita/';
        
    case 'MU00108396'
        addpath(genpath('/home/localadmin/Documents/MATLAB'));
        saveFolder = '/mnt/syncitium/Daisuke/cuesaccade_data';
        rootFolder = '/mnt/physio/Monash Data/Joanita/2021/cuesaccade_data/';
end

%saveFigFolder = [saveFolder, '/20230504'];
%mkdir(saveFigFolder);

%% recorded data
animal = 'hugo'; %'m1899' 'andy' 'ollie' 
year = '2021';

saveFigFolder = fullfile(saveServer, '20230619',year,animal);
mkdir(saveFigFolder);


dataType = 0;%0: each channel, 1: all channels per day
fitIt = 1;

[loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);

% to obtain index of specified month&date&channel
% thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
%     regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','12December','14','*_ch13*')))));
thisdata = [];

if isempty(thisdata)
    thisdata = 1:length(channels);
end

%% omit data
% no saccade response
% low spontaneous firing
% low number of successful trials

% parameters
load(fullfile(saveServer,'param20230405.mat'),'param');
ncDirs = length(param.cardinalDir);
%param.lagRange(2,:)=[-1 0.5];

psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);

ng = [];
previousDate = [];
for idata = thisdata
    try
        % datech = [years{idata} filesep months{idata} filesep dates{idata} filesep num2str(channels{idata})];
        datech = [months{idata} filesep dates{idata} filesep num2str(channels{idata})];
        disp(datech);

        saveSuffix = [animal replace(datech,filesep,'_') ];%'_cue'];

        thisDate = [months{idata} '_' dates{idata}];
        if sum(strcmp(thisDate, {'06June_06','06June_11','06June_09'}))>0
            %june11  Sample points must be unique.
            %june09
            %june06 weird blank period in time around 500-600s
            continue;
        end
        saveFolder = fullfile(saveServer, year,animal);%17/6/23
        if ~exist(saveFolder, 'dir')
            mkdir(saveFolder);
        end
        saveName = fullfile(saveFolder, [saveSuffix '.mat']);

        load(loadNames{idata},'ephysdata','dd');
        
        spk_all = ephysdata.spikes.spk;
        if ~isempty(spk_all)
            %% concatenate across trials
            [spk_all_cat, t_cat] = concatenate_spk(spk_all,  dd.eye);
            clear spk_all
            mFiringRate = length(spk_all_cat)/(t_cat(end)-t_cat(1)); %spks/s
        else
            mFiringRate = 0;
        end
        clear ephysdata

        if mFiringRate < 5
            disp(['skipped as mFiringRate<5']);
            save(saveName,'mFiringRate');
            continue;
        end


        %% prepare behavioral data (common across channels per day)
        eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
        if  ~exist(eyeName, 'file') %~strcmp(thisDate, previousDate)

            [eyeData_rmotl_cat, catEvTimes, t_tr, onsets_cat,meta_cat,blinks,outliers] ...
                = processEyeData(dd.eye, dd, param);
            %             [pspec_parea,faxis_parea] = pmtm(eyeData_rmotl_cat.parea, 10, ...
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

            save(eyeName,'eyeData_rmotl_cat','catEvTimes',...
                'onsets_cat','meta_cat','blinks','outliers','t_tr',...
                'startSaccNoTask', 'endSaccNoTask', ...
                'saccDirNoTask', 'dirIndexNoTask');
            close all
            %% prepare predictor variables
            t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
            predictorInfo = preparePredictors(dd, eyeData_rmotl_cat, t_r, param, catEvTimes);
            save(fullfile(saveFolder,['predictorInfo_' animal thisDate '.mat']), 'predictorInfo');
        else
            disp('loading eye/predictor data');
            %load(fullfile(saveFolder,['predictorInfo_' animal thisDate '.mat']), 'predictorInfo');
            load(eyeName,'eyeData_rmotl_cat','catEvTimes',...
                'onsets_cat','meta_cat','blinks','outliers','t_tr',...
                'startSaccNoTask', 'endSaccNoTask', ...
                'saccDirNoTask', 'dirIndexNoTask');
            t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
            predictorInfo = preparePredictors(dd, eyeData_rmotl_cat, t_r, param, catEvTimes);
            save(fullfile(saveFolder,['predictorInfo_' animal thisDate '.mat']), 'predictorInfo');
            load(fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']));
         end


        if fitIt
            %% obtain kernels!
            disp('fit kernels')
            [trIdx_r] = retrieveTrIdx_r(t_cat, t_r, t_tr);
            [predicted_all, predicted, PSTH_f, kernelInfo] = fitPSTH_cv(spk_all_cat, ...
                predictorInfo.t_r, param.predictorNames,   predictorInfo.predictors_r, predictorInfo.npredVars,...
                param.psth_sigma, param.kernelInterval, param.lagRange, param.ridgeParams, trIdx_r);

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



            %% Figure for target onset response
            [f, cellclassInfo] = showTonsetResp(t_r, y_r, catEvTimes, dd, psthNames, ...
                startSaccNoTask, saccDirNoTask, param, [-0.5 0.5]);
            cellclassInfo.datech = datech;
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

            %% save results

            save(saveName, 'PSTH_f','predicted_all', 'predicted','kernelInfo','t_r','cellclassInfo','param','mFiringRate','t_cat');
                % 'avgfOnsetResp', 'avgCueResp', 'winSamps_fc', ...
                % 'avgTOnsetByCue','winSamps_sacc', 'singleSaccResp', 'sortedSaccLabels',...
                %);
            %'pspec_psth','pspec_parea','faxis_psth','faxis_parea');
            %         clear spk_all dd kernel kernel_x kernel_y psth_all mDir seDir mDir_pred seDir_pred

            %previousDate = thisDate;
        end

    catch err
        disp(err);
        ng = [ng idata];
        close all;
    end
end