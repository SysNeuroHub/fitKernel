
if ispc
addpath(genpath('C:/Users/dshi0006/git'))
saveFolder = '//storage.erc.monash.edu/shares/R-MNHS-Syncitium/Shared/Daisuke/cuesaccade_data';
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
    regexptranslate('wildcard','12December\09\*_ch1.mat'))));

%% omit data
% no saccade response
% low spontaneous firing
% low number of successful trials

% parameters
param.marginSize = 50;%40; %frames
param.psth_sigma = .05;%0.1; %0.025;%0.02;%0.01;%0.05;%[s] %gaussian filter
param.dt_r = 0.02; %for securing memory for kernel fitting %0.025;%0.5;
param.lagRange = [-.5 .5];%[-.4 .3];%[-1 1];%[-10 20]; %temporal window for kernel estimation [s]
param.ridgeParams = 100;%[0 1e-1 1 1e2 1e3]; %10
% visualize = 0;
param.predictorNames = {'vision','eyespeed','pdiam','blink','reward'}; %eyeposition
param.figTWin = [-0.5 0.5]; %temporal window for peri-event traces [s]
param.respWin = [0.03 0.25]; %temporal window to compute direction tuning
param.pareaTh = 3;
param.pareaDiffTh = 5;
param.cutoffFreq = 0.1;
param.evName = 'tOnset';%'cOnset';
%param.minSaccInterval = 0.5; %29/1/22
cardinalDir = linspace(0,360,9); %direction for saccade and target
param.cardinalDir = cardinalDir(1:end-1);
ncDirs = length(param.cardinalDir);

previousDate = [];
for idata = 365:length(channels) %1061;%865;%
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
        eyeName = fullfile(saveFolder,['eyeCat_' thisDate '.mat']);
        
        disp('loading eye/predictor data');
        load(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');
        load(fullfile(saveFolder,['eyeCat_' thisDate '.mat']));
        
        
        saveName = fullfile(saveFolder, [saveSuffix '.mat']);
        load(saveName, 'PSTH_f','kernelInfo','predicted_all','predicted');
        
        [~, pdiam] = getPupilDiameter(eyeData_rmotl_cat);
        pdiam_r = interp1(eyeData_rmotl_cat.t, pdiam, predictorInfo.t_r)';
          
        
        %% triggered by cueOnsets
        cOnset = catEvTimes.cOnset;
        validEvents = find(~isnan(cOnset));
        cOnsetTimes = cOnset(validEvents);
        tgtDir = getTgtDir(dd.targetloc(validEvents), param.cardinalDir);
        
        [avgConsetResp, winSamps_conset, singleConsetResp, sortedConsetLabels, uniqueConsetLabels] ...
            = eventLockedAvg(cat(1,PSTH_f',predicted_all, predicted, pdiam_r), ...
            predictorInfo.t_r, cOnsetTimes, tgtDir, param.figTWin);
         psthNames = cat(2,{'psth','predicted_all'},param.predictorNames,'pdiam ori');
        
        nvars = size(avgConsetResp,2);
        figure('position',[0 0 400 1000]);
        for ivar = 1:nvars
            subplot(nvars, 1, ivar);
            imagesc(winSamps_conset, param.cardinalDir, squeeze(avgConsetResp(:,ivar,:)));
            set(gca, 'ytick',param.cardinalDir);
            ylabel(psthNames{ivar});
            mcolorbar(gca,.5);
        end
        xlabel('time from cOnset [s]');
        screen2png(['cOnset_' saveSuffix]);
        close all;
        
        
        
        %% detect microsaccade Engbert and Kliegel vis res
       msaccParam.lambda=6;
       msaccParam.minDur = 0.006;%[s] Engbert used 0.012
       %msaccParam.minData = 150;
       
       msaccTimeIdx = detectSaccTimes(eyeData_rmotl_cat.x,eyeData_rmotl_cat.y,...
           eyeData_rmotl_cat.dt,msaccParam.lambda, msaccParam.minDur);
       msaccTimes = [eyeData_rmotl_cat.t(msaccTimeIdx(:,1)) eyeData_rmotl_cat.t(msaccTimeIdx(:,2))]; 
       
       
        %% response to micro-saccades (outside of task)
        tOnset = catEvTimes.tOnset;
        cOnset = catEvTimes.cOnset;
        validEvents = intersect(find(~isnan(tOnset)), find(~isnan(cOnset)));
        tOnset_val = tOnset(validEvents);
        cOnset_val = cOnset(validEvents);
        
        barrierDur = 2;% [s] 1
        tOnset_trace = event2Trace(t_cat, tOnset(~isnan(tOnset)), barrierDur); %19/2/22
        cOnset_trace = event2Trace(t_cat, cOnset(~isnan(cOnset)), barrierDur); %19/2/22
        tcOnset_trace = event2Trace(t_cat, [tOnset_val; cOnset_val], barrierDur); %19/2/22
        blink_trace = event2Trace(t_cat, [catEvTimes.blinkStartTimes; catEvTimes.blinkEndTimes],barrierDur);
        outlier_trace = event2Trace(t_cat, [catEvTimes.outlierStartTimes; catEvTimes.outlierEndTimes],barrierDur);
        
        excEventT_cat = (tOnset_trace + cOnset_trace + tcOnset_trace + blink_trace + outlier_trace > 0); %28/1/22
        
        [startSaccNoTask, endSaccNoTask] = selectSaccades(msaccTimes(:,1), ...
            msaccTimes(:,2), t_cat, excEventT_cat,[], diff(param.figTWin));%param.minSaccInterval);
        
        [saccDirNoTask, dirIndexNoTask] = getSaccDir(startSaccNoTask, endSaccNoTask, ...
            eyeData_rmotl_cat, param.cardinalDir);
        
       
        [avgSaccResp, winSamps_sacc, singleSaccResp, sortedSaccLabels, uniqueSaccLabels] ...
            = eventLockedAvg(cat(1,PSTH_f',predicted_all, predicted, pdiam_r), ...
            predictorInfo.t_r, startSaccNoTask, saccDirNoTask, param.figTWin);
        
        nvars = size(avgSaccResp,2);
        figure('position',[0 0 400 1000]);
        for ivar = 1:nvars
            subplot(nvars, 1, ivar);
            imagesc(winSamps_sacc, param.cardinalDir, squeeze(avgSaccResp(:,ivar,:)));
            set(gca, 'ytick',param.cardinalDir);
            ylabel(psthNames{ivar});
            mcolorbar(gca,.5);           
        end
        xlabel('time from saccade onset [s]');
            
        
        screen2png(['./saccOn_' saveSuffix]);
        close all;
        
        
        
        %% save results
         
        save(saveName, 'msaccParam','msaccTimes', 'avgSaccResp', 'winSamps_sacc', ...
            'singleSaccResp', 'sortedSaccLabels', 'uniqueSaccLabels', ...
            'avgConsetResp', 'winSamps_conset', 'singleConsetResp', ...
            'sortedConsetLabels', 'uniqueConsetLabels','-append');
        clear spk_all dd kernel kernel_x kernel_y psth_all mDir seDir mDir_pred seDir_pred
        
        previousDate = thisDate;
        
    end
end
