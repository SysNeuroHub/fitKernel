
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
    regexptranslate('wildcard','12December\09\*_ch13.mat'))));

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
for idata = 390:length(channels) %1061;%865;%
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
    
    saveName = fullfile(saveFolder, [saveSuffix '.mat']);

    load(loadNames{idata},'ephysdata','dd');
    
    if dataType == 0
        spk_all = ephysdata.spikes.spk;
        chName = ['_ch' num2str(ephysdata.spikes.chanIds)];
        clear ephysdata
    end
    
    
    if ~isempty(spk_all)
%        load(saveName, 'mFiringRate');
%         nTrials = length(dd.eye);
%         fs_eye = median([dd.eye.fs]);
%         eyeData = dd.eye;
%         
%         %% concatenate across trials
         [spk_all_cat, t_cat] = concatenate_spk(spk_all, {dd.eye.t});
         clear spk_all eyeData dd
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
        
        
        %% detect microsaccade Engbert and Kliegel vis res
       msaccParam.lambda=6;
       msaccParam.minDur = 0.005;%[s] Engbert used 0.012
       %NG 0.003 
       %OK 0.009 0.012 with 1% cut off but kernel is not too different
       msaccParam.minData = 150;
       
       msaccTimeIdx = detectSaccTimes(eyeData_rmotl_cat.x,eyeData_rmotl_cat.y,...
           eyeData_rmotl_cat.dt,msaccParam.lambda, msaccParam.minDur);
        
          %% select only fixation period          
        theseTrials = intersect(find(~isnan(onsets_cat.fOnset)), find(~isnan(onsets_cat.cOnset)));
        fOnset_nonan = onsets_cat.fOnset(theseTrials);
        cOnset_nonan = onsets_cat.cOnset(theseTrials); %12/2/22
        fixTrace = event2Trace(predictorInfo.t_r, [fOnset_nonan cOnset_nonan-0.5]);
        saccTrace = event2Trace(predictorInfo.t_r, [catEvTimes.saccadeStartTimes catEvTimes.saccadeEndTimes], 1);
        blinkTrace = event2Trace(predictorInfo.t_r, [catEvTimes.blinkStartTimes catEvTimes.blinkEndTimes], 1);
        otlTrace = event2Trace(predictorInfo.t_r, [catEvTimes.outlierStartTimes catEvTimes.outlierEndTimes], 1);
        msaccTrace = event2Trace(predictorInfo.t_r, [eyeData_rmotl_cat.t(msaccTimeIdx(:,1)) eyeData_rmotl_cat.t(msaccTimeIdx(:,2))],1);
        includeIdx = intersect(find(fixTrace),  find(1-saccTrace));
        includeIdx = intersect(includeIdx, find(1-blinkTrace));
        includeIdx = intersect(includeIdx, find(1-otlTrace));
        includeIdx = intersect(includeIdx, find(1-msaccTrace));
        
        if length(includeIdx) <  msaccParam.minData %3s
            disp(num2str(length(includeIdx)));
            disp('too few data for fitting');
            kernelInfo_selected = [];
            save(saveName, 'kernelInfo_selected', 'msaccTimeIdx','msaccParam','-append');
            continue;
        end
        %this period can still include saccade - after fixation cue start
        %till pupil focus on the cue
        predictors_s = interp1(predictorInfo.t_r(includeIdx), ...
            predictorInfo.predictors_r(17,includeIdx),predictorInfo.t_r)';
        
%         predictors_s = predictorInfo.predictors_r(17,:);
%         excludeIdx = setdiff(1:length(predictorInfo.t_r), includeIdx);
%         predictors_s(excludeIdx)=nan;
        
%         plot(predictorInfo.t_r, predictorInfo.predictors_r(17,:));%original pdiam over time
%         hold on
%         plot(predictorInfo.t_r, predictors_s, '.'); %only around fixation onset
%         vbox(catEvTimes.saccadeStartTimes, catEvTimes.saccadeEndTimes,gca);
        
        
        %% obtain kernels!
        disp('fit kernels')
        %param.lagRange = repmat([-0.5 0.5],[18,1]);
        %param.lagRange(1:8,1)=0;
        [predicted_all, PSTH_f, kernelInfo_all] = fitPSTH(spk_all_cat, ...
            predictorInfo.t_r, predictorInfo.predictors_r(17,:), param.psth_sigma, param.lagRange, param.ridgeParams);

        [~, ~, kernelInfo_selected] = fitPSTH(spk_all_cat, ...
            predictorInfo.t_r, predictors_s, param.psth_sigma, param.lagRange, param.ridgeParams);

        plot(kernelInfo_all.tlags, kernelInfo_all.kernel, ...
            kernelInfo_selected.tlags, kernelInfo_selected.kernel);
        legend('all period', 'fixation period');
        xlabel('Time lag [s]');
        screen2png(['compKernels' saveSuffix]);
close all;
        
        %% save results
        
        save(saveName, 'kernelInfo_selected', 'msaccTimeIdx','msaccParam','-append');
            %'pspec_psth','pspec_parea','faxis_psth','faxis_parea');
        clear spk_all dd kernel kernel_x kernel_y psth_all mDir seDir mDir_pred seDir_pred
        
        
        previousDate = thisDate;
        
    end
end
