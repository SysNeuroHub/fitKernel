addpath(genpath('C:\Users\dshi0006\git'))
% load('C:\Users\dshi0006\Downloads\hugo_oephysdata_ch23.mat', ...
%     'ch','dd','ephysdata');

saveFolder = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\Daisuke\cuesaccade_data';

%% recorded data
animal = 'hugo';
rootFolder = '\\storage.erc.monash.edu.au\shares\R-MNHS-Physio\SysNeuroData\Monash Data\Joanita\2021/cuesaccade_data/';
dataType = 0;%0: each channel, 1: all channels per day


[loadNames, months, dates, channels] = getMonthDateCh(animal, rootFolder);

% to obtain index of specified month&date&channel
%find(1-cellfun(@isempty, regexp(loadNames, regexptranslate('wildcard','12December\15\*_ch26.mat'))))

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
param.predictorNames = {'vision','eyeposition','pdiam','blink'};
param.figTWin = [-0.5 0.5]; %temporal window for peri-event traces [s]
param.respWin = [0.03 0.25]; %temporal window to compute direction tuning
param.pareaTh = 3;
param.pareaDiffTh = 5;
param.cutoffFreq = 0.1;
param.evName = 'tOnset';%'cOnset';
cardinalDir = linspace(0,360,9); %direction for saccade and target
param.cardinalDir = cardinalDir(1:end-1);
ncDirs = length(param.cardinalDir);

previousDate = [];
for idata = 865%95:length(channels) %1061;%865;%
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
%         if ~strcmp(thisDate, previousDate)
%             [eyeData_cat, onsets_cat, meta_cat] = concatenate_eye(eyeData, dd);
%             t_tr={eyeData.t};
%             
%             
%             %% detect and interpolate blinks
%             %removing outliers helps for larger kernels
%             [eyeData_rmblk_cat, blinks] = removeBlinksEDF(eyeData_cat, meta_cat, param.marginSize,1);
%             close;
%             
%             %% remove parea outliers based on diff
%             [eyeData_rmotl_cat, outliers] = removePareaOutliers(eyeData_rmblk_cat, ...
%                 param.marginSize, param.pareaTh, param.pareaDiffTh);
%             screen2png(['rmOtl_' saveSuffix]);
%             close;
%             
%             clear eyeData_cat eyeData_rmblk_cat eyeData
%             
%             [pspec_parea,faxis_parea] = pmtm(eyeData_rmotl_cat.parea, 10, ...
%                 length(eyeData_rmotl_cat.parea), fs_eye);%slow
%             
%             save(fullfile(saveFolder,['eyeCat_' thisDate '.mat']), 'eyeData_rmotl_cat',...
%                 'onsets_cat','meta_cat','blinks','outliers','pspec_parea','faxis_parea','t_tr');
%             
%             %% prepare predictor variables
%             % make this part a function?
%              t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
%            
%             predictorInfo = preparePredictor(t_r, param.predictorNames, eyeData_rmotl_cat, ...
%                 onsets_cat, blinks, param);
% 
%             save(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');
%         else
            disp('loading eye/predictor data');
            load(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');
            load(fullfile(saveFolder,['eyeCat_' thisDate '.mat']), 'eyeData_rmotl_cat',...
                'onsets_cat','meta_cat','blinks','outliers','pspec_parea','faxis_parea','t_tr');
             t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
%         end
        
        
        
        %% select only fixation period
        theseTrials = intersect(find(~isnan(onsets_cat.fOnset)), find(~isnan(onsets_cat.tOnset)));
        fOnset_nonan = onsets_cat.fOnset(theseTrials);
        tOnset_nonan = onsets_cat.tOnset(theseTrials);
        trace_rising = zeros(length(t_r),1);
        trace_falling = zeros(length(t_r),1);
        for itr = 1:length(fOnset_nonan)
            [~, thisRise] = min(abs(t_r - fOnset_nonan(itr)));
            [~, thisFall] = min(abs(t_r - tOnset_nonan(itr)));
            if thisFall <= thisRise
                continue;
            end
            trace_rising(thisRise) = 1;
            trace_falling(thisFall) = 1;
        end
        includeIdx = find(cumsum(trace_rising)-cumsum(trace_falling) >= 1);
    
        predictors_s = interp1(t_r(includeIdx), predictorInfo.predictors_r(17,includeIdx),t_r)';
        
        %% obtain kernels!
        disp('fit kernels')
        [~, ~, kernelInfo_all] = fitPSTH(spk_all_cat, ...
            predictorInfo.t_r, predictorInfo.predictors_r(17,:), param.psth_sigma, param.lagRange, param.ridgeParams);
        
        [~, ~, kernelInfo_selected] = fitPSTH(spk_all_cat, ...
            predictorInfo.t_r, predictors_s, param.psth_sigma, param.lagRange, param.ridgeParams);
        
        
        plot(kernelInfo_all.tlags, kernelInfo_all.kernel, ...
            kernelInfo_selected.tlags, kernelInfo_selected.kernel);
        legend('all period', 'fixation period');
        xlabel('Time lag [s]');
       
        screen2png(['compKernels' saveSuffix]);
            
        
    
    end
end
