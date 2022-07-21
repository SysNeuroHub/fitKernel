%created from resp_cueConditions.m

if ispc
    addpath(genpath('C:/Users/dshi0006/git'))
    rootFolder = '//storage.erc.monash.edu.au/shares/R-MNHS-Physio/SysNeuroData/Monash Data/Joanita/2021/cuesaccade_data/';
    %     saveFolder = '//storage.erc.monash.edu/shares/R-MNHS-Syncitium/Shared/Daisuke/cuesaccade_data';
    %     saveFigFolder = '//storage.erc.monash.edu/shares/R-MNHS-Syncitium/Shared/Daisuke/cuesaccade_data/20220617/';
    saveFolder = 'E:/tmp/cuesaccade_data';
    saveFigFolder = 'E:/tmp/cuesaccade_data/20220622/';
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
    regexptranslate('wildcard','09September\03\*_ch31.mat'))));

%% omit data
% no saccade response
% low spontaneous firing
% low number of successful trials

% parameters
param.th_frate = 5;
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
param.times = 1e-3*(0:30:390);%[0 50 100 150 200 250 300 350 400];

%
param.tOnRespWin = [0 0.5];
param.baseWin = [-0.1 0];
param.Pth = 0.01; %0.05

psthIdx = 1;
allMdlIdx = 2;
visionIdx = 3;
eyevelIdx = 4;

previousDate = [];
for idata = thisdata%1:length(channels) %1061;%865;%
    datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
    disp(datech);
    
    saveSuffix = [animal replace(datech,'/','_')];
    
    thisDate = [months{idata} '_' dates{idata}];
   
    loadName = fullfile(saveFolder, [saveSuffix '.mat']);
    saveName = fullfile(saveFolder, ['cellclassInfo_' saveSuffix '.mat']);
    saveFigName = fullfile(saveFigFolder, ['cellclassFig_' saveSuffix '.fig']);
    if exist(loadName,'file')>0
            
       load(loadNames{idata},'dd');
         
        %% prepare behavioral data (common across channels per day)
         
        disp('loading eye/predictor data');
        load(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');
        load(fullfile(saveFolder,['eyeCat_' thisDate '.mat']),'catEvTimes',...
            'blinks','outliers','eyeData_rmotl_cat');
        t_cat = eyeData_rmotl_cat.t;
        
        
        %% load neural data 
        load(loadName, 'PSTH_f','predicted_all','predicted');
        
        
        %% triggered by tOnsets
         onset = catEvTimes.(param.evName);
        
       %validEvents = intersect(find(~isnan(onset)), find(dd.cueOn==icue-1));
        validEvents = intersect(find(~isnan(onset)), find(~isnan(catEvTimes.cOnset)));
        %only use trials when the choices were registered.
        %this is a temporary fix as my current algorithm assumes stimuli were NOT
        %presented, causing no visual response in the model
        
        onsetTimes = onset(validEvents);
        tgtDir = getTgtDir(dd.targetloc(validEvents), param.cardinalDir);
        
        [~,dirIdx] = intersect(param.cardinalDir, unique(tgtDir));
        
        [avgOnsetResp, winSamps, singleOnsetResp, ...
            sortedOnsetLabels, uniqueOnsetLabels] ...
            = eventLockedAvg(cat(1,PSTH_f', predicted_all, predicted), ...
            predictorInfo.t_r, onsetTimes, tgtDir, param.figTWin);
     
       
        respTidx = intersect(find(param.tOnRespWin(1)<=winSamps), ...
            find(param.tOnRespWin(2)>=winSamps));
               
       
        tonsetRespAmp = characteriseResp(singleOnsetResp, ...
            winSamps, param.tOnRespWin, param.baseWin, 'mean');
        %trials x kinds
        
        prefDir = 180/pi*circ_mean(tgtDir'/180*pi, tonsetRespAmp(:,psthIdx));
        [~, prefDirIdx] = min(abs(circ_dist(pi/180*param.cardinalDir, pi/180*prefDir)));
        
        [PtonsetResp] = signrank(tonsetRespAmp(:,psthIdx));
        %[~,prefDir] = max(mean(avgOnsetResp(:,psthIdx,respTidx),3));%from alignMtxDir.m
        theseTrials = find(tgtDir == param.cardinalDir(prefDirIdx));%trials with cell's preferred direction
        %if numel(theseTrials)<=1
        %    continue; %cannot do the following stats
        %end
        
        [PtonsetResp_paired] = signrank(tonsetRespAmp(theseTrials,visionIdx), tonsetRespAmp(theseTrials,eyevelIdx));
        
        
        
        mtOnsetResp = squeeze(mean(singleOnsetResp(theseTrials,:,:)));%avg response to preferred direction
        setOnsetResp = 1/sqrt(numel(theseTrials))*squeeze(std(singleOnsetResp(theseTrials,:,:)));%avg response to preferred direction
 
        mtOnsetRespAmp = mean(tonsetRespAmp(theseTrials,[psthIdx visionIdx eyevelIdx]),1);
        
        
        %% triggered by saccade onsets (outside of the task)
        %from fitPSTH_test.m
        %load(saveName, 'singleSaccResp'); %event x saccade direction x time
        tOnset = catEvTimes.tOnset;
        cOnset = catEvTimes.cOnset;
        validEvents = intersect(find(~isnan(tOnset)), find(~isnan(cOnset)));
        tOnset = tOnset(validEvents);
        cOnset = cOnset(validEvents);
        
        tcOnset_trace = event2Trace(t_cat, [tOnset; cOnset], 2*0.5);
        excEventT_cat = (tcOnset_trace + blinks + outliers > 0); %28/1/22
        
        test = load(fullfile(saveFolder,['eyeCat_' thisDate '.mat']), 'startSaccNoTask');
        if isempty(fieldnames(test))
            [startSaccNoTask, endSaccNoTask] = selectSaccades(catEvTimes.saccadeStartTimes, ...
                catEvTimes.saccadeEndTimes, t_cat, excEventT_cat);%param.minSaccInterval);
            %<slow
            
            [saccDirNoTask, dirIndexNoTask] = getSaccDir(startSaccNoTask, endSaccNoTask, ...
                eyeData_rmotl_cat, param.cardinalDir);
            %<slow
            
            save(fullfile(saveFolder,['eyeCat_' thisDate '.mat']), 'startSaccNoTask', 'endSaccNoTask', ...
                'saccDirNoTask', 'dirIndexNoTask','-append');
        else
            load(fullfile(saveFolder,['eyeCat_' thisDate '.mat']), 'startSaccNoTask', 'endSaccNoTask', ...
                'saccDirNoTask', 'dirIndexNoTask');
        end
        
        [avgSaccResp, ~, singleSaccResp, sortedSaccLabels, uniqueSaccLabels] ...
            = eventLockedAvg(cat(1,PSTH_f',predicted_all, predicted), ...
            predictorInfo.t_r, startSaccNoTask, saccDirNoTask, param.figTWin);
        
        %use the same saccade direction to the one used for tOnset
        theseSaccTrials = (saccDirNoTask == param.cardinalDir(prefDirIdx));
        msaccResp = squeeze(nanmean(singleSaccResp(theseSaccTrials,:,:)));%avg response to preferred direction
         sesaccResp = 1/sqrt(numel(theseSaccTrials))*squeeze(nanstd(singleSaccResp(theseSaccTrials,:,:),1));
        
         saccRespAmp = characteriseResp(singleSaccResp, ...
            winSamps, param.tOnRespWin, param.baseWin, 'mean');
         msaccRespAmp = mean(saccRespAmp(theseSaccTrials,psthIdx),1);
       [PsaccResp] = signrank(saccRespAmp(:,psthIdx));
   
        
        %% categorise the cell
        unitClass = getCellClass(PtonsetResp, PtonsetResp_paired, mtOnsetRespAmp(2:3), ...
            PsaccResp, param.Pth);
        
        
        %% save results
        cellclassInfo.unitClass = unitClass;
        cellclassInfo.PtonsetResp = PtonsetResp;
        cellclassInfo.PtonsetResp_paired = PtonsetResp_paired;
        cellclassInfo.mtOnsetRespAmp = mtOnsetRespAmp;
        cellclassInfo.mtOnsetResp = mtOnsetResp;
        cellclassInfo.npreftonsetTrials = numel(theseTrials);
        cellclassInfo.PsaccResp = PsaccResp;
        cellclassInfo.msaccResp = msaccResp;
        cellclassInfo.msaccRespAmp = msaccRespAmp;
        cellclassInfo.nprefSaccTrials = numel(theseSaccTrials);
        cellclassInfo.winSamps = winSamps;
        cellclassInfo.datech = datech;
        
        save(saveName, 'cellclassInfo');
        
        medianSaccDelay = nanmedian(catEvTimes.cOnset-catEvTimes.tOnset);
        
        %% visualize the result
        figure('position',[0 0 1000 1000]);
        ax(1) = subplot(211);
        %plot(winSamps, mtOnsetResp([allMdlIdx visionIdx eyevelIdx],:));hold on;
        boundedline(winSamps, mtOnsetResp(psthIdx,:), setOnsetResp(psthIdx,:),'k', 'linewidth',2);
        hold on;
        boundedline(winSamps, mtOnsetResp(allMdlIdx,:), setOnsetResp(allMdlIdx,:),'b', 'transparency', 0.5);
        boundedline(winSamps, mtOnsetResp(visionIdx,:), setOnsetResp(visionIdx,:),'r', 'transparency', 0.5);
        boundedline(winSamps, mtOnsetResp(eyevelIdx,:), setOnsetResp(eyevelIdx,:),'c', 'transparency', 0.5);
        %vbox(param.baseWin(1), param.baseWin(2))
        %vbox(param.tOnRespWin(1), param.tOnRespWin(2),[],[.7 1 .7]);
        vline(medianSaccDelay);
        ylabel('tOnset');   
        title(num2str(unitClass));
        
        ax(2) = subplot(212);
        boundedline(winSamps, msaccResp(psthIdx,:), sesaccResp(psthIdx,:),'k', 'linewidth',2);
        hold on;
        boundedline(winSamps, msaccResp(allMdlIdx,:), sesaccResp(allMdlIdx,:),'b', 'transparency', 0.5);
        boundedline(winSamps, msaccResp(visionIdx,:), sesaccResp(visionIdx,:),'r', 'transparency', 0.5);
        boundedline(winSamps, msaccResp(eyevelIdx,:), sesaccResp(eyevelIdx,:),'c', 'transparency', 0.5);
        linkaxes(ax(:),'y');
        %vbox(param.baseWin(1), param.baseWin(2))
        %vbox(param.tOnRespWin(1), param.tOnRespWin(2),[],[.7 1 .7]);
        vline(0);
        xlim([winSamps(1)-medianSaccDelay winSamps(end)-medianSaccDelay])
        ylabel('saccade(outside task)');
        legend('observed','all mdl','vision','eye velocity','location','northwest')
        
        saveas(gcf,saveFigName);
        close;
        
        previousDate = thisDate;
        
    end
end
