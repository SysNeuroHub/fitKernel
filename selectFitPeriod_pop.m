
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

kernelall_pop = [];
kernelsel_pop = [];
ich=1;
datechs = cell(1);
for idata = 1:length(channels) %1061;%865;%
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
    try
        load(saveName, 'kernelInfo','kernelInfo_selected', 'msaccTimeIdx','msaccParam');
    catch err
        continue;
    end
    if ~isempty(kernelInfo_selected)
        tlags = kernelInfo.tlags;
        kernelall_pop = cat(2, kernelall_pop, kernelInfo.kernel(:,17));
        kernelsel_pop = cat(2, kernelsel_pop, kernelInfo_selected.kernel);
        
        datechs{ich} = datech;
        ich=ich+1;
    end
end

subplot(311)
plot(tlags, mean(kernelall_pop,2), 'b',...
    tlags, mean(kernelsel_pop,2),'r');
title(['#ch: ' num2str(length(datechs)) ', lambda: ' num2str(msaccParam.lambda) ...
    ', minDur: ' num2str(msaccParam.minDur) ', minData: ' num2str(msaccParam.minData)]);
legend('all period', 'fixation period wo saccade');

subplot(312)
plot(tlags, kernelall_pop, 'color',[.5 .5 .5]);hold on;
plot(tlags, mean(kernelall_pop,2), 'b');

subplot(313)
plot(tlags, kernelsel_pop, 'color',[.5 .5 .5]);hold on;
plot(tlags, mean(kernelsel_pop,2), 'r');

xlabel('Time lag [s]');
screen2png(['compKernels_pop']);
close all;
