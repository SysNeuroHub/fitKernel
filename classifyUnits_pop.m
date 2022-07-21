%created from resp_cueConditions.m

if ispc
    addpath(genpath('C:/Users/dshi0006/git'))
    saveFolder = 'E:/tmp/cuesaccade_data';
    rootFolder = '//storage.erc.monash.edu.au/shares/R-MNHS-Physio/SysNeuroData/Monash Data/Joanita/2021/cuesaccade_data/';
    saveFigFolder = '//storage.erc.monash.edu/shares/R-MNHS-Syncitium/Shared/Daisuke/cuesaccade_data/20220617/';
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
    regexptranslate('wildcard','09September\01\*_ch27.mat'))));

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
        

psthIdx = 1;
allMdlIdx = 2;
visionIdx = 3;
eyevelIdx = 4;

ii=1;
for idata = 1:length(channels) 
    datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
    disp(datech);
    
    saveSuffix = [animal replace(datech,'/','_')];
    
    thisDate = [months{idata} '_' dates{idata}];
   
    saveName = fullfile(saveFolder, ['cellclassInfo_' saveSuffix '.mat']);
    if exist(saveName,'file')>0
        load(saveName, 'cellclassInfo');
        cellclassInfo_pop(ii) = cellclassInfo;
        datech_pop{ii} = datech;
        ii = ii+1;
        clear cellclassInfo;
    end
end


unitClass = [cellclassInfo_pop.unitClass];
[nunits,edges]=histcounts(unitClass, 0.5+0:6);

mtOnsetRespAmp = [];
nmtOnsetRespAmp = [];
msaccamps = [];
% load(loadName,'winSamps_sacc')
% winSamps = winSamps_sacc;
for ii = 1:numel(cellclassInfo_pop)
    %         tmp = cellclassInfo_pop(ii).mtOnsetResp([psthIdx visionIdx eyevelIdx],:);
    %         tmp = reshape(tmp, [1 size(tmp,1) size(tmp,2)]);
    %         amps = characteriseResp(tmp, ...
    %             winSamps, param.tOnRespWin, param.baseWin);
    amps = cellclassInfo_pop(ii).mtOnsetRespAmp;
    nmtOnsetRespAmp(ii,:) = [amps(2)/amps(1) amps(3)/amps(1)];
    mtOnsetRespAmp(ii) = amps(1);
    
    %tmp = cellclassInfo_pop(ii).msaccResp(psthIdx,:);
    %tmp = reshape(tmp, [1 size(tmp,1) size(tmp,2)]);
    %saccamps(ii) = characteriseResp(tmp, ...
    %    winSamps, param.tOnRespWin, param.baseWin);
    msaccamps(ii) = cellclassInfo_pop(ii).msaccRespAmp;
    
%     unitClass(ii) = getCellClass(PtonsetResp, PtonsetResp_paired, ...
%         mtOnsetRespAmp, PsaccResp, Pth);
    
end

figure('position',[0 0 1700 500]);
subplot(131);
histaxis = -8:1:15;%0:2:40;%0:5:60;
histogram(mtOnsetRespAmp,histaxis);
hold on
histogram(mtOnsetRespAmp(unitClass<=3), histaxis);
vline(0);
xlabel('response to tOnset [Hz]');
ylabel('#units')

subplot(132);
hist3(nmtOnsetRespAmp(unitClass<4,:),{-0.5:0.1:2,-0.5:0.1:2},...
    'edgecolor','interp','CdataMode','auto');
colormap(gca,1-gray);
view(2)
axis equal tight;
squareplot
set(gca,'tickdir','out');
[~,ax] = mcolorbar(gca,0.5);
ylabel(ax,'#units');
grid off
xlabel('vision response');ylabel('eye velocity response');

subplot(133);
pie(nunits,{'visual dominant','eye dominant','visual+eye','eye outside task','not responsive'});
colormap(gca,'default');

saveas(gcf, 'classifyUnits_pop.fig');
screen2png('classifyUnits_pop');
