
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

constPrefDir = 1; %11/3/22

ii=1;
avgOnsetResp_pop=[]; theseDatech = [];
for idata = 1:length(channels) 
    datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
    disp(datech);
    
    saveSuffix = [animal replace(datech,'/','_')];
      saveName = fullfile(saveFolder, [saveSuffix '.mat']);
      
    %% load results
    try
        load(saveName, 'avgOnsetResp','fitPars','fitErr');        
        avgOnsetResp_pop(:,:,:,:,ii) = avgOnsetResp;
        % target direction x varType x time x wo/w cue x trigger type x channels
        fitPars_pop(:,:,:,:,ii) = fitPars;
        % param x varType x time x w/wo cue
        %param: [Dp(peak angle), Rp(peak amp), Rn(2nd peak amp = 0), Ro(baseline), sigma(width)]
        fitErr_pop(:,:,:,ii) = fitErr;
        
        theseDatech{ii} = datech;
        ii=ii+1;
    catch err
        continue;
    end
end

param.times = 1e-3*(0:30:390);
taxis = param.times(1:end-1)+.5*mean(diff(param.times));
psthNames = cat(2,{'psth','psth_f','predicted_all'},...
    param.predictorNames,'pdiam ori');
winSamps = winSamps_conset;
tgtTimes = intersect(find(winSamps>0.03), find(winSamps<0.25));
% times = 1e-3*[0 50 100 150 200 250 300 350 400];

save('resp_cueConditions_pop.mat','avgOnsetResp_pop','fitPars_pop','fitErr_pop',...
    'theseDatech','taxis','winSamps','tgtTimes','psthNames');

%% fitting parameter w/wo cue
fitPars_pop = fitPars_pop([1 2 4 5],:,:,:,:);
for ipsth = 1:9
    figure('position',[0 0 1400 1000]);
    for ipar=1:4
        switch ipar
            case 1
                parName = 'Dp(peak angle)';
            case 2
                parName = 'Rp(peak amp)';
            case 3
                parName = 'Ro(baseline)';
            case 4
                parName = 'sigma(width)';
        end
        for it = 1:length(taxis)
            woCue = squeeze(fitPars_pop(ipar, ipsth, it, 1,:));
            wCue = squeeze(fitPars_pop(ipar, ipsth, it, 2,:));
            subplot(4,length(taxis),(ipar-1)*length(taxis)+it);
            plot(woCue, wCue,'.');
            if ipar==4
                xlabel('woCue');
            end
            if it==1
                ylabel([parName ' wCue']);
            end
            if ipar==1
                title(taxis(it));
            end
        end
    end
    squareplots;
    marginplots;
    screen2png(['pars' parName '_' psthNames{ipsth}]);
    close;
end

%% triggered by c/tOnsets
for itrig = 1:2
    switch itrig
        case 1
            onsetName = 'cOnset';
        case 2
            onsetName = 'tOnset';
    end
    
    
  
    figure('position',[0 0 1400 1000]);
    for icue = 1:3
        nvars = size(avgOnsetResp_pop,2);
        for ivar = 1:nvars
            if icue <3
                thisData = squeeze(nanmean(avgOnsetResp_pop(:,ivar,:,icue,itrig,:),6));
            elseif icue ==3
                thisData = squeeze(diff(nanmean(avgOnsetResp_pop(:,ivar,:,:,itrig,:),6),1,4));
            end
            subplot(nvars, 3, 3*(ivar-1)+icue);
            imagesc(winSamps, param.cardinalDir, thisData);
            set(gca, 'ytick',param.cardinalDir);
            if icue==1
                ylabel(psthNames{ivar});
            end
            mcolorbar(gca,.5);
            if ivar == 1
                if icue==1
                    title('wo cue');
                elseif icue ==2
                    title('w cue');
                elseif icue == 3
                    title('wcue-wocue');
                end
            end
        end
    end
    xlabel(['time from' onsetName ' [s]']);
    screen2png([onsetName '_all']);
    close all;
end

%% trace per saccade direction
itrig = 2;%tOnset
for idir=1:length(param.cardinalDir)
    nch = size(avgOnsetResp_pop,6);
    figure('position',[0 0 1400 1000]);
    nvars = size(avgOnsetResp_pop,2);
    tiledlayout('flow');
    axeshandle=[];
    for ivar = 1:nvars
        mData = squeeze(nanmean(avgOnsetResp_pop(idir,ivar,:,:,itrig,:),6));
        eData = 1/sqrt(nch)*squeeze(nanstd(avgOnsetResp_pop(idir,ivar,:,:,itrig,:),[],6));
        %subplot(nvars, 3, 3*(ivar-1)+icue);
        axeshandle(ivar) = nexttile;
        %imagesc(winSamps_onset, param.cardinalDir, thisData);
        errorbar(winSamps,mData(:,1),eData(:,1));
        hold on
        errorbar(winSamps,mData(:,2),eData(:,2));
        ylabel(psthNames{ivar});
        if ivar == 1
            legend('wo cue','w cue','location','northwest');
            title(num2str(param.cardinalDir(idir)));
        end
        vline(0);
        grid minor
    end
    % marginplots(axeshandle);%not compatible w tiledChartLayout
    xlabel(['time from' onsetName ' [s]']);
    screen2png([onsetName '_all_dir' num2str(param.cardinalDir(idir))]);
    close all;
end


%% limit channels by their preferred direction
itrig = 2;
onsetName = 'tOnset';
[~,prefDir] = max(mean(squeeze(avgOnsetResp_pop(:,1,tgtTimes,1,itrig,:)),2));
prefDir = squeeze(prefDir);

%nonanCh = find(~isnan(sum(sum(sum(sum(sum(avgOnsetResp_pop,1),2),3),4),5)));

for jdir = [1 5] %pref dir
    theseChannels = find(prefDir == jdir);
    for idir=[1 5] %stim dir
        nch = size(avgOnsetResp_pop,6);
        figure('position',[0 0 1400 1000]);
        nvars = size(avgOnsetResp_pop,2);
        tiledlayout('flow');
        axeshandle=[];
        for ivar = 1:nvars
            mData = squeeze(nanmean(avgOnsetResp_pop(idir,ivar,:,:,itrig,theseChannels),6));
            eData = 1/sqrt(nch)*squeeze(nanstd(avgOnsetResp_pop(idir,ivar,:,:,itrig,theseChannels),[],6));
            %subplot(nvars, 3, 3*(ivar-1)+icue);
            axeshandle(ivar) = nexttile;
            %imagesc(winSamps_onset, param.cardinalDir, thisData);
            errorbar(winSamps,mData(:,1),eData(:,1));
            hold on
            errorbar(winSamps,mData(:,2),eData(:,2));
            ylabel(psthNames{ivar});
            if ivar == 1
                legend('wo cue','w cue','location','northwest');
                title(['stimDir: ' num2str(param.cardinalDir(idir)) ...
                    ' prefDir: ' num2str(param.cardinalDir(jdir))]);
            end
            vline(0);
            grid minor
        end
        % marginplots(axeshandle);%not compatible w tiledChartLayout
        xlabel(['time from' onsetName ' [s]']);
        screen2png([onsetName '_all_StimDir' num2str(param.cardinalDir(idir)) ...
            '_PrefDir' num2str(param.cardinalDir(jdir))]);
        close all;
    end
end


%% resp aligned to preferred saccde direction
%preferred direction per channel
%[~, prefDirs] = max(squeeze(mean(mean(avgOnsetResp_pop(:,1,tgtTimes,:,2,:),3),4))); 
[~, prefDirs] = max(squeeze(mean(avgOnsetResp_pop(:,1,tgtTimes,1,2,:),3))); 
prefDirs = repmat(prefDirs, [size(avgOnsetResp_pop,2)*size(avgOnsetResp_pop,4)*size(avgOnsetResp_pop,5),1]);
prefDirs = prefDirs(:);

% target direction x varType x time x wo/w cue x trigger type x channels
pavgOnsetResp = permute(avgOnsetResp_pop, [3 1 2 4 5 6]);
[centeredDir, centeredData]  = alignMtxDir(pavgOnsetResp, tgtTimes, param.cardinalDir, prefDirs);
centeredData = permute(centeredData, [2 3 1 4 5 6]);



itrig = 2;
figure('position',[0 0 1400 1000]);
for icue = 1:3
    nvars = size(avgOnsetResp_pop,2);
    for ivar = 1:nvars
        if icue <3
            thisData = squeeze(nanmean(centeredData(:,ivar,:,icue,itrig,:),6));
        elseif icue ==3
            thisData = squeeze(diff(nanmean(centeredData(:,ivar,:,:,itrig,:),6),1,4));
        end
        subplot(nvars, 3, 3*(ivar-1)+icue);
        imagesc(winSamps, centeredDir, thisData);
        set(gca, 'ytick', centeredDir);
        if icue==1
            ylabel(psthNames{ivar});
        end
        mcolorbar(gca,.5);
        if ivar == 1
            if icue==1
                title('wo cue');
            elseif icue ==2
                title('w cue');
            elseif icue == 3
                title('wcue-wocue');
            end
        end
    end
end
xlabel(['time from' onsetName ' [s]']);
screen2png([onsetName '_all_centered_w']);
close all;

%% tuning curves 
ncolumns = length(times)-1;
itrig = 2;
figure('position',[0 0 1400 1000]);
for it = 1:length(times)-1
    tidx = intersect(find(winSamps>times(it)), find(winSamps<times(it+1)));

    nvars = size(avgOnsetResp_pop,2);
    for ivar = 1:nvars
        meanData = squeeze(nanmean(mean(centeredData(:,ivar,tidx,:,itrig,:),3),6));
        seData = 1/sqrt(size(centeredData,6))*squeeze(nanstd(mean(centeredData(:,ivar,tidx,:,itrig,:),3),[],6));
        
        subplot(nvars, ncolumns, ncolumns*(ivar-1)+it);
        errorbar(repmat(centeredDir',[1 2]), meanData, seData);
        if ivar==1
            title(['after tonset' num2str(times(it)) '-' ...
                num2str(times(it+1))]);
        end
        if it==1
            ylabel(psthNames{ivar});
        end
        %axis padded
    end
end
marginplots;
screen2png([onsetName '_all_tuningCurves']);
close all;


%% fit gaussian per channel
for it = 1:length(times)-1
    tidx = intersect(find(winSamps>times(it)), find(winSamps<times(it+1)));
    ravgOnsetResp_pop(:,:,it,:,:,:) = mean(avgOnsetResp_pop(:,:,tidx,:,:,:),3);
% target direction x varType x time x wo/w cue x trigger type x channels
end
[nDir, nVar, nT, nWwo, nTrig, nChannels] = size(ravgOnsetResp_pop);


addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\Matteobox');
pars = [];
err = [];
for it = 1:nT
    for ivar = 1:nVar
        for iWwo = 1:nWwo
            for iTrig = nTrig
                for ich = 1:nChannels
                    %cache = reshape(ravgOnsetResp_pop(:,:,it,:,:,ich), [nDir, nVar*nWwo*nTrig]);
                    [pars(:,ivar, it, iWwo, iTrig,ich), err(ivar, it, iWwo, iTrig,ich)] ...
                        = fitori(param.cardinalDir', ravgOnsetResp_pop(:,ivar, it, iWwo, iTrig,ich));
                end
            end
        end
    end
end

