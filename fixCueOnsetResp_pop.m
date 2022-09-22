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
    regexptranslate('wildcard','09September\01\*_ch25.mat'))));

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
avgfOnsetResp_pop = [];
navgfOnsetResp_pop = [];
avgCueResp_pop = [];
navgCueResp_pop = [];

kk = [];
for idata = 1:length(channels)
    %     %rerun: 68
    %     %03March/23; 03March/25; 03March/26; 03March/29; 03March/30; 04April/08; 04April/09; 06June/14; 06June/18; %07July/28
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
    %
    %
    %     load(loadNames{idata});
    %
    %     if dataType == 0
    %         spk_all = ephysdata.spikes.spk;
    %         chName = ['_ch' num2str(ephysdata.spikes.chanIds)];
    %         clear ephysdata
    %     end
    %
    %
    %     if ~isempty(spk_all)
    %
    %
    %         nTrials = length(dd.eye);
    %         fs_eye = median([dd.eye.fs]);
    %         eyeData = dd.eye;
    %
    %         %% concatenate across trials
    %         [spk_all_cat, t_cat] = concatenate_spk(spk_all, {dd.eye.t});
    %         clear spk_all
    %         mFiringRate = length(spk_all_cat)/(t_cat(end)-t_cat(1)); %spks/s
    %         clear spk_all_cat t_cat
    saveName = fullfile(saveFolder, [saveSuffix '.mat']);
    if ~exist(saveName,'file')
        kk = [kk idata];
        continue;
    end
    load(saveName,'mFiringRate');
    if  ~exist('mFiringRate','var') || mFiringRate < 5
        kk = [kk idata];
        disp([chName 'skipped as mFiringRate<5']);
        continue;
    end
    
    
    load(saveName, 'avgfOnsetResp', 'avgCueResp','winSamps_fc');
    
    if ~exist('avgfOnsetResp','var') || sum(isnan(avgfOnsetResp(:)))>0
        %cue condition did not exist
        continue;
    end
    
    navgfOnsetResp = avgfOnsetResp/max(avgfOnsetResp(:));
    navgCueResp = avgCueResp/max(avgCueResp(:));
    
    avgfOnsetResp_pop = cat(4,avgfOnsetResp_pop, avgfOnsetResp);
    navgfOnsetResp_pop = cat(4,navgfOnsetResp_pop, navgfOnsetResp);
    
    avgCueResp_pop = cat(4,avgCueResp_pop, avgCueResp);
    navgCueResp_pop = cat(4,navgCueResp_pop, navgCueResp);
    
    clear avgfOnsetResp avgCueResp 
    
end

mavgCueResp = squeeze(mean(avgCueResp_pop,4));
mnavgCueResp = squeeze(mean(navgCueResp_pop,4));

mavgfOnsetResp = squeeze(mean(avgfOnsetResp_pop,4));
mnavgfOnsetResp = squeeze(mean(navgfOnsetResp_pop,4));

crange = prctile(mavgfOnsetResp(:),[.1 99.9]);
ax(1)=subplot(131);
imagesc(winSamps_fc, 1:6, squeeze(mavgfOnsetResp(1,:,:)));
set(gca,'yticklabel',psthNames);
caxis(crange);
mcolorbar;
ylabel('without cue');

ax(2)=subplot(132);
imagesc(winSamps_fc, 1:6, squeeze(mavgfOnsetResp(2,:,:)));
caxis(crange);
mcolorbar;
xlabel('time from fixation onset');
ylabel('with cue');

ax(3)=subplot(133);
imagesc(winSamps_fc, 1:6, mavgCueResp);
mcolorbar;
xlabel('time from cue onset');

figure('position',[0 0 500 1000]);
for ivar = 1:6
    subplot(6,1,ivar);
    plot(winSamps_fc, squeeze(mavgfOnsetResp(:,ivar,:)));
    axis padded;
    ylabel(psthNames{ivar});
end
legend('wo cue','w cue');
xlabel('time from fixation onset');