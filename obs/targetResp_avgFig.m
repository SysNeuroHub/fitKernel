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

%[loadNames, months, dates, channels] = getMonthDateCh(animal, rootFolder);

% to obtain index of specified month&date&channel
thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
    regexptranslate('wildcard','12December\13\*_ch13.mat'))));

previousDate = [];
for idata = thisdata%1:length(channels) %1061;%865;%
    datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
    disp(datech);
    
    saveSuffix = [animal replace(datech,'/','_')];
    
    thisDate = [months{idata} '_' dates{idata}];
    
     saveName = fullfile(saveFolder, [saveSuffix '.mat']);
        
        
    if exist(saveName,'file')

        load(saveName, 'PSTH_f','predicted_all', 'predicted','kernelInfo',...
         't_r','mFiringRate','param','winTgtSamps', 'singleTgtResp');
       
     tgtResp = squeeze(mean(singleTgtResp,1));

        psthNames = cat(2,{'psth','predicted_all'}, param.predictorNames);
        
subplot(211);
plot(winTgtSamps, tgtResp(1,:), 'linewidth',2,'color','k');
vline(0);

subplot(212);
plot(winTgtSamps, tgtResp(3:end,:), 'linewidth',2);
hold on
plot(winTgtSamps, tgtResp(2,:), 'linewidth',2,'color','k');
linksubaxes('xy');
vline(0);
marginplots;

screen2png(['tOnset_avgFig_' saveSuffix]);
close 

        previousDate = thisDate;
        
    end
end
