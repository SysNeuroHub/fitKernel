%% compute 

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
[loadNames, months, dates, channels] = getMonthDateCh(animal, rootFolder);

% to obtain index of specified month&date&channel
thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
    regexptranslate('wildcard','09September\21\*_ch27.mat'))));

ii=1;
previousDate = [];
for idata = thisdata%1:length(channels) %1061;%865;%
    datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
    disp(datech);
    
    saveSuffix = [animal replace(datech,'/','_')];
    
    thisDate = [months{idata} '_' dates{idata}];
    
    saveName = fullfile(saveFolder, [saveSuffix '.mat']);
    
    
    if exist(saveName,'file')
        
        load(saveName, 'PSTH_f','predicted_all', 'predicted','kernelInfo',...
            't_r','mFiringRate','param','winTgtSamps', 'singleTgtResp',...
            'avgOnsetResp','winSamps_conset');
        
        %tgtResp = squeeze(mean(singleTgtResp,1));% originally used for my presentation
        tgtResp = squeeze(mean(avgOnsetResp(:,:,:,1)));
        winTgtSamps = winSamps_conset;
        
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
        %
        %         screen2png(['tOnset_avgFig_' saveSuffix]);
        %         close
        
        
        observed = tgtResp(1,winTgtSamps>0);
        allMdl = tgtResp(3,winTgtSamps>0);
        visionMdl = tgtResp(4,winTgtSamps>0);
        eyeMdl = tgtResp(5,winTgtSamps>0);
        
        mobserved = mean(tgtResp(1,winTgtSamps<=0&winTgtSamps>-0.1));
        mallMdl = mean(tgtResp(3,winTgtSamps<=0&winTgtSamps>-0.1));
        mvisionMdl = mean(tgtResp(4,winTgtSamps<=0&winTgtSamps>-0.1));
        meyeMdl = mean(tgtResp(5,winTgtSamps<=0&winTgtSamps>-0.1));

        %version1: 
        %         visResp(ii) = max(visionMdl) - mvisionMdl;
        %         eyeResp(ii) = max(eyeMdl) - meyeMdl;
        %         obsResp(ii) = max(observed) - mobserved;
        
        %version2: take into account of negative resp. eg. 746
        visResp(ii) = max(abs(visionMdl - mvisionMdl));
        eyeResp(ii) = max(abs(eyeMdl - meyeMdl));
        obsResp(ii) = max(abs(observed - mobserved));

        %version3: without abs
        
        datech_pop{ii} = datech;
        
        ii = ii+1;
        previousDate = thisDate;
        
    end
end

neyeResp = (eyeResp./obsResp);
nvisResp =  visResp./obsResp;

subplot(211); plot(eyeResp, visResp, '.');
xlabel('eye velocity amp');
ylabel('vision amp');
squareplot;
axis padded;

subplot(212); plot(neyeResp, nvisResp, '.');
squareplot;
axis padded;
xlabel('eye velocity amp/observed amp');
ylabel('vision amp/observed amp');
title(['corr: ' num2str(corr(neyeResp', nvisResp'))]);

figure;
histogram(obsResp);
xlabel('observed max response [Hz]');
ylabel('#units');
