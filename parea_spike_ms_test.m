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

cutoffFreqs_ms = [1e-2 1e-1 1e0];
bpFreq = [5 15];

[loadNames, months, dates, channels] = getMonthDateCh(animal, rootFolder);

% to obtain index of specified month&date&channel
thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
    regexptranslate('wildcard','12December\16\*_ch5.mat'))));

previousDate = [];
for idata = thisdata%14:length(channels) %1061;%865;%
    datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
    disp(datech);
    
    saveSuffix = [animal replace(datech,'/','_')];
    
    thisDate = [months{idata} '_' dates{idata}];
    
    saveName = fullfile(saveFolder, [saveSuffix '.mat']);
    
    
    if exist(saveName,'file')
        
        load(saveName, 'PSTH_f','predicted_all', 'predicted','kernelInfo',...
            't_r','mFiringRate','param');
        
               
        
        %% prepare behavioral data (common across channels per day)
        eyeName = fullfile(saveFolder,['eyeCat_' thisDate '.mat']);
        if  ~strcmp(thisDate, previousDate)
            load(eyeName);
           % t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
            load(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');
        end
        
        %pdiam_r = predictorInfo.predictors_r(17,:);%NG!!!
        %pdiam_r = predictorInfo.predictors_r(17,:).^2';%somehow this works
        parea_r = interp1(eyeData_rmotl_cat.t,eyeData_rmotl_cat.parea, t_r);
        cutoffFreqs_ms = [1e-3 5e-3 1e-2 1e-1];
       [corrCoef_parea_spk] = parea_spike_ms(PSTH_f, parea_r, t_r, cutoffFreqs_ms);
%        screen2png(fullfile(saveFolder, ['_parea_spk_' saveSuffix '.mat']));
       screen2png(fullfile(saveFolder, ['parea_spk_' saveSuffix '.png']));
        
       
        order = 3;
        fs = 1/median(diff(t_r));
        Wn = bpFreq/(fs/2);
        [b,a]=butter(order, Wn, 'bandpass');
        signal = PSTH_f;
        ntotFrames = length(signal);
        signal_c = filtfilt(b,a,cat(1,flipud(signal), ...
            signal, flipud(signal)));
        alphaOsc = signal_c(ntotFrames+1:2*ntotFrames);
        analytic = hilbert(alphaOsc);
        alphaAmp =abs(analytic);
       
       [corrCoef_parea_alpha] = parea_spike_ms(alphaAmp, parea_r, t_r, cutoffFreqs_ms);
       screen2png(fullfile(saveFolder, ['parea_alpha_' saveSuffix '.png']));
       
       close all;
       
        %% save results
        save(saveName, 'corrCoef_parea_spk','corrCoef_parea_alpha','cutoffFreqs_ms','-append');
        
        previousDate = thisDate;
        
    end
end
