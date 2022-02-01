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

paramPupilOn.marginSacc = 0.2;%s
paramPupilOn.cutoffFreq = 2;%1;
paramPupilOn.durTh = 0.25;%0.5;
paramPupilOn.sizeTh = 5;
bpFreq = [5 15];%[hz] %alpha oscillation of PSTH

[loadNames, months, dates, channels] = getMonthDateCh(animal, rootFolder);

% to obtain index of specified month&date&channel
% thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
%     regexptranslate('wildcard','09September\01\*_ch27.mat'))));

previousDate = [];
for idata = 865%1:length(channels) %1061;%865;%
    datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
    disp(datech);
    
    saveSuffix = [animal replace(datech,'/','_')];
    
    thisDate = [months{idata} '_' dates{idata}];
    
    saveName = fullfile(saveFolder, [saveSuffix '.mat']);
    
    
    if exist(saveName,'file')
        
        load(saveName, 'PSTH_f','predicted_all', 'predicted','kernelInfo',...
            't_r','mFiringRate','param');
        
        %% extract alpha amplitude in PSTH
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
        
        %sanity check
        subplot(211);
        plot(t_r, PSTH_f, t_r, alphaOsc);
        subplot(212);
        plot(t_r, alphaAmp);
        linksubaxes('x');
        
        
        %% prepare behavioral data (common across channels per day)
        eyeName = fullfile(saveFolder,['eyeCat_' thisDate '.mat']);
        if  ~strcmp(thisDate, previousDate)
            load(eyeName);
            t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
            load(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');
        end
        
        t_cat = eyeData_rmotl_cat.t; %correct?
        
        periSaccadeTrace = event2Trace(t_cat, [catEvTimes.saccadeStartTimes ...
            catEvTimes.saccadeEndTimes], paramPupilOn.marginSacc);
        [periSaccadeStart, periSaccadeEnd] = trace2Event(periSaccadeTrace, t_cat);
        
        theseTrials = intersect(find(~isnan(catEvTimes.tOnset)), find(~isnan(catEvTimes.cOnset)));
        periStimTrace = event2Trace(t_cat, [catEvTimes.tOnset(theseTrials) catEvTimes.cOnset(theseTrials)], paramPupilOn.marginSacc);
        periStimTimes = trace2Event(periStimTrace, t_cat);
        
        
        otlTrace = event2Trace(t_cat,[catEvTimes.outlierStartTimes catEvTimes.outlierEndTimes]);
        blkTrace = event2Trace(t_cat, [catEvTimes.blinkStartTimes catEvTimes.blinkEndTimes]);
        saccTrace = event2Trace(t_cat, [periSaccadeStart periSaccadeEnd]);
        fixSuccess = intersect(find(~isnan(catEvTimes.fOnset)), find(~isnan(catEvTimes.tOnset)));
        fixTrace = event2Trace(t_cat, [catEvTimes.fOnset(fixSuccess) catEvTimes.tOnset(fixSuccess)]);
        %stimTrace = event2Trace(t_cat,[ periStimStart periStimEnd]);
        excludeTrace = (otlTrace+blkTrace+saccTrace+(1-fixTrace)>0);
        excludeTimes = trace2Event(excludeTrace, t_cat);
        
        %% detect pupil dilation/constriction
        [dlStartTimes, dlEndTimes, csStartTimes, csEndTimes] = ...
            detectPupilOnsets(eyeData_rmotl_cat, paramPupilOn.cutoffFreq, ...
            paramPupilOn.durTh, paramPupilOn.sizeTh, excludeTimes);
        
        
        %sanity check
        %         dlTimes = event2Trace(t_cat, [dlStartTimes dlEndTimes]);
        %         csTimes = event2Trace(t_cat, [csStartTimes csEndTimes]);
        %         subplot(311);
        %         plot(t_cat, fixTrace, 'g',t_cat, saccTrace,'c',t_cat,blkTrace,'k',t_cat, otlTrace,'y');
        %         grid on;
        %         subplot(312);
        %         plot(t_cat, excludeTrace);grid on;
        %         subplot(313);
        %         plot(t_cat, dlTimes,'r',t_cat, csTimes,'b','linewidth',2);
        %         grid on
        %         linksubaxes('x');
        %         marginplots;
        
        %% combine all events
        if isempty(dlStartTimes) || isempty(csStartTimes)
            disp('no pupil dilation/constriction detected');
            continue;
        end
        pupilTimes = [dlStartTimes; csStartTimes];
        pupilTypes = [ones(length(dlStartTimes),1);
            2*ones(length(csStartTimes),1)];
        pupilLabels = {'dilation st','constriction st'};
        
        %% spike/prdiction response to pupil dilation/constriction
        [avgPupilResp, winSamps, singlePupilResp, sortedPupilLabels, uniquePupilLabels] ...
            = eventLockedAvg(cat(1,PSTH_f',alphaOsc', alphaAmp', predicted_all, predicted, predictorInfo.predictors_r(17,:)), ...
            predictorInfo.t_r, pupilTimes, pupilTypes, param.figTWin);
        
        sePupilResp = [];
        for ip = 1:length(pupilLabels)
            thesePupilEvs = find(pupilTypes == ip);
            sePupilResp(ip,:,:) = 1/sqrt(length(thesePupilEvs))*std(singlePupilResp(thesePupilEvs, :,:),[],1);
        end
        psthNames = cat(2,{'psth','alpha','alphaAmp','predicted_all'},param.predictorNames ,{'pdiam'});
        
        nvars = size(avgPupilResp,2);
        figure('position',[0 0 600 1000]);
        ax=[];
        for ivar = 1:nvars
            ax(ivar)=subplot(nvars, 1, ivar);
            %imagesc(winSamps, param.cardinalDir, squeeze(avgDlResp(:,ivar,:)));
            errorbar(repmat(winSamps',[1 length(pupilLabels)]), ...
                squeeze(avgPupilResp(:,ivar,:))',squeeze(sePupilResp(:,ivar,:))');
            
            %set(gca, 'ytick',1:4,'yticklabel',pupilLabels);
            %xlabel('time from pupil onset [s]');
            title(psthNames{ivar});
            if ivar==1
                title(['dilation: ' num2str(length(dlStartTimes)) ', constriction:' num2str(length(csStartTimes))])
            end
        end
        %linksubaxes('y',ax(1:nvars-1));
        marginplots;
        %legend(pupilLabels,'location','best');
        mlegend(pupilLabels,gca,'northoutside')
        
        screen2png(fullfile(saveFolder,['pupilOn_f_' saveSuffix]));
        close all;
        
        %% save results
        save(saveName, 'pupilTimes','pupilLabels','pupilTypes',...
            'avgPupilResp', 'sePupilResp','winSamps', 'singlePupilResp',...
            'alphaAmp','alphaOsc','bpFreq','paramPupilOn''-append');
        %'pspec_psth','pspec_parea','faxis_psth','faxis_parea');
        clear spk_all dd kernel kernel_x kernel_y psth_all mDir seDir mDir_pred seDir_pred
        
        previousDate = thisDate;
        
    end
end
