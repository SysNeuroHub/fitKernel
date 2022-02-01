%% use fixation period?
%% extract alpha power from PSTH



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

marginSacc = 0.2;%s
cutoffFreq = 1;
durTh = 0.5;
sizeTh = 5;

[loadNames, months, dates, channels] = getMonthDateCh(animal, rootFolder);

% to obtain index of specified month&date&channel
thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
    regexptranslate('wildcard','09September\01\*_ch27.mat'))));

previousDate = [];
for idata = 1:length(channels) %1061;%865;%
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
            t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
            load(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');
        end
        
        t_cat = eyeData_rmotl_cat.t; %correct?
        
        periSaccadeTrace = event2Trace(t_cat, [catEvTimes.saccadeStartTimes ...
            catEvTimes.saccadeEndTimes], marginSacc);
        [periSaccadeStart, periSaccadeEnd] = trace2Event(periSaccadeTrace, t_cat);
        
        theseTrials = intersect(find(~isnan(catEvTimes.tOnset)), find(~isnan(catEvTimes.cOnset)));
        periStimTrace = event2Trace(t_cat, [catEvTimes.tOnset(theseTrials) catEvTimes.cOnset(theseTrials)], marginSacc);
        [periStimStart, periStimEnd] = trace2Event(periStimTrace, t_cat);
        
        
        otlTrace = event2Trace(t_cat,[catEvTimes.outlierStartTimes catEvTimes.outlierEndTimes]);
        blkTrace = event2Trace(t_cat, [catEvTimes.blinkStartTimes catEvTimes.blinkEndTimes]);
        saccTrace = event2Trace(t_cat, [periSaccadeStart periSaccadeEnd]);
        stimTrace = event2Trace(t_cat,[ periStimStart periStimEnd]);
        excludeTrace = (otlTrace+blkTrace+saccTrace+stimTrace>0);
        excludeTimes = trace2Event(excludeTrace, t_cat);
        
        %% detect pupil dilation/constriction
        [dlStartTimes, dlEndTimes, csStartTimes, csEndTimes] = ...
            detectPupilOnsets(eyeData_rmotl_cat, cutoffFreq, durTh, sizeTh, excludeTimes);

        %% combine all events
        pupilTimes = [dlStartTimes; csStartTimes];
        pupilTypes = [ones(length(dlStartTimes),1);
            2*ones(length(csStartTimes),1)];
        pupilLabels = {'dilation st','constriction st'};
        
        %% spike/prdiction response to pupil dilation/constriction
        
        [avgPupilResp, winSamps, singlePupilResp, sortedPupilLabels, uniquePupilLabels] ...
            = eventLockedAvg(cat(1,PSTH_f',predicted_all, predicted, predictorInfo.predictors_r(17,:)), ...
            predictorInfo.t_r, pupilTimes, pupilTypes, param.figTWin);
        
        sePupilResp = [];
        for ip = 1:length(pupilLabels)
            thesePupilEvs = find(pupilTypes == ip);
            sePupilResp(ip,:,:) = std(singlePupilResp(thesePupilEvs, :,:),[],1);
        end
        psthNames = cat(2,{'psth','predicted_all'},param.predictorNames ,{'pdiam'});

        nvars = size(avgPupilResp,2);
        figure('position',[0 0 400 1000]);
        ax=[];
        for ivar = 1:nvars
            ax(ivar)=subplot(nvars, 1, ivar);
            %imagesc(winSamps, param.cardinalDir, squeeze(avgDlResp(:,ivar,:)));
            plot(winSamps, squeeze(avgPupilResp(:,ivar,:)));
            
            %set(gca, 'ytick',1:4,'yticklabel',pupilLabels);
            %xlabel('time from pupil onset [s]');
            title(psthNames{ivar});
            if ivar==1
                title(['dilation: ' num2str(length(dlStartTimes)) ', constriction:' num2str(length(csStartTimes))])
            end
        end
        linksubaxes('y',ax(1:nvars-1));
        marginplots;
        %legend(pupilLabels,'location','best');
        mlegend(pupilLabels,gca,'northoutside')
        
        screen2png(fullfile(saveFolder,['pupilOn_' saveSuffix]));
        close all;
        
        %% save results
       save(saveName, 'pupilTimes','pupilLabels','pupilTypes',...
           'avgPupilResp', 'winSamps', 'singlePupilResp','-append');
            %'pspec_psth','pspec_parea','faxis_psth','faxis_parea');
        clear spk_all dd kernel kernel_x kernel_y psth_all mDir seDir mDir_pred seDir_pred
        
        previousDate = thisDate;
        
    end
end
