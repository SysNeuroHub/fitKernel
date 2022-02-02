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

       load(loadNames{idata}, 'dd');
        
        %% prepare behavioral data (common across channels per day)
        eyeName = fullfile(saveFolder,['eyeCat_' thisDate '.mat']);
        if  ~strcmp(thisDate, previousDate)             
            load(eyeName);
            t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
            load(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');
        end
        
        t_cat = eyeData_rmotl_cat.t; %correct?
        
        
        %% target onsets
        dirs = unique(dd.targetloc);
        
        successTr = find(dd.successTrials);

        
        %% combine all events
        %         pupilTimes = [dlStartTimes; csStartTimes];
        %         pupilTypes = [ones(length(dlStartTimes),1);
        %             2*ones(length(csStartTimes),1)];
        %pupilLabels = {'dilation st','constriction st'};
        
        
        [~,minDirIdx] = arrayfun(@(x)(min(abs(circ_dist(pi/180*x, pi/180*param.cardinalDir)))), dd.targetloc);

        %% spike/prdiction response to targetOnsets
        [avgTgtResp, winTgtSamps, singleTgtResp, sortedTgtLabels, uniqueTgtLabels] ...
            = eventLockedAvg(cat(1,PSTH_f',predicted_all, predicted), ...
            predictorInfo.t_r, catEvTimes.tOnset(successTr), minDirIdx(successTr), param.figTWin);
        
        psthNames = cat(2,{'psth','predicted_all'}, param.predictorNames);
        
        nvars = size(avgTgtResp,2);
        figure('position',[0 0 700 1000]);
        ax=[];
        for ivar = 1:nvars
            ax(ivar)=subplot(nvars, 1, ivar);
            imagesc(winTgtSamps, param.cardinalDir, squeeze(avgTgtResp(:,ivar,:)));
            %plot(winSamps, squeeze(avgTgtResp(:,ivar,:)));
            if ivar==1
                original = squeeze(avgTgtResp(:,ivar,:));
                crange = prctile(original(:),[5 99]);
            end
            
            set(gca, 'ytick',param.cardinalDir,'yticklabel',param.cardinalDir);
            set(gca,'tickdir','out');
            %xlabel('time from pupil onset [s]');
            title(psthNames{ivar});
            caxis(crange);
            if ivar==1
                mcolorbar(gca,.5);
            end
        end
        
        screen2png(fullfile(saveFolder,['tgtOn_' saveSuffix]));
        close all;
        
        
        %% direction tuning of each time point
        directions = pi/180*repmat(param.cardinalDir(minDirIdx(successTr))',...
            [1, size(singleTgtResp,2) size(singleTgtResp,3)]);
        %singleTgtResp(isnan(singleTgtResp))=0;
        cvar = circ_var(directions, singleTgtResp,1);
        cmean = circ_mean(directions, singleTgtResp,1);
        
        subplot(211);
        imagesc(winTgtSamps, 1:nvars, squeeze(cmean));
        set(gca, 'ytick',1:nvars,'yticklabel',psthNames);
        title('circular mean')
        mcolorbar(gca,.5);
        
        subplot(212);
        imagesc(winTgtSamps, 1:nvars, squeeze(cvar));
        set(gca, 'ytick',1:nvars,'yticklabel',psthNames);
        title('circular variance');
        mcolorbar(gca,.5);
        xlabel('Time from target onset [s]');
        
        screen2png(fullfile(saveFolder,['cmeanvar_' saveSuffix]));
        close all;
      
        %% save results
       save(saveName,...
           'avgTgtResp', 'winTgtSamps', 'singleTgtResp','cvar','cmean','-append');
            %'pspec_psth','pspec_parea','faxis_psth','faxis_parea');
        
        
        previousDate = thisDate;
        
    end
end
