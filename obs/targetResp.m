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
    regexptranslate('wildcard','08August\13\*_ch19.mat'))));

previousDate = [];
for idata = thisdata;%1:length(channels) %1061;%865;%
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
        
        successTr = intersect(find(dd.successTrials), find(~isnan(catEvTimes.tOnset))); %2/2/22
           
        [~,minDirIdx] = arrayfun(@(x)(min(abs(circ_dist(pi/180*x, pi/180*param.cardinalDir)))), dd.targetloc);

        %% spike/prdiction response to targetOnsets
        [avgTgtResp, winTgtSamps, singleTgtResp, sortedTgtLabels, uniqueTgtLabels] ...
            = eventLockedAvg(cat(1,PSTH_f',predicted_all, predicted), ...
            predictorInfo.t_r, catEvTimes.tOnset(successTr), minDirIdx(successTr), param.figTWin);
        
        psthNames = cat(2,{'psth','predicted_all'}, param.predictorNames);
        tgtTimes = intersect(find(winTgtSamps>0.03), find(winTgtSamps<0.25));
        pavgTgtResp = permute(avgTgtResp, [3 1 2]);
        [centeredDir, cavgTgtResp]  = alignMtxDir(pavgTgtResp, tgtTimes, param.cardinalDir);
        
        nvars = size(avgTgtResp,2);
        figure('position',[0 0 400 1000]);
        ax=[];
        for ivar = 1:nvars

            %thisImage = squeeze(avgTgtResp(:,ivar,:));
            thisImage = squeeze(cavgTgtResp(:,:,ivar))';

            ax(ivar)=subplot(nvars, 1, ivar);
            %imagesc(winTgtSamps, param.cardinalDir, thisImage);
            imagesc(winTgtSamps, centeredDir, thisImage);
            
            %plot(winSamps, squeeze(avgTgtResp(:,ivar,:)));
            %             if ivar==1
            %                 original = squeeze(avgTgtResp(:,ivar,:));
            %                 crange = prctile(original(:),[5 99]);
            %             end
            
            %set(gca, 'ytick',param.cardinalDir,'yticklabel',param.cardinalDir);
            set(gca, 'ytick',centeredDir,'yticklabel',centeredDir);
            set(gca,'tickdir','out');
            %xlabel('time from pupil onset [s]');
            title(psthNames{ivar});
            %caxis(crange);
            caxis(prctile(thisImage(:),[5 99]));
            %if ivar==1
                mcolorbar(gca,.5);
            %end
            xlim([0 param.figTWin(2)]);
        end
        
        screen2png(fullfile(saveFolder,['tgtOn_' saveSuffix]));
        close all;
        
        
        %% direction tuning of each time point
        directions = pi/180*repmat(param.cardinalDir(minDirIdx(successTr))',...
            [1, size(singleTgtResp,2) size(singleTgtResp,3)]);
        %singleTgtResp(isnan(singleTgtResp))=0;
        
        singleTgtResp_m = singleTgtResp;
        singleTgtResp_m(:,3:end,:) = singleTgtResp_m(:,3:6,:) + mFiringRate;
        cvar = squeeze(circ_var(directions, singleTgtResp_m,1)); %[variables x time]
        cmean = squeeze(circ_mean(directions, singleTgtResp_m,1)); %[variables x time]
        
        %% temporal sequence of each kernels wrt to measured
        %option1:distance between preferred direction
        %dist_cmean = circ_dist(cmean(1,:),cmean(2,:));
        %
        %option2: similarity of tuning curve
        %for tt = 1:size(avgTgtResp,3)
        %    R=corrcoef(avgTgtResp(:,1,tt),avgTgtResp(:,4,tt));
        %    corr(tt) = R(1,2);
        %end
        
        subplot(211);
        imagesc(winTgtSamps, 1:nvars, squeeze(cmean));
        set(gca, 'ytick',1:nvars,'yticklabel',psthNames);
        title('circular mean');
        colormap(gca,'hsv');
        xlim([0 param.figTWin(2)]);
        mcolorbar(gca,.5);
        
        subplot(212);
        imagesc(winTgtSamps, 1:nvars, squeeze(cvar));
        set(gca, 'ytick',1:nvars,'yticklabel',psthNames);
        title('circular variance');
        mcolorbar(gca,.5);
        xlim([0 param.figTWin(2)]);
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
