
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
    regexptranslate('wildcard','12December/09/*_ch1.mat'))))

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
param.predictorNames = {'vision','eyespeed','pdiam','blink'}; %eyeposition
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

previousDate = [];
for idata = 916:length(channels) %1061;%865;%
    %886 NG in getSaccDir
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
    
    
    load(loadNames{idata});
    
    if dataType == 0
        spk_all = ephysdata.spikes.spk;
        chName = ['_ch' num2str(ephysdata.spikes.chanIds)];
        clear ephysdata
    end
    
    
    if ~isempty(spk_all)
        
        
        nTrials = length(dd.eye);
        fs_eye = median([dd.eye.fs]);
        eyeData = dd.eye;
        
        %% concatenate across trials
        [spk_all_cat, t_cat] = concatenate_spk(spk_all, {dd.eye.t});
        clear spk_all
        mFiringRate = length(spk_all_cat)/(t_cat(end)-t_cat(1)); %spks/s
        if mFiringRate < 5
            disp([chName 'skipped as mFiringRate<5']);
            continue;
        end
        
        
        %% prepare behavioral data (common across channels per day)
        eyeName = fullfile(saveFolder,['eyeCat_' thisDate '.mat']);
        if  ~strcmp(thisDate, previousDate)
            load(eyeName, 'eyeData_rmotl_cat',...
                'onsets_cat','meta_cat','blinks','outliers','pspec_parea','faxis_parea','t_tr');
            
            t_r = (eyeData_rmotl_cat.t(1):param.dt_r:eyeData_rmotl_cat.t(end))';
            
            %check startsacc and endsacc (to be done in concatenate_eye)
            if length(meta_cat.STARTSACC) ~= length(meta_cat.ENDSACC)%41
                nSaccs = min(length(meta_cat.STARTSACC), length(meta_cat.ENDSACC));
                ngIdx=find(meta_cat.ENDSACC(1:nSaccs)-meta_cat.STARTSACC(1:nSaccs)<0);
                
                if isempty(ngIdx)
                    meta_cat.STARTSACC = meta_cat.STARTSACC(1:nSaccs);
                    meta_cat.ENDSACC = meta_cat.ENDSACC(1:nSaccs);
                    %else
                    % FILL ME?
                end
            end
            
            assert(isempty(find(meta_cat.ENDSACC-meta_cat.STARTSACC<0)));
            okSacc = find(meta_cat.ENDSACC-meta_cat.STARTSACC>0);
            meta_cat.STARTSACC = meta_cat.STARTSACC(okSacc);
            meta_cat.ENDSACC = meta_cat.ENDSACC(okSacc);
            %             tic
            %             [eyeData_cat, onsets_cat, meta_cat] = concatenate_eye(eyeData, dd);%takes 90s
            %             t0=toc
            %             t_tr={eyeData.t};
            %
            %
            %             %% detect and interpolate blinks
            %             %removing outliers helps for larger kernels
            %             tic
            %             [eyeData_rmblk_cat, blinks] = removeBlinksEDF(eyeData_cat, meta_cat, param.marginSize,1);
            %             close;
            %             t1=toc
            %
            %             %% remove parea outliers based on diff
            %             tic
            %             [eyeData_rmotl_cat, outliers] = removePareaOutliers(eyeData_rmblk_cat, ...
            %                 param.marginSize, param.pareaTh, param.pareaDiffTh);
            %             screen2png(['rmOtl_' saveSuffix]);
            %             close;
            %             t2=toc
            
            %% detect saccades
            [startSacc, endSacc] = selectSaccades(meta_cat.STARTSACC, meta_cat.ENDSACC,...
                t_cat, (blinks+outliers>0));    % exculde blink and outlier periods
            
            catEvTimes = onsets_cat;
            [catEvTimes.blinkStartTimes, catEvTimes.blinkEndTimes] = trace2Event(blinks, t_cat);
            [catEvTimes.outlierStartTimes, catEvTimes.outlierEndTimes] = trace2Event(outliers, t_cat);
            catEvTimes.saccadeStartTimes = startSacc;
            catEvTimes.saccadeEndTimes = endSacc;
            %catEvents.t
            %catTraces.saccadeDirs < no need to be here
            save(eyeName, 'catEvTimes','-append');
            clear  meta_cat
            
            %             [pspec_parea,faxis_parea] = pmtm(eyeData_rmotl_cat.parea, 10, ...
            %                 length(eyeData_rmotl_cat.parea), fs_eye);%slow
            
            
            %% prepare predictor variables
            predictorInfo = preparePredictors(dd, eyeData_rmotl_cat, t_r, param, catEvTimes);
            save(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');
        else
            disp('loading eye/predictor data');
            load(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');
            load(fullfile(saveFolder,['eyeCat_' thisDate '.mat']));
        end
        
        
        
        
        %% obtain kernels!
        disp('fit kernels')
        %param.lagRange = repmat([-0.5 0.5],[18,1]);
        %param.lagRange(1:8,1)=0;
        [predicted_all, PSTH_f, kernelInfo] = fitPSTH(spk_all_cat, ...
            predictorInfo.t_r, predictorInfo.predictors_r, param.psth_sigma, param.lagRange, param.ridgeParams);

        predicted = zeros(predictorInfo.nPredictors, length(predictorInfo.t_r));
        for ivar = 1:predictorInfo.nPredictors
            if ivar==1
                theseVarIdx = 1:predictorInfo.npredVars(1);
            else
                theseVarIdx = sum(predictorInfo.npredVars(1:ivar-1))+1:sum(predictorInfo.npredVars(1:ivar));
            end
            
            predicted(ivar, :) = predictXs(predictorInfo.t_r, predictorInfo.predictors_r(theseVarIdx,:), ...
                kernelInfo.intercept, kernelInfo.kernel(:,theseVarIdx), param.lagRange);
        end
        
        figure('position',[0 0 1000 500]);
        subplot(1,2,1);
        plot(predictorInfo.t_r, PSTH_f, 'color',[.5 .5 .5]);hold on
        plot(predictorInfo.t_r, predicted_all, 'linewidth',2);
        %xlim([1935 1966]);
        xlim([100 200])
        legend('recorded','fitted');
        xlabel('time [s]'); ylabel('firing rate [Hz]');
        
        title(['expval: ' num2str(kernelInfo.expval), ', R: ' num2str(kernelInfo.corrcoef)]);
        
        a2=subplot(3,2,2);
        thisIm = kernelInfo.kernel(:,1:predictorInfo.npredVars(1))';
        crange = prctile(abs(thisIm(:)),99);
        imagesc(kernelInfo.tlags, param.cardinalDir, thisIm);
        caxis([-crange crange]);
        set(gca,'ytick',param.cardinalDir);
        xlabel('time from targetOnset [s]');
        mcolorbar(a2,.5);
        
        a3=subplot(3,2,4);
        thisIm = kernelInfo.kernel(:,predictorInfo.npredVars(1)+1:sum(predictorInfo.npredVars(1:2)))';
        crange = prctile(abs(thisIm(:)),99);
        imagesc(kernelInfo.tlags,param.cardinalDir, thisIm);
        caxis([-crange crange]);
        set(gca,'ytick',param.cardinalDir);
        xlabel('time from eye movement [s]');
        mcolorbar(a3,.5);
        
        a4=subplot(3,2,6);
        plot(kernelInfo.tlags, kernelInfo.kernel(:,sum(predictorInfo.npredVars(1:2))+1:end)');
        xlabel('time from pupil dilation/blink [s]');
        axis tight;
        
        %set(gca,'ytick',predictorInfo.predictors_r);
        screen2png(['kernels_' saveSuffix]);
        close;
        
        
        %% spectral analysis
        %         disp('spectral analysis');
        %         figure('position',[0 0 1000 500]);
        %         subplot(121);
        %         loglog(faxis_parea, pspec_parea);
        %         grid on
        %         axis tight
        %         xlabel('frequency [Hz]'); ylabel('psd');
        %         title('parea');
        %
        %         [pspec_psth,faxis_psth] = pmtm(PSTH_f, 10, length(PSTH_f), 1/param.dt_r);%slow
        %         subplot(122);
        %         semilogy(faxis_psth, pspec_psth);
        %         grid on
        %         axis tight
        %         xlabel('frequency [Hz]'); ylabel('psd');
        %         title('psth');
        %         screen2png(['pspec_parea_psth' saveSuffix]);
        %         close;
        
        
        %% decompose back into individual trials
        eyeData_rmotl_tr = decompose_eye(eyeData_rmotl_cat, t_tr);
        
        
        %% upsampling
        psth_predicted_all = interp1(predictorInfo.t_r, predicted_all, eyeData_rmotl_cat.t, 'linear');
        psth_predicted = zeros(predictorInfo.nPredictors, length(eyeData_rmotl_cat.t));
        for ivar = 1:predictorInfo.nPredictors
            psth_predicted(ivar,:) = interp1(predictorInfo.t_r, predicted(ivar,:), eyeData_rmotl_cat.t, 'linear');
        end
        PSTH_f_upsample = interp1(predictorInfo.t_r, PSTH_f, eyeData_rmotl_cat.t, 'linear');
        % psth_predicted =  predicted_slow+predicted_fast+predicted_fast_x+predicted_fast_y;
        
        
        %% triggered by external events TOBE REMOVED
        psth_predicted_all_tr = cell(nTrials, 1);
        psth_predicted_tr = cell(nTrials, predictorInfo.nPredictors);
        psth_denoised_tr = cell(nTrials,1);
        psth_tr = cell(nTrials,1);
        headidx = 1;
        for itr = 1:nTrials
            nFrames = length(t_tr{itr});
            
            psth_predicted_all_tr{itr} = psth_predicted_all(headidx:headidx+nFrames-1);
            psth_tr{itr} = PSTH_f_upsample(headidx:headidx+nFrames-1); %11/1/22
            psth_denoised_tr{itr} = psth_tr{itr} - psth_predicted_all_tr{itr};
            
            for ivar = 1:predictorInfo.nPredictors
                psth_predicted_tr{itr,ivar} = psth_predicted(ivar, headidx:headidx+nFrames-1);
            end
            headidx = headidx+nFrames;
        end
        
        psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);
        psth_all = cat(2,psth_tr, psth_predicted_all_tr, psth_predicted_tr); %cell(#trials, 2+nPredictors)
        
        %% avg trials
        %         [f, psth_snippet, parea_snippet, dist_snippet, taxis_snippet] ...
        %             = pupilFigureAvgSingle(dd, eyeData_rmotl_tr, psth_all, param.evName, param.figTWin);
        [f, psth_snippet, pdiam_snippet, dist_snippet, taxis_snippet] ...
            = pupilFigure(dd, eyeData_rmotl_tr, psth_all, param.evName, param.figTWin);
        legend(psthNames(2:end),'location','northwest');
        screen2png(['pupilPsth_' param.evName saveSuffix], f);
        close;
        
        
        %% response to saccades
        tOnset = catEvTimes.tOnset;
        cOnset = catEvTimes.cOnset;
        validEvents = intersect(find(~isnan(tOnset)), find(~isnan(cOnset)));
        tOnset = tOnset(validEvents);
        cOnset = cOnset(validEvents);
        
        tcOnset_trace = event2Trace(t_cat, [tOnset; cOnset], 2*0.5);
        excEventT_cat = (tcOnset_trace + blinks + outliers > 0); %28/1/22
        
        [startSaccNoTask, endSaccNoTask] = selectSaccades(catEvTimes.saccadeStartTimes, ...
            catEvTimes.saccadeEndTimes, t_cat, excEventT_cat);%param.minSaccInterval);
        
        [saccDirNoTask, dirIndexNoTask] = getSaccDir(startSaccNoTask, endSaccNoTask, ...
            eyeData_rmotl_cat, param.cardinalDir);
        
        [avgSaccResp, winSamps_sacc, singleSaccResp, sortedSaccLabels, uniqueSaccLabels] ...
            = eventLockedAvg(cat(1,PSTH_f',predicted_all, predicted), ...
            predictorInfo.t_r, startSaccNoTask, saccDirNoTask, param.figTWin);
        
        nvars = size(avgSaccResp,2);
        figure('position',[0 0 400 1000]);
        for ivar = 1:nvars
            subplot(nvars, 1, ivar);
            imagesc(winSamps_sacc, param.cardinalDir, squeeze(avgSaccResp(:,ivar,:)));
            set(gca, 'ytick',param.cardinalDir);
            xlabel('time from saccade onset [s]');
            ylabel(psthNames{ivar});
            mcolorbar(gca,.5);
            
        end
        
        screen2png(['saccOn_' saveSuffix]);
        close;
        
        
        %% TODO response to pupil dilation / constiction
        
        %% individual trials
        theseTimes = intersect(find(taxis_snippet > param.respWin(1)), find(taxis_snippet < param.respWin(2)));
        msnippet = squeeze(mean(psth_snippet(theseTimes, :, :), 1));
        Rmsnippet = corrcoef(msnippet);
        
        fsnippet = figure('position',[0 0 1900 1400]);
        ax4 = [];
        for itype = 1:length(psthNames)
            ax4(itype)=subplot(2,length(psthNames),itype);
            thisData = squeeze(psth_snippet(:,:,itype))';
            imagesc(taxis_snippet, 1:size(psth_snippet,2), thisData);
            %ylim([200 300]);
            title(psthNames{itype});
            if itype==1
                ylabel('trial');
                xlabel(['time from ' param.evName]);
            end
            crange = prctile(thisData(:),[1 99]);
            if diff(crange)==0
                crange = [crange(1) crange(1)+1];
            end
            caxis(crange);
            mcolorbar(ax4(itype), .2);
            
            ax4(itype+length(psthNames))=subplot(2,length(psthNames),itype+length(psthNames));
            plot(msnippet(:,itype), msnippet(:,1),  '.');
            axis square;
            title(['vs ' psthNames{itype} ', R: ' num2str(Rmsnippet(1,itype))]);
            xlabel(psthNames{itype})
        end
        screen2png(['indtrials_' saveSuffix]);
        close;
        
        
        %% compare direction tuning (trig by dd.tOnset) before/after denoising
        disp('direction tuning');
        dirs = unique(dd.targetloc);
        
        %         theseTr = find(dd.successTrials);
        %         [mDir, winSamps_sacc, singleDirResp, sortedDirLabels, uniqueDirLabels] ...
        %             = eventLockedAvg(cat(1,PSTH_f',predicted_all), ...
        %             predictorInfo.t_r, catEvTimes.tOnset(theseTr), pi/180*dd.targetloc(theseTr), param.figTWin);
        %         cvar = circ_var(sortedDirLabels, singleDirResp);
        %         cmean = circ_mean(sortedDirLabels, singleDirResp);
        
        [mDir, sdDir, cmean, cvar, nTrials_dtune] = dirTuning(psth_tr, eyeData_rmotl_tr, dd, param.respWin);
        [mDir_pred, sdDir_pred, cmean_pred, cvar_pred ] = dirTuning(psth_predicted_all_tr, eyeData_rmotl_tr, dd, param.respWin);
        seDir = sdDir./sqrt(nTrials_dtune);
        seDir_pred = sdDir_pred./sqrt(nTrials_dtune);
        
        figure('position',[0 0 1000 400]);
        errorbar(dirs/180*pi, mDir, seDir, 'b');hold on
        errorbar(dirs/180*pi, mDir_pred, seDir_pred, 'g');
        legend('observed','fitted','location','southeast');
        tname = sprintf('direction tuning %d-%d[ms]\nobserved cvar: %.2f, predicted cvar: %.2f',...
            1e3*param.respWin(1), 1e3*param.respWin(2),cvar, cvar_pred);
        title(tname);
        
        screen2png(['dirTuning_' saveSuffix]);
        close all
        
        
        %% save results
        saveName = fullfile(saveFolder, [saveSuffix '.mat']);
        
        save(saveName, 'PSTH_f','predicted_all', 'predicted','kernelInfo','psth_all','mDir',...
            'seDir','mDir_pred','seDir_pred','cvar','cmean', 't_r',...
            'cvar_pred','cmean_pred','mFiringRate','param',...
            'winSamps_sacc', 'singleSaccResp', 'sortedSaccLabels');
            %'pspec_psth','pspec_parea','faxis_psth','faxis_parea');
        clear spk_all dd kernel kernel_x kernel_y psth_all mDir seDir mDir_pred seDir_pred
        
        previousDate = thisDate;
        
    end
end