
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
    regexptranslate('wildcard','04April\08\*_ch25.mat'))));

%% omit data
% no saccade response
% low spontaneous firing
% low number of successful trials

% parameters
param.th_frate = 5;
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
param.times = 1e-3*(0:30:390);%[0 50 100 150 200 250 300 350 400];

previousDate = [];
for idata = 154%:length(channels) %1061;%865;%
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
        
%         %% concatenate across trials
%         [spk_all_cat, t_cat] = concatenate_spk(spk_all, {dd.eye.t});
%         clear spk_all
%         mFiringRate = length(spk_all_cat)/(t_cat(end)-t_cat(1)); %spks/s
%         if mFiringRate < param.th_frate
%             disp([chName 'skipped as mFiringRate<' num2str(param.th_frate)]);
%             continue;
%         end
        
        
        %% prepare behavioral data (common across channels per day)
        eyeName = fullfile(saveFolder,['eyeCat_' thisDate '.mat']);
        
        disp('loading eye/predictor data');
        load(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');
        load(fullfile(saveFolder,['eyeCat_' thisDate '.mat']));
        
        
        saveName = fullfile(saveFolder, [saveSuffix '.mat']);
        load(saveName, 'PSTH_f','kernelInfo','predicted_all','predicted');
        
        [~, pdiam] = getPupilDiameter(eyeData_rmotl_cat);
        pdiam_r = interp1(eyeData_rmotl_cat.t, pdiam, predictorInfo.t_r)';
        
        
        param.bandpassFreq = [7 13];
        [PSTH_ff, PSTH_amp, PSTH_phase] = getAmpPhase(PSTH_f, 1/param.dt_r, param.bandpassFreq);
        
        
        
        %% triggered by tOnsets
        psthNames = cat(2,{'psth','psth_f','predicted_all'},...
            param.predictorNames,'pdiam ori');
        
       
        onset = catEvTimes.(param.evName);
        
        avgOnsetResp = nan(length(param.cardinalDir),length(psthNames),51,2);
        fitPars = []; fitErr = []; fittedOnsetResp_m = [];
        singleOnsetResp_m = [];
        avgOnsetResp_m = [];
        for icue = 1:2
            validEvents = intersect(find(~isnan(onset)), find(dd.cueOn==icue-1));
            %< include successful/unsucessful saccades
            
            onsetTimes = onset(validEvents);
            tgtDir = getTgtDir(dd.targetloc(validEvents), param.cardinalDir);
            
            [~,dirIdx]=intersect(param.cardinalDir, unique(tgtDir));
            [avgOnsetResp(dirIdx,:,:,icue), winSamps_onset, singleOnsetResp, ...
                sortedOnsetLabels{icue}, uniqueOnsetLabels] ...
                = eventLockedAvg(cat(1,PSTH_f',PSTH_ff', predicted_all, predicted, pdiam_r), ...
                predictorInfo.t_r, onsetTimes, tgtDir, param.figTWin);
            
            %% resample in time
            for it = 1:length(param.times)-1
                tidx = intersect(find(winSamps_onset>param.times(it)), ...
                    find(winSamps_onset<param.times(it+1)));
                singleOnsetResp_m{icue}(:,:,it) = mean(singleOnsetResp(:,:,tidx),3);
                avgOnsetResp_m(:,:,it,icue) = mean(avgOnsetResp(:,:,tidx,icue),3);
                % target direction x varType x time x wo/w cue x trigger type x channels
            end
        end
        
        minAvgOnsetResp_m = squeeze(min(mean(avgOnsetResp_m,4),[],1));
        % varType x time
        
        %% fit curve
        for icue = 1:2
            for it = 1:length(param.times)-1
                for ivar = 1:length(psthNames)
                    
                    if it==1
                        initPars = [];
                    else
                        initPars = fitPars(:,ivar, it-1, icue);
                    end
                    %TODO: replace w fitResponse
                    [fitPars(:,ivar, it, icue), fitErr(ivar, it, icue)] ...
                        = fitoriWrapped(sortedOnsetLabels{icue}, singleOnsetResp_m{icue}(:,ivar,it),...
                        [], [nan nan 0 minAvgOnsetResp_m(ivar,it) nan],'',20, initPars);
                    %[Dp(peak angle), Rp(peak amp), Rn(2nd peak amp = 0), Ro(baseline), sigma(width)]
                    fittedOnsetResp_m(:,ivar,it,icue) ...
                        = orituneWrapped(fitPars(:,ivar,it,icue), param.cardinalDir);
                end
            end
        end
        
        avgOnsetResp_m = [];
        for it = 1:length(param.times)-1
            tidx = intersect(find(winSamps_onset>param.times(it)), find(winSamps_onset<param.times(it+1)));
            avgOnsetResp_m(:,:,it,:) = mean(avgOnsetResp(:,:,tidx,:),3);
            % target direction x varType x time x wo/w cue x trigger type x channels
        end
        
        %% inspect fitting result
        figure('position',[0 0 1400 1000]);
        for ifit = 1:2
            for icue = 1:2
                nvars = size(avgOnsetResp_m,2);
                for ivar = 1:nvars
                    if (ifit==1)
                        thisData = squeeze(avgOnsetResp_m(:,ivar,:,icue));
                        if icue == 1
                            irow=1;
                        elseif icue == 2
                            irow=3;
                        end
                    elseif (ifit==2)
                        thisData = squeeze(fittedOnsetResp_m(:,ivar,:,icue));
                        if icue == 1
                            irow=2;
                        elseif icue == 2
                            irow=4;
                        end
                    end
                    
                    
                    subplot(nvars, 4, 4*(ivar-1) + irow);
                    imagesc(param.times(1:end-1)+.5*mean(diff(param.times)), param.cardinalDir, thisData);
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
                        end
                    end
                end
            end
        end
        xlabel(['time from' param.evName ' [s]']);
        screen2png([param.evName '_' saveSuffix '_fit']);
        close all;
        
        
        %% w vs wo cue
        colormap('parula');
        figure('position',[0 0 1400 1000]);
        for icue = 1:3
            nvars = size(avgOnsetResp,2);
            for ivar = 1:nvars
                if icue <3
                    thisData = squeeze(avgOnsetResp(:,ivar,:,icue));
                elseif icue ==3
                    thisData = squeeze(diff(avgOnsetResp(:,ivar,:,:),1,4));
                end
                subplot(nvars, 3, 3*(ivar-1)+icue);
                imagesc(winSamps_onset, param.cardinalDir, thisData);
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
        xlabel(['time from' param.evName ' [s]']);
        screen2png([param.evName '_' saveSuffix]);
        close all;
        
        
        %% save results
        
        save(saveName, 'avgOnsetResp','avgOnsetResp_m','fitPars','fitErr',...
            'fittedOnsetResp_m','PSTH_ff','PSTH_amp','PSTH_phase','-append');
        clear spk_all dd kernel kernel_x kernel_y psth_all mDir seDir mDir_pred seDir_pred
        
        previousDate = thisDate;
        
    end
end
