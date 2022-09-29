
if ispc
    %saveFolder = '//storage.erc.monash.edu/shares/R-MNHS-Syncitium/Shared/Daisuke/cuesaccade_data';
    %rootFolder = '//storage.erc.monash.edu.au/shares/R-MNHS-Physio/SysNeuroData/Monash Data/Joanita/2021/cuesaccade_data/';
    saveFolder = 'E:/tmp/cuesaccade_data';
    saveFigFolder = [saveFolder, '/20220722'];
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

% parameters
load('E:\tmp\cuesaccade_data\param20220719','param');
psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);

previousDate = [];
for idata = 1:length(channels) %1061;%865;%
    datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
    disp(datech);
    
    saveSuffix = [animal replace(datech,'/','_')];
    
    thisDate = [months{idata} '_' dates{idata}];
    
    saveName = fullfile(saveFolder, [saveSuffix '.mat']);
    load(saveName,'mFiringRate');
    if mFiringRate < 5
        disp(['skipped as mFiringRate<5']);
        continue;
    end
    
    load(loadNames{idata});
    
    
        %% prepare behavioral data (common across channels per day)
        eyeName = fullfile(saveFolder,['eyeCat_' thisDate '.mat']);
        
        disp('loading eye/predictor data');
        load(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');
        load(fullfile(saveFolder,['eyeCat_' thisDate '.mat']));
        
        
        load(saveName, 'PSTH_f','kernelInfo','predicted_all','predicted');
        
        [~, pdiam] = getPupilDiameter(eyeData_rmotl_cat);
        pdiam_r = interp1(eyeData_rmotl_cat.t, pdiam, predictorInfo.t_r)';
        
        
        param.bandpassFreq = [7 13];
        [PSTH_ff, PSTH_amp, PSTH_phase] = getAmpPhase(PSTH_f, 1/param.dt_r, param.bandpassFreq);
        
        
        
        %% triggered by tOnsets
       
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
        
        
        avgOnsetResp_m = [];
        for it = 1:length(param.times)-1
            tidx = intersect(find(winSamps_onset>param.times(it)), find(winSamps_onset<param.times(it+1)));
            avgOnsetResp_m(:,:,it,:) = mean(avgOnsetResp(:,:,tidx,:),3);
            % target direction x varType x time x wo/w cue x trigger type x channels
        end
        
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
