

%% recorded data
animal = 'hugo';
rootFolder = '\\storage.erc.monash.edu.au\shares\R-MNHS-Physio\SysNeuroData\Monash Data\Joanita\2021/cuesaccade_data/';
saveFolder = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\Daisuke\cuesaccade_data';

excDur = 0.25; %[s]
respWin = [0 0.2]; %[s]

[loadNames, months, dates, channels] = getMonthDateCh(animal, rootFolder);

thisdata = find(1-cellfun(@isempty, regexp(loadNames, regexptranslate('wildcard','11November\26\*_ch28.mat'))));

for idata = thisdata%1:length(channels) %1061;%865;%

datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
disp(datech);

saveSuffix = [animal replace(datech,'/','_')];

thisDate = [months{idata} '_' dates{idata}];

    
    saveName = fullfile(saveFolder, [saveSuffix '.mat']);
    
    if ~exist(saveName,'file')
        continue;
    end
    
    load(saveName, 'mFiringRate',...
        'param','eyeData_rmotl_cat');
    
    %% raw data
    load(loadNames{idata}, 'ephysdata','dd');
    
    nTrials = length(dd.eye);
    fs_eye = median([dd.eye.fs]);
    eyeData = dd.eye;
    
    %% concatenate across trials
    spk_all = ephysdata.spikes.spk;
    [spk_all_cat, t_cat] = concatenate_spk(spk_all, {dd.eye.t});
    clear spk_all
    
    disp('loading eye/predictor data');
    load(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');
    load(fullfile(saveFolder,['eyeCat_' thisDate '.mat']), ...
        'onsets_cat','meta_cat','blinks','outliers','t_tr');
    
    
    %% get saccade onset times for each direction
    cardinalDir = param.cardinalDir;
    startSacc = meta_cat.STARTSACC;
    endSacc = meta_cat.ENDSACC;
    if length(startSacc) ~= length(endSacc)%41
        nSaccs = min(length(startSacc), length(endSacc));  
        ngIdx=find(endSacc(1:nSaccs)-startSacc(1:nSaccs)<0);
        
        if isempty(ngIdx)
            startSacc = startSacc(1:nSaccs);
            endSacc = endSacc(1:nSaccs);
        %else
        % FILL ME
        end
    end
    
    assert(isempty(find(endSacc-startSacc<0)));
    okSacc = find(endSacc-startSacc>0);
    startSacc = startSacc(okSacc);
    endSacc = endSacc(okSacc);
    
    x = eyeData_rmotl_cat.x;
    y = eyeData_rmotl_cat.y;
    t = eyeData_rmotl_cat.t;
    eyeRad = atan2(y, x); %[-pi pi]
    dist = sqrt(y.^2+x.^2);
    
    %sanity check
    plot(t, dist);
    hold on
    plot(startSacc,0,'ro');
    plot(endSacc,0,'gx');
    
    minDirIdx = zeros(length(t),1);
    for tt = 1:length(t)
        [~,minDirIdx(tt)] = min(abs(circ_dist(eyeRad(tt), pi/180*cardinalDir)));
    end
    
    dirIndex = zeros(length(startSacc),1);
    for isacc = 1:length(startSacc)
        tsnippet = intersect(find(t>startSacc(isacc)), find(t<endSacc(isacc)));
        dirIndex(isacc) = mode(minDirIdx(tsnippet));
    end
    
    
    %% times to omit from the analysis
    tOnset = onsets_cat.tOnset;
    cOnset = onsets_cat.cOnset;
    validEvents = intersect(find(~isnan(tOnset)), find(~isnan(cOnset)));
    tOnset_nonan = tOnset(validEvents);
    cOnset_nonan = cOnset(validEvents);
    excEventT_cat = zeros(length(t_cat),1);
    for itr = 1:length(validEvents)
        evStartT = min(tOnset_nonan(itr)-excDur, cOnset_nonan(itr)-excDur);
        excStartT = max(t_cat(1), evStartT);
        evEndT = max(tOnset_nonan(itr)+excDur, cOnset_nonan(itr)+excDur);
        excEndT = min(t_cat(end), evEndT);
        theseTimes = intersect(find(t_cat>=excStartT), find(t_cat<=excEndT));
        
        excEventT_cat(theseTimes) = 1;
    end
    
     
    %% omit saccades near target/cue onsets
    startSaccTidx = [];
    for isacc = 1:length(startSacc)
        [~,startSaccTidx(isacc)] = min(abs(t_cat - startSacc(isacc)));
    end
    [~,startSaccIdx] = setdiff(startSaccTidx, find(excEventT_cat));
    startSaccT = t_cat(startSaccTidx(startSaccIdx));
    saccDir = param.cardinalDir(dirIndex(startSaccIdx));
    
    %% obtain concatenated psth and its predictions
    % PSTH_r = getPSTH(spk_all_cat, predictorInfo.t_r);
    % dt_r = median(diff(predictorInfo.t_r));
    % PSTH_f = filtPSTH(PSTH_r, dt_r, param.psth_sigma, 2);%causal
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
    
    %% recorded and predicted psth aligned to saccade onsets
    [avgSaccResp, winSamps_sacc, singleSaccResp, sortedLabels, uniqueLabels] ...
        = eventLockedAvg(cat(1,PSTH_f',predicted_all, predicted), predictorInfo.t_r, startSaccT, saccDir, param.figTWin);
    
    
    %% behavior signals aligned to saccade onsets
    [avgBhvResp, ~, singleBhvResp] ...
        = eventLockedAvg(predictorInfo.predictors_r, predictorInfo.t_r, startSaccT, saccDir, param.figTWin);
    
    
    %% individual response to saccades
    %predictorInfo
    idir = 3
   theseEvents = find(sortedLabels == param.cardinalDir(idir));
    imagesc(squeeze(singleBhvResp(theseEvents, 8+idir,:))');
    
    %% avg across saccades
    
    %  plot(winSamps_sacc,squeeze(avgBhvResp(7,:,:))');
    
    
    
    
    varName = [{'recorded'},{'predictedAll'}, param.predictorNames(:)'];
    nvars = size(avgSaccResp,2);
    figure('position',[0 0 400 1000]);
    for ivar = 1:nvars
        subplot(nvars, 1, ivar);
        imagesc(winSamps_sacc, param.cardinalDir, squeeze(avgSaccResp(:,ivar,:)));
        set(gca, 'ytick',param.cardinalDir);
        xlabel('time from saccade onset [s]');
        ylabel(varName{ivar});
        mcolorbar(gca,.5);
        
    end
    
    screen2png(['saccOn_' saveSuffix]);
    close;
    
    
    
    %% direction selectivity
    theseTimes = intersect(find(winSamps_sacc >= respWin(1)), find(winSamps_sacc <= respWin(2)));
    mResp = squeeze(mean(avgSaccResp(:,:,theseTimes),3));
    mDir_sacc = circ_mean(pi/180*param.cardinalDir', mResp(:,1));
    mDir_sacc_pred = circ_mean(pi/180*param.cardinalDir', mResp(:,4));
    
    figure;
    subplot(211);
    bins = [param.cardinalDir 360]-22.5;%hack
    histogram(saccDir, bins,'orientation','horizontal');
    ylabel('saccade direction [deg]');
    xlabel('#saccades');
    axis ij;
    
    subplot(212);
    plot(mResp(:,1),param.cardinalDir, 'b');hold on;
    plot(mResp(:,4),param.cardinalDir ,'g');
    legend(varName{1},varName{4},'location','northeastoutside');
    title(['mDir ' num2str(mDir_sacc) '(recorded), ' num2str(mDir_sacc_pred) '(pred), time ' ...
        num2str(respWin(1)) ' - ' num2str(respWin(2)) 's']);
    ylabel('saccade direction [deg]');
    xlabel('firing rate [/s]');
    axis ij;
    screen2png(['dirTun_sacc_' saveSuffix]);
    close
    
    
    save(saveName, 'mDir_sacc', 'mDir_sacc_pred', 'avgSaccResp', 'winSamps_sacc','-append');
end


