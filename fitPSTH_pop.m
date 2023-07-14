%% get ready
%addpath(genpath('C:\Users\dshi0006\git'))
%setenv('COMPUTERNAME', 'MU00011697');
[saveServer, rootFolder] = getReady();




%% recorded data
animal = 'hugo';
dataType = 0;%0: each channel, 1: all channels per day

load(fullfile(saveServer,'param20230405.mat'),'param');
param.mfiringRateTh = 5;
param.expvalTh = 3;%0;
param.ntargetTrTh = 200;
param.ptonsetRespTh = 0.05;

%alldata = [];
mFiringRate_pop = [];
kernel_pop = [];
expval_pop = [];
corrcoef_pop = [];
corrcoef_pred_spk_pop = [];
PtonsetResp_pop = [];
PsaccResp_pop = [];
expval_ind_pop = [];
expval_tgt_pop = [];
corr_avgtgt_pop = [];
expval_avgtgt_pop = [];
ntargetTrials_pop = [];
ntotTrials_pop = [];
id_pop = [];
Rsqadj_pop = [];
for yy = 3
    switch yy
        case 1
            year = '2021';
        case 2
            year = '2022';
        case 3
            year = '2023';
    end
    [loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);
    
    % to obtain index of specified month&date&channel
    %find(1-cellfun(@isempty, regexp(loadNames, regexptranslate('wildcard','12December\09\*_ch32.mat'))))
    
    %% omit data
    % no saccade response
    % low spontaneous firing
    % low number of successful trials
    
    nData = length(channels);
    %dataByYear = [];%struct(length(channels),1);
    %datech_pop = nan(length(channels),1);
    %corrcoef_pred_spk_pop = nan(length(channels),1);
    for idata = 1:length(channels)
        datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
        thisid = [animal '/' year '/' datech];
        disp(thisid);
        
        saveSuffix = [animal replace(datech,'/','_')];
        
        thisDate = [months{idata} '_' dates{idata}];
        
        saveFolder = fullfile(saveServer, year,animal);%17/6/23
        saveName = fullfile(saveFolder, [saveSuffix '.mat']);
        
        if exist(saveName, 'file')
            
            %datech_pop{idata} = datech;
            S = load(saveName, 'PSTH_f','predicted_all', 'predicted', ...
                'kernelInfo','t_r','cellclassInfo','param','mFiringRate','t_cat',...
                'dds');
            
            if isfield(S,'kernelInfo')
                mFiringRate_pop = cat(2,mFiringRate_pop,S.mFiringRate);
                kernel_pop = cat(3,kernel_pop,S.kernelInfo.kernel);
                expval_pop = cat(2,expval_pop,S.kernelInfo.expval);
                corrcoef_pop = cat(2,corrcoef_pop,S.kernelInfo.corrcoef);
                
                R = corrcoef(S.PSTH_f, S.predicted_all);
                corrcoef_pred_spk_pop = cat(2, corrcoef_pred_spk_pop, R(1,2));
                %id_pop = cat(2,id_pop, thisid);
                id_pop{numel(id_pop)+1} = thisid;
                
                eyeName = fullfile(saveFolder,[
                    'eyeCat_' animal thisDate '.mat']);
                load(fullfile(saveFolder,['predictorInfo_' animal thisDate '.mat']), ...
                    'predictorInfo');
                load(eyeName,'eyeData_rmotl_cat','catEvTimes');
                
                %                 if ~isfield(S,'dds')
                %                     load(loadNames{idata},'dd'); %slow to read
                %                     dds = getCueSaccadeSmall(dd);
                %                     save(saveName, 'dds','-append');
                %                    S.dds = dd;
                %                 end
                dd = S.dds;
                
                
                %% Rsq adj of subjset of variables
                nsub=3;
                Rsqadj = zeros(nsub,1);
                for jj = 1:nsub
                    switch jj
                        case 1 %full model
                            tgtGroups = 1:5;
                        case 2 %omit eye speed
                            tgtGroups = setxor(1:5, 2);
                        case 3 %omit eye position
                            tgtGroups = setxor(1:5, 3);
                    end
                    
                    [Rsqadjusted,rr,r0] = fitSubset(S.PSTH_f, predictorInfo, ...
                        tgtGroups, param);
                    
                    Rsqadj(jj) = Rsqadjusted;
                end
                Rsqadj_pop = [Rsqadj_pop Rsqadj];
                
                %% response to target
                PtonsetResp_pop = [PtonsetResp_pop S.cellclassInfo.PtonsetResp];
                
                %% response to saccade
                PsaccResp_pop = [PsaccResp_pop S.cellclassInfo.PsaccResp];
                
                %% explained variance per kernel
                expval = zeros(size(S.predicted,2)+1,1);
                expval(1,1) = getExpVal(S.PSTH_f-mean(S.PSTH_f), ...
                    S.predicted_all-mean(S.predicted_all));
                for ivar = 1:size(S.predicted,2)
                    expval(ivar+1,1) = getExpVal(S.PSTH_f-mean(S.PSTH_f), ...
                        S.predicted(:,ivar)-mean(S.predicted(:,ivar)));
                end
                expval_ind_pop = [expval_ind_pop expval];
                
                %% explained variance for target response
                expval_tgt(1,1) = getExpVal_tgt(S.PSTH_f, S.predicted_all, catEvTimes, S.t_r, [0 0.5]);
                expval_tgt(2:6,1) = getExpVal_tgt(S.PSTH_f, S.predicted, catEvTimes, S.t_r, [0 0.5]);
               
                expval_tgt_pop = [expval_tgt_pop expval_tgt];
                
                %% explained variance for target response averaged across trials
                [expval_avgtgt(1,1), corr_avgtgt(1,1)] = getExpVal_avgtgt(S.PSTH_f, S.predicted_all, ...
                    catEvTimes, S.t_r, [0 0.5], param.cardinalDir, dd);                
                [expval_avgtgt(2:6,1), corr_avgtgt(2:6,1)] = getExpVal_avgtgt(S.PSTH_f, S.predicted, ...
                    catEvTimes, S.t_r, [0 0.5], param.cardinalDir, dd);
                expval_avgtgt_pop = [expval_avgtgt_pop expval_avgtgt];
                corr_avgtgt_pop = [corr_avgtgt_pop corr_avgtgt];

                %% number of target trials
                ntargetTrials_pop = [ntargetTrials_pop sum(~isnan(catEvTimes.tOnset))];
                
                %% number of total trials
                ntotTrials_pop = [ntotTrials_pop numel(catEvTimes.tOnset)];
                
                
                %retrieve just once
                tlags = S.kernelInfo.tlags;
                clear S
            end
        end
    end
    %dataByYear = dataByYear(~isnan(dataByYear));
    %alldata = [alldata dataByYear(:)];
end

% save('fitPSTH_pop20220202','avgPupilResp_pop', '-append');
kernel_pop = squeeze(kernel_pop);
save(fullfile(saveServer,['fitPSTH_pop20230713' animal]),'mFiringRate_pop','kernel_pop','expval_pop','corrcoef_pop',...
    'corrcoef_pred_spk_pop','id_pop','ntotTrials_pop','ntargetTrials_pop','param',...
    'PsaccResp_pop','PtonsetResp_pop','expval_ind_pop','expval_tgt_pop','tlags',...
    'corr_avgtgt_pop','expval_avgtgt_pop','Rsqadj_pop');

%% apply inclusion critetia
% [okunits, mfiringRateOK, expvalOK, ntargetTrOK, ptonsetRespOK] ...
%     = inclusionCriteria(mFiringRate_pop, expval_pop, ntargetTrials_pop, PtonsetResp_pop, param);
% [okunits, mfiringRateOK, expvalOK, ntargetTrOK, ptonsetRespOK] ...
%     = inclusionCriteria(mFiringRate_pop, expval_tgt_pop(1,:), ntargetTrials_pop, PtonsetResp_pop, param);
[okunits, mfiringRateOK, expvalOK, ntargetTrOK, ptonsetRespOK] ...
    = inclusionCriteria(mFiringRate_pop, expval_ind_pop(1,:), ntargetTrials_pop, PtonsetResp_pop, param);

kernel_pop = kernel_pop(:,okunits);
expval_ind_pop = expval_ind_pop(:,okunits);
expval_tgt_pop = expval_tgt_pop(:,okunits);
expval_avgtgt_pop = expval_avgtgt_pop(:,okunits);
corr_avgtgt_pop = corr_avgtgt_pop(:,okunits);
id_pop = id_pop(okunits);
mFiringRate_pop = mFiringRate_pop(okunits);
Rsqadj_pop = Rsqadj_pop(:,okunits);


%% selected units
% theseIDs = {'hugo/2021/09September/01/25',...
%     'hugo/2021/11November/16/6',...
%     'hugo/2021/12December/14/13'};
% theseIDs = {'hugo/2021/09September/01/25',...
%     'hugo/2021/11November/02/18',...
%     'hugo/2022/07July/29/19'};
theseIDs = {'hugo/2022/03March/10/20',... %eye speed driven
    'hugo/2022/07July/29/19'}; %eye position driven

[~, selectedIDs] = intersect(id_pop, theseIDs);

%% histogram of individual explained variance
figure('position',[ 1120         454         787         500]);
for ii = 1:size(expval_ind_pop,1)
    ax(ii) = subplot(size(expval_ind_pop,1),2,2*ii-1);
    histogram(expval_ind_pop(ii,:),[-10:1:20]);
    if ii==1
        ylabel('all');
        title('expval ind pop');
    else
        ylabel(param.predictorNames{ii-1});
    end
    
    ax(ii) = subplot(size(expval_ind_pop,1),2,2*ii);
    histogram(expval_tgt_pop(ii,:),[-10:1:30]);
    if ii==1
        title('expval tgt pop');
    end
end
xlabel('Explained variance [%]')
savePaperFigure(gcf,['expval_' animal]);

%% scatter plot of 
fig = showScatterTriplets(Rsqadj_pop, ...
    {'full mdl','wo eye speed','wo eye pos'}, [0 .5], selectedIDs);
squareplots;
screen2png(['Rsqadj_' animal]);

%% scatter plot of individual explained variances
fig = showScatterTriplets(expval_ind_pop(2:4,:), ...
    param.predictorNames, [-4 20], selectedIDs);
screen2png(['expval_all_a' animal]);

fig = showScatterTriplets(100*expval_ind_pop(2:4,:)./expval_ind_pop(1,:), ...
    param.predictorNames, [-100 200], selectedIDs);
screen2png(['expval_all_r' animal]);

%% explained variance between kernels
fig = showScatterTriplets(expval_tgt_pop(2:4,:), ...
    param.predictorNames, []);%, selectedIDs);
screen2png(['expval_tgt_a_' animal]);close;

fig = showScatterTriplets(100*expval_tgt_pop(2:4,:)./expval_tgt_pop(1,:), ...
    param.predictorNames, [-100 200]);%, selectedIDs);
screen2png(['expval_tgt_r_' animal]);close;

%% correlation on avg response between kernels
% fig = showScatterTriplets(corr_tgt_pop(2:4,:), ...
%     param.predictorNames, [], selectedIDs);
% screen2png(['corr_tgt_a_' animal]);close;
% 
% fig = showScatterTriplets(100*corr_tgt_pop(2:4,:)./corr_tgt_pop(1,:), ...
%     param.predictorNames, [-20 100], selectedIDs);
% screen2png(['corr_tgt_r_' animal]);close;

%% explained variance on avg response between kernels
fig = showScatterTriplets(expval_avgtgt_pop(2:4,:), ...
    param.predictorNames, [], selectedIDs);
screen2png(['expval_avgtgt_a_' animal]);close;

fig = showScatterTriplets(100*expval_avgtgt_pop(2:4,:)./expval_avgtgt_pop(1,:), ...
    param.predictorNames, [-200 100], selectedIDs);
screen2png(['expval_avgtgt_r_' animal]);close;

%% correlation on avg response between kernels
fig = showScatterTriplets(corr_avgtgt_pop(2:4,:), ...
    param.predictorNames, [], selectedIDs);
screen2png(['corr_avgtgt_a_' animal]);close;

fig = showScatterTriplets(100*corr_avgtgt_pop(2:4,:)./corr_avgtgt_pop(1,:), ...
    param.predictorNames, [-20 100], selectedIDs);
screen2png(['corr_avgtgt_r_' animal]);close;


%% show average kernel before centering
[f, kernel_avg] = showKernel3(kernel_pop, tlags, param.cardinalDir, 0);
savePaperFigure(f,['avgKernel_' animal]);


%% centerring by preferred direction
tgtRange = [0.05 0.15; 0.03 0.25; -0.1 0.1];
[f, kernel_centered_avg] = showKernel3(kernel_pop, tlags, param.cardinalDir, 1, tgtRange);
savePaperFigure(f,['avgKernel_centered_' animal]);


%% preferred direction across 3 kernels
%plot3(prefdir{1},prefdir{2},prefdir{3},'.');
subplot(131);
%sigUnits = (prefdirPval{1}<0.05) & (prefdirPval{2}<0.05); %EMPTY
plot(prefdir{1}, prefdir{2},'.');
%plot(prefdir{1}(sigUnits), prefdir{2}(sigUnits),'b.');
[rho, pval] = circ_corrcc(prefdir{1}*pi/180,prefdir{2}*pi/180);
title(['rho:' num2str(rho) ', pval:' num2str(pval)])
xlabel('vision'); ylabel('eye speed');
axis equal square;
subplot(132);
plot(prefdir{2},prefdir{3},'.');
[rho, pval] = circ_corrcc(prefdir{2}*pi/180,prefdir{3}*pi/180);
title(['rho:' num2str(rho) ', pval:' num2str(pval)])
xlabel('eye speed'); ylabel('eye position');
axis equal square;
subplot(133);
plot(prefdir{3},prefdir{1},'.');
[rho, pval] = circ_corrcc(prefdir{3}*pi/180,prefdir{1}*pi/180);
title(['rho:' num2str(rho) ', pval:' num2str(pval)])
xlabel('eye position'); ylabel('vision');
axis equal square;

screen2png(['prefDirCorr_' animal]);
% circular correlation


%% pupil dilation/constriction
pupilLabels = {'dilation st','constriction st'};
psthNames_pupil = cat(2,{'psth','alpha','alphaAmp','predicted_all'},param.predictorNames ,{'pdiam'});
nChannels = size(avgPupilResp_pop,4);
mPupilResp = squeeze(mean(avgPupilResp_pop,4));
sePupilResp = 1/sqrt(nChannels)*squeeze(std(avgPupilResp_pop,[],4));
nvars = size(avgPupilResp_pop,2);
figure('position',[0 0 1000 1000]);
ax=[];
for ivar = 1:nvars
    ax(ivar)=subplot(5, 2, ivar);
    %imagesc(winSamps, param.cardinalDir, squeeze(avgDlResp(:,ivar,:)));
    errorbar(repmat(winSamps_sacc',[1 length(pupilLabels)]), ...
        squeeze(mPupilResp(:,ivar,:))',squeeze(sePupilResp(:,ivar,:))');
    
    %set(gca, 'ytick',1:4,'yticklabel',pupilLabels);
    %xlabel('time from pupil onset [s]');
    title(psthNames_pupil{ivar});
    xlim([-0.3 0.3]);
    vline(0);
end
%linksubaxes('y',ax(1:nvars-1));
legend(pupilLabels);
marginplots;
screen2png('pupilOn_pop');

%% saccade resp. recorded vs predicted
psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);
avgSacc = permute(squeeze(nanmean(avgSaccResp_pop,4)),[3 1 2]);
crange = prctile(avgSacc(:),[1 99]);
for ipred = 1:size(avgSacc,3)
    subplot(7,1,ipred);
    imagesc(winSamps_sacc, param.cardinalDir, squeeze(avgSacc(:,:,ipred))');
    ylabel(psthNames{ipred})
    %caxis(crange);
    vline(0);
    mcolorbar(gca,.5);
end
xlabel('time from saccade onset [s]');
screen2png('Sacc_pop');
savePaperFigure(gcf,'Sacc_pop');


%% saccade resp aligned to preferred saccde direction. recorded vs predicted
tgtTimes = intersect(find(winSamps>0.03), find(winSamps<0.25));
psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);
pavgSaccResp = permute(avgSaccResp_pop, [3 1 4 2]);
[centeredDir, centeredData_sacc]  = alignMtxDir(pavgSaccResp, tgtTimes, param.cardinalDir);
%[time x direction x channels x predictors]
avgCentSacc = squeeze(nanmean(centeredData_sacc,3));
crange = prctile(avgCentSacc(:),[1 99]);
for ipred = 1:size(avgCentSacc,3)
    subplot(6,1,ipred);
    imagesc(winSamps, centeredDir, squeeze(avgCentSacc(:,:,ipred))');
    ylabel(psthNames{ipred})
    caxis(crange);
    vline(0);
    mcolorbar(gca,.5);
end
xlabel('time from saccade onset [s]');
screen2png('centeredSacc_pop');
savePaperFigure(gcf,'centeredSacc_pop');

% 
% %% target resp aligned to preferred target direction
% pavgTgtResp = permute(avgTgtResp_pop, [3 1 4 2]);
% [centeredDir, centeredData_tgt]  = alignMtxDir(pavgTgtResp, tgtTimes, param.cardinalDir);
% %[time x direction x channels x predictors]
% avgCentTgt = squeeze(nanmean(centeredData_tgt,3));
% crange = prctile(avgCentTgt(:),[1 99]);
% for ipred = 1:size(avgCentTgt,3)
%     subplot(6,1,ipred);
%     imagesc(winSamps, centeredDir, squeeze(avgCentTgt(:,:,ipred))');
%     ylabel(psthNames{ipred})
%     caxis(crange);
%     vline(0);
%     mcolorbar(gca,.5);
% end
% xlabel('time from target onset [s]');
% screen2png('centeredTgt_pop');
% savePaperFigure(gcf,'centeredTgt_pop');
% 
% 
% %% kernel aligned to preferred direction
% %centerBin = 4;
% %centeredDir = 180/pi*circ_dist(pi/180*param.cardinalDir, pi/180*param.cardinalDir(centerBin));
% [centeredDir, centeredData_vis]  = alignMtxDir(kernel_pop(:,1:8,:), tgtTimes, param.cardinalDir);
% [~, centeredData_eye]  = alignMtxDir(kernel_pop(:,9:16,:), tgtTimes, param.cardinalDir);
% 
% 
% figure('position',[680   276   846   702]);
% subplot(221);
% thisImage = mean(centeredData_vis,3)';
% crange = max(abs(thisImage(:)));
% imagesc(kerneltlags, centeredDir, thisImage);
% vline(0);
% caxis([-crange crange]);
% mcolorbar(gca,.5);
% set(gca,'ytick', centeredDir);
% xlabel('time from target onset [s]');
% ylabel('relative target direction [deg]');
% 
% subplot(223);
% thisImage = mean(centeredData_eye,3)';
% imagesc(kerneltlags, centeredDir, thisImage);
% crange = max(abs(thisImage(:)));
% caxis([-crange crange]);
% vline(0);
% mcolorbar(gca,.5);
% set(gca,'ytick',centeredDir);
% xlabel('time from eye movement [s]');
% ylabel('relative eye direction [deg]');
% 
% subplot(424);
% plot(kerneltlags, squeeze(kernel_pop(:,17,:)), 'color',[.5 .5 .5]);
% hold on;
% plot(kerneltlags, mean(kernel_pop(:,17,:),3), 'linewidth',2);
% xlabel('time from pupil dilation [s]');
% axis tight
% vline(0);
% 
% subplot(422);
% plot(kerneltlags, squeeze(kernel_pop(:,18,:)), 'color',[.5 .5 .5]);
% hold on;
% plot(kerneltlags, mean(kernel_pop(:,18,:),3), 'linewidth',2);
% xlabel('time from blink onset [s]');
% axis tight
% vline(0);
% 
% subplot(426);
% plot(kerneltlags, squeeze(kernel_pop(:,19,:)), 'color',[.5 .5 .5]);
% hold on;
% plot(kerneltlags, mean(kernel_pop(:,19,:),3), 'linewidth',2);
% xlabel('time from reward on [s]');
% axis tight
% vline(0);
% 
% subplot(428);
% plot(kerneltlags, squeeze(kernel_pop(:,20,:)), 'color',[.5 .5 .5]);
% hold on;
% plot(kerneltlags, mean(kernel_pop(:,20,:),3), 'linewidth',2);
% xlabel('time from punish onset [s]');
% axis tight
% vline(0);
% 
% screen2png('centeredkernel_pop');
% savePaperFigure(gcf,'centeredkernel_pop');



% %% explained variance, correlation
% subplot(211);
% plot(mFiringRate_pop, expval_pop, '.');
% xlabel('mean firing rate [Hz]');
% ylabel('explained variance [%]');
% axis square;
% marginplot;
% 
% subplot(212);
% plot(mFiringRate_pop, corrcoef_pop, '.');
% xlabel('mean firing rate [Hz]');
% ylabel('correlation coef');
% axis square
% marginplot;
% screen2png('mFiringRate_corrcoef_expVar_pop');
% savePaperFigure(gcf,'mFiringRate_corrcoef_expVar_pop');


% %% powerspectrum
% subplot(121);
% semilogy(faxis_common, pspec_parea_pop, 'color',[.5 .5 .5]);
% hold on
% semilogy(faxis_common, mean(pspec_parea_pop,2), 'linewidth',2);
% xlim([0 25]);
% xlabel('Frequency [Hz]'); ylabel('parea PSD');
%
% subplot(122);
% semilogy(faxis_common, pspec_psth_pop, 'color',[.5 .5 .5]);
% hold on
% semilogy(faxis_common, mean(pspec_psth_pop,2), 'linewidth',2);
% xlabel('Frequency [Hz]'); ylabel('psth PSD');
% xlim([0 25]);
%
% screen2png('powerspectrum_pop');


%% direction tuning DOITAGAIN
% subplot(121);
% plot(cmean_pop(1,), cmean_pred_pop, '.');
% squareplot;
% marginplot;
% title('Preferred direction (cirular mean)');
% xlabel('measured [rad]');
% ylabel('fitted [rad]');
%
% subplot(122);
% plot(cvar_pop, cvar_pred_pop, '.');
% squareplot;
% marginplot;
% title('Tuning width (cirular variance)');
% xlabel('measured [rad]');
% ylabel('fitted [rad]');
%
% screen2png('cmean_var_pop');

% %% alpha power 
% corrbins = -1:.04:1;
% for ifreq = 1:4
%     subplot(4,2,2*ifreq-1);
%     histogram(corrcoef_parea_spk_pop(ifreq,:), corrbins);
%     [~,Pspk(ifreq)]=ttest(corrcoef_parea_spk_pop(ifreq,:));
%     vline(0);
%     title(['ttest p:' num2str(Pspk(ifreq))]);
%     
%     subplot(4,2,2*ifreq);
%     histogram(corrcoef_parea_alpha_pop(ifreq,:), corrbins);
%     [~,Palpha(ifreq)]=ttest(corrcoef_parea_alpha_pop(ifreq,:));
%     title(['ttest p:' num2str(Palpha(ifreq))]);
%     vline(0);
% end
% marginplots;
% savePaperFigure(gcf,'corrcoef_parea-spk-alpha');

