addpath(genpath('C:\Users\dshi0006\git'))
% load('C:\Users\dshi0006\Downloads\hugo_oephysdata_ch23.mat', ...
%     'ch','dd','ephysdata');

saveFolder = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\Daisuke\cuesaccade_data';

%% recorded data
animal = 'hugo';
rootFolder = '\\storage.erc.monash.edu.au\shares\R-MNHS-Physio\SysNeuroData\Monash Data\Joanita\2021/cuesaccade_data/';
dataType = 0;%0: each channel, 1: all channels per day


[loadNames, months, dates, channels] = getMonthDateCh(animal, rootFolder);

% to obtain index of specified month&date&channel
%find(1-cellfun(@isempty, regexp(loadNames, regexptranslate('wildcard','12December\09\*_ch32.mat'))))

%% omit data
% no saccade response
% low spontaneous firing
% low number of successful trials

% parameters
% param.marginSize = 50;%40; %frames
% param.psth_sigma = .05;%0.1; %0.025;%0.02;%0.01;%0.05;%[s] %gaussian filter
% param.dt_r = 0.02; %for securing memory for kernel fitting %0.025;%0.5;
% param.lagRange = [-.5 .5];%[-.4 .3];%[-1 1];%[-10 20]; %temporal window for kernel estimation [s]
% param.ridgeParams = 100;%[0 1e-1 1 1e2 1e3]; %10
% % visualize = 0;
% param.predictorNames = {'vision','eyeposition','pdiam','blink'};
% param.figTWin = [-0.5 0.5]; %temporal window for peri-event traces [s]
% param.respWin = [0.03 0.25]; %temporal window to compute direction tuning
% param.pareaTh = 3;
% param.pareaDiffTh = 5;
% param.cutoffFreq = 0.1;
% param.evName = 'tOnset';%'cOnset';
% cardinalDir = linspace(0,360,9); %direction for saccade and target
% param.cardinalDir = cardinalDir(1:end-1);
% ncDirs = length(param.cardinalDir);

idataIdx = 1;

for idata = 1:length(channels) %1061;%865;%
    datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
    disp(datech);
    
    saveSuffix = [animal replace(datech,'/','_')];
    
    thisDate = [months{idata} '_' dates{idata}];
    
    saveName = fullfile(saveFolder, [saveSuffix '.mat']);
    
    if exist(saveName, 'file')
        
        datech_pop{idataIdx} = datech;
        load(saveName);%, 'kernelInfo','psth_all','mDir',...
%             'seDir','mDir_pred','seDir_pred','cvar','cmean', ...
%             'cvar_pred','cmean_pred','mFiringRate',...
%             'param','eyeData_rmotl_cat',...
%             'pspec_psth','pspec_parea','faxis_psth','faxis_parea',...
%             'corrCoef_parea_alpha','corrCoef_parea_spk');
        
%         avgPupilResp_pop(:,:,:,idataIdx) = avgPupilResp;%[dl/cs x variables x time] pupilOnsetsResp.m

        kerneltlags = kernelInfo.tlags;

        %faxis_common = 1:50; %[Hz]
        %pspec_psth_interp = interp1(faxis_psth, pspec_psth, faxis_common);
        %pspec_parea_interp = interp1(faxis_parea, pspec_parea, faxis_common);
        
        mFiringRate_pop(idataIdx) = mFiringRate;
        kernel_pop(:,:,idataIdx) = kernelInfo.kernel;
        expval_pop(idataIdx) = kernelInfo.expval;
        corrcoef_pop(idataIdx) = kernelInfo.corrcoef;
        %pspec_psth_pop(:,idataIdx) = pspec_psth_interp;
        %pspec_parea_pop(:,idataIdx) = pspec_parea_interp;
        cmean_pop(:,:,idataIdx) = cmean;
        %cmean_pred_pop(idataIdx) = cmean_pred;
        cvar_pop(:,:,idataIdx) = cvar;
        %cvar_pred_pop(idataIdx) = cvar_pred;
        %corrcoef_parea_alpha_pop(:,idataIdx) = corrCoef_parea_alpha;
        %corrcoef_parea_spk_pop(:,idataIdx) = corrCoef_parea_spk;
        %avgTgtResp_pop(:,:,:,idataIdx) = avgTgtResp;
        
        R = corrcoef(PSTH_f, predicted_all);
        corrcoef_pred_spk_pop(idataIdx) = R(1,2); 
        avgSaccResp = [];
        for idir = 1:length(param.cardinalDir)
            theseTr = find(sortedSaccLabels == param.cardinalDir(idir));
            avgSaccResp(idir,:,:) = mean(singleSaccResp(theseTr,:,:),1);
        end
        %FIX: avgSaccResp includes nans
        
        avgSaccResp_pop(:,:,:,idataIdx) = avgSaccResp; %[saccDir, kernel?, time, channel]
        idataIdx = idataIdx+1;
        clear mFiringRate kernelInfo
    end 
end
% save('fitPSTH_pop20220202','avgPupilResp_pop', '-append');
save('fitPSTH_pop20220209','mFiringRate_pop','kernel_pop','expval_pop','corrcoef_pop',...
    'cmean_pop','cvar_pop','datech_pop','avgSaccResp_pop','avgSaccResp');



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


%% target resp aligned to preferred target direction
pavgTgtResp = permute(avgTgtResp_pop, [3 1 4 2]);
[centeredDir, centeredData_tgt]  = alignMtxDir(pavgTgtResp, tgtTimes, param.cardinalDir);
%[time x direction x channels x predictors]
avgCentTgt = squeeze(nanmean(centeredData_tgt,3)); 
crange = prctile(avgCentTgt(:),[1 99]);
for ipred = 1:size(avgCentTgt,3)
    subplot(6,1,ipred);
    imagesc(winSamps, centeredDir, squeeze(avgCentTgt(:,:,ipred))');
    ylabel(psthNames{ipred})
    caxis(crange);
    vline(0);
    mcolorbar(gca,.5);
end
xlabel('time from target onset [s]');
screen2png('centeredTgt_pop');
savePaperFigure(gcf,'centeredTgt_pop');


%% kernel aligned to preferred direction
%centerBin = 4;
%centeredDir = 180/pi*circ_dist(pi/180*param.cardinalDir, pi/180*param.cardinalDir(centerBin)); 
[centeredDir, centeredData_vis]  = alignMtxDir(kernel_pop(:,1:8,:), tgtTimes, param.cardinalDir);
[~, centeredData_eye]  = alignMtxDir(kernel_pop(:,9:16,:), tgtTimes, param.cardinalDir);


figure('position',[680   276   846   702]);
subplot(221);
thisImage = mean(centeredData_vis,3)';
crange = max(abs(thisImage(:)));
imagesc(kerneltlags, centeredDir, thisImage);
vline(0);
caxis([-crange crange]);
mcolorbar(gca,.5);
set(gca,'ytick', centeredDir);
xlabel('time from target onset [s]');
ylabel('relative target direction [deg]');

subplot(223);
thisImage = mean(centeredData_eye,3)';
imagesc(kerneltlags, centeredDir, thisImage);
crange = max(abs(thisImage(:)));
caxis([-crange crange]);
vline(0);
mcolorbar(gca,.5);
set(gca,'ytick',centeredDir);
xlabel('time from eye movement [s]');
ylabel('relative eye direction [deg]');

subplot(424);
plot(kerneltlags, squeeze(kernel_pop(:,17,:)), 'color',[.5 .5 .5]);
hold on;
plot(kerneltlags, mean(kernel_pop(:,17,:),3), 'linewidth',2);
xlabel('time from pupil dilation [s]');
axis tight
vline(0);

subplot(422);
plot(kerneltlags, squeeze(kernel_pop(:,18,:)), 'color',[.5 .5 .5]);
hold on;
plot(kerneltlags, mean(kernel_pop(:,18,:),3), 'linewidth',2);
xlabel('time from blink onset [s]');
axis tight
vline(0);

subplot(426);
plot(kerneltlags, squeeze(kernel_pop(:,19,:)), 'color',[.5 .5 .5]);
hold on;
plot(kerneltlags, mean(kernel_pop(:,19,:),3), 'linewidth',2);
xlabel('time from reward on [s]');
axis tight
vline(0);

subplot(428);
plot(kerneltlags, squeeze(kernel_pop(:,20,:)), 'color',[.5 .5 .5]);
hold on;
plot(kerneltlags, mean(kernel_pop(:,20,:),3), 'linewidth',2);
xlabel('time from punish onset [s]');
axis tight
vline(0);

screen2png('centeredkernel_pop');
savePaperFigure(gcf,'centeredkernel_pop');



%% explained variance, correlation
subplot(211);
plot(mFiringRate_pop, expval_pop, '.');
xlabel('mean firing rate [Hz]');
ylabel('explained variance [%]');
axis square;
marginplot;

subplot(212);
plot(mFiringRate_pop, corrcoef_pop, '.');
xlabel('mean firing rate [Hz]');
ylabel('correlation coef');
axis square
marginplot;
screen2png('mFiringRate_corrcoef_expVar_pop');
savePaperFigure(gcf,'mFiringRate_corrcoef_expVar_pop');


%% powerspectrum
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

%%
corrbins = -1:.04:1;
for ifreq = 1:4
    subplot(4,2,2*ifreq-1);
    histogram(corrcoef_parea_spk_pop(ifreq,:), corrbins);
    [~,Pspk(ifreq)]=ttest(corrcoef_parea_spk_pop(ifreq,:));
    vline(0);
    title(['ttest p:' num2str(Pspk(ifreq))]);

    subplot(4,2,2*ifreq);
    histogram(corrcoef_parea_alpha_pop(ifreq,:), corrbins);
    [~,Palpha(ifreq)]=ttest(corrcoef_parea_alpha_pop(ifreq,:));
    title(['ttest p:' num2str(Palpha(ifreq))]);
    vline(0);
end
marginplots;
savePaperFigure(gcf,'corrcoef_parea-spk-alpha');

