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

for idata = 940;%439;%1:length(channels) %1061;%865;%
    datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
    disp(datech);
    
    saveSuffix = [animal replace(datech,'/','_')];
    
    thisDate = [months{idata} '_' dates{idata}];
    
    %     load(fullfile(saveFolder,['predictorInfo_' thisDate '.mat']), 'predictorInfo');
    %     load(fullfile(saveFolder,['eyeCat_' thisDate '.mat']), 'eyeData_rmotl_cat',...
    %         'onsets_cat','meta_cat','blinks','outliers','pspec_parea','faxis_parea','t_tr');
    
    saveName = fullfile(saveFolder, [saveSuffix '.mat']);
    
    if exist(saveName, 'file')
        
        datech_pop{idataIdx} = datech;
        load(saveName, 'kernelInfo','psth_all','mDir',...
            'seDir','mDir_pred','seDir_pred','cvar','cmean', ...
            'cvar_pred','cmean_pred','mFiringRate',...
            'param','eyeData_rmotl_cat',...
            'pspec_psth','pspec_parea','faxis_psth','faxis_parea');
        
        
        faxis_common = 1:50; %[Hz]
        pspec_psth_interp = interp1(faxis_psth, pspec_psth, faxis_common);
        pspec_parea_interp = interp1(faxis_parea, pspec_parea, faxis_common);
        
        mFiringRate_pop(idataIdx) = mFiringRate;
        kernel_pop(:,:,idataIdx) = kernelInfo.kernel;
        expval_pop(idataIdx) = kernelInfo.expval;
        corrcoef_pop(idataIdx) = kernelInfo.corrcoef;
        pspec_psth_pop(:,idataIdx) = pspec_psth_interp;
        pspec_parea_pop(:,idataIdx) = pspec_parea_interp;
        cmean_pop(idataIdx) = cmean;
        cmean_pred_pop(idataIdx) = cmean_pred;
        cvar_pop(idataIdx) = cvar;
        cvar_pred_pop(idataIdx) = cvar_pred;
        idataIdx = idataIdx+1;
        clear mFiringRate kernelInfo pspec_psth pspec_parea mDir mDir_pred
    end 
end
save('fitPSTH_pop','mFiringRate_pop','kernel_pop','expval_pop','corrcoef_pop',...
    'pspec_psth_pop','pspec_parea_pop','cmean_pop','cmean_pred_pop','cvar_pop',...
    'cvar_pred_pop','datech_pop','faxis_common');

%% kernel
kerneltlags = kernelInfo.tlags;
subplot(221);
thisImage = mean(kernel_pop(:,1:8,:),3)';
crange = max(abs(thisImage(:)));
imagesc(kerneltlags, param.cardinalDir, thisImage);
caxis([-crange crange]);
set(gca,'ytick',param.cardinalDir);
xlabel('time from target onset [s]');
ylabel('target direction [deg]');

subplot(223);
thisImage = mean(kernel_pop(:,9:16,:),3)';
imagesc(kerneltlags, param.cardinalDir, thisImage);
crange = max(abs(thisImage(:)));
caxis([-crange crange]);
set(gca,'ytick',param.cardinalDir);
xlabel('time from eye movement [s]');
ylabel('eye direction [deg]');

subplot(224);
plot(kerneltlags, squeeze(kernel_pop(:,17,:)), 'color',[.5 .5 .5]);
hold on;
plot(kerneltlags, mean(kernel_pop(:,17,:),3), 'linewidth',2);
xlabel('time from pupil dilation [s]');
axis tight

subplot(222);
plot(kerneltlags, squeeze(kernel_pop(:,18,:)), 'color',[.5 .5 .5]);
hold on;
plot(kerneltlags, mean(kernel_pop(:,18,:),3), 'linewidth',2);
xlabel('time from blink onset [s]');
axis tight

screen2png('kernel_pop');


%% kernel aligned to preferred direction
centerBin = 4;
centeredDir = 180/pi*circ_dist(pi/180*param.cardinalDir, pi/180*param.cardinalDir(centerBin)); 
tgtTimes = intersect(find(kerneltlags>0.03), find(kerneltlags<0.25));
centeredData_vis = zeros(size(kernel_pop(:,1:8,:)));
for idata = 1:size(kernel_pop,3)
    thisData = kernel_pop(:,1:8,idata);
    [~,prefDir] = max(mean(thisData(tgtTimes,:,:),1));
    centeredData_vis(:,:,idata) = circshift(thisData, centerBin - prefDir, 2);
end
centeredData_eye = zeros(size(kernel_pop(:,1:8,:)));
for idata = 1:size(kernel_pop,3)
    thisData = kernel_pop(:,9:16,idata);
    [~,prefDir] = max(mean(thisData(tgtTimes,:,:),1));
    centeredData_eye(:,:,idata) = circshift(thisData, centerBin - prefDir, 2);
end

subplot(221);
thisImage = mean(centeredData_vis,3)';
crange = max(abs(thisImage(:)));
imagesc(kerneltlags, centeredDir, thisImage);
caxis([-crange crange]);
set(gca,'ytick', centeredDir);
xlabel('time from target onset [s]');
ylabel('target direction relative to preferred [deg]');

subplot(223);
thisImage = mean(centeredData_eye,3)';
imagesc(kerneltlags, param.cardinalDir, thisImage);
crange = max(abs(thisImage(:)));
caxis([-crange crange]);
set(gca,'ytick',centeredDir);
xlabel('time from eye movement [s]');
ylabel('eye direction relative to preferred [deg]');

screen2png('centeredkernel_pop');



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


%% powerspectrum
subplot(121);
semilogy(faxis_common, pspec_parea_pop, 'color',[.5 .5 .5]);
hold on
semilogy(faxis_common, mean(pspec_parea_pop,2), 'linewidth',2);
xlim([0 25]);
xlabel('Frequency [Hz]'); ylabel('parea PSD');

subplot(122);
semilogy(faxis_common, pspec_psth_pop, 'color',[.5 .5 .5]);
hold on
semilogy(faxis_common, mean(pspec_psth_pop,2), 'linewidth',2);
xlabel('Frequency [Hz]'); ylabel('psth PSD');
xlim([0 25]);

screen2png('powerspectrum_pop');


%% direction tuning
subplot(121);
plot(cmean_pop, cmean_pred_pop, '.');
squareplot;
marginplot;
title('Preferred direction (cirular mean)');
xlabel('measured [rad]');
ylabel('fitted [rad]');

subplot(122);
plot(cvar_pop, cvar_pred_pop, '.');
squareplot;
marginplot;
title('Tuning width (cirular variance)');
xlabel('measured [rad]');
ylabel('fitted [rad]');

screen2png('cmean_var_pop');

