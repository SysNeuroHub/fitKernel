%% shot tonset resp w/wo cue of selected functional cell types, detected in pickUnitsByClass.m

[saveServer, rootFolder] = getReady();
load(fullfile(saveServer,'param20230405.mat'),'param');

dataDir = '/mnt/syncitium/Daisuke/cuesaccade_data/figPSTH_pop20231026hugo/';

animal = 'hugo';%'ollie';%'andy';% 'andy' '


%% load gain info <> avg tgt resp [HACK]
load(fullfile(dataDir, 'fitPSTH_pop20231026hugo.mat'),...
    'gainInfo_pop','id_pop','kernel_pop');


%% load kernel info from tgt units
gainKernelName = fullfile(dataDir, 'gainkernel_splt_all.mat');
if exist(gainKernelName, 'file')
    load(gainKernelName);
else
    id_all = cell(1);
    kernel_splt_all = cell(1,10);
    prefDir_all = [];
    iii=1;
    for yyy = 1:3
        switch yyy
            case 1
                year = '2021';
            case 2
                year = '2022';
            case 3
                year = '2023';
        end
        [loadNames, months, dates, channels] = getMonthDateCh(animal, year, rootFolder);
        saveName_splt_pop = [];
        for idata = 1:numel(channels)
            datech = [months{idata} filesep dates{idata} filesep num2str(channels{idata})];
            disp([year ': ' num2str(idata) '/' num2str(numel(channels)) ', ' datech ]);

            saveFolder = fullfile(saveServer, year,animal);%17/6/23
            saveSuffix = [animal replace(datech,filesep,'_') '_linear_rReg'];
            saveName_splt = fullfile(saveFolder, [saveSuffix '_splitPredictor.mat']);
            if exist(saveName_splt,'file')
                try
                    thisKernelInfo = load(saveName_splt, 'kernelInfo');
                    kernel_splt_all(iii,1:10) = thisKernelInfo.kernelInfo.kernel;
                    prefDir_kernel(iii)

                    thisid = [animal '/' year '/' datech];

                    %% load avg resp to tgt
                    [~, thisIDidx] = intersect(id_pop, thisid);
                    if ~isempty(thisIDidx)
                        avgResp_all(:,:,:,iii) = squeeze(gainInfo_pop(thisIDidx).avgTonsetByCue(:,1,:,:));
                        prefDir_all(iii) = gainInfo_pop(thisIDidx).prefDir;
                    end

                    id_all{iii} = thisid;
                    iii = iii+1;
                catch err
                    disp(err);
                    continue;
                end
            end
        end
    end
        winSamps = gainInfo_pop(1).winSamps;
    tlags = thisKernelInfo.kernelInfo.tlags;
    save(gainKernelName, 'kernel_splt_all','avgResp_all','prefDir_all','id_all','winSamps','tlags');
end


%% tgt units
load(fullfile(dataDir, 'pickUnitsByClass.mat'),"funcClass",'nUnits');
for itgt = 1:5
    switch itgt
        case 1
            tgtModality = 'vision';
            tgtID = funcClass.id_v;
        case 2
            tgtModality = 'eyeSpeed';
            tgtID = funcClass.id_es;
        case 3
            tgtModality = 'eyePosition';
            tgtID = funcClass.id_ep;
        case 4
            tgtModality = 'integrator';
            tgtID = funcClass.id_i;
        case 5
            tgtModality = 'all';
            tgtID = funcClass.id_all;
    end

    %thisid = [animal '/' year '/' datech];

    %tgtUnits_s = cellfun(@(y)(replace(y,filesep,'_')), tgtUnits);

    [~,tgtIDidx] = intersect(id_all, tgtID);
    %[~,tgtIDidx_g] = intersect(id_pop, tgtID);

    %kernel_pop: {units kernelType}
    kernel_splt_selected = kernel_splt_all(tgtIDidx,:)';

    avgResp_selected = avgResp_all(:,:,:,tgtIDidx);
    prefDir_selected = prefDir_all(tgtIDidx);


    %% show average kernel before centering
    [f, kernel_avg] = showKernelByCue(kernel_splt_selected, tlags, param.cardinalDir, 0);
    screen2png( ['avgKernel_' animal '_' tgtModality]);


    %% centerring by preferred direction
    tgtRange = [0.05 0.15; ... %vision
        0.03 0.25; ... %eye speed
        -0.1 0.1]; %eye position
    [f, kernel_centered_avg, prefKernelDir] = showKernelByCue(kernel_splt_selected, tlags, param.cardinalDir, 1, tgtRange);
    screen2png(['avgKernel_centered_' animal '_' tgtModality ]);


    %% show individual resp
mResp = squeeze(mean(avgResp_selected,4));
seResp = squeeze(ste(avgResp_selected,4));

dir0 = 1;
dir180 = 5;

figure('position',[0 0 1800 1200]);
ax3(1)=subplot(321);
boundedline(winSamps, squeeze(mResp(dir180,:,1)), squeeze(seResp(dir180,:,1)),'k','transparency',.5);
boundedline(winSamps, squeeze(mResp(dir0,:,1)), squeeze(seResp(dir0,:,1)),'transparency',.5);
title('wo cue'); grid on;

ax3(2)=subplot(322);
boundedline(winSamps, squeeze(mResp(dir180,:,2)), squeeze(seResp(dir180,:,2)),'k','transparency',.5);
boundedline(winSamps, squeeze(mResp(dir0,:,2)), squeeze(seResp(dir0,:,2)),'transparency',.5);
title('w cue');legend('180deg','0deg','location','northwest');
linkaxes(ax3); grid on;

    %% show avg response to tgt stim before centering
    for icue = 1:2
        ax1(icue)=subplot(3,2,icue+2);
        imagesc(winSamps, param.cardinalDir, squeeze(mean(avgResp_selected(:,:,icue,:),4)));
        hline;vline;
        if icue==1
            ylabel('target direction');
            title('wo cue');
        elseif icue==2
            title('w cue');
        end
    end
    linkcaxes(ax1);
    mcolorbar;
        
    %% show avg response to tgt stim after centering
    avgResp_selected_centered = [];
    for iunit = 1:size(avgResp_selected,4)
        [avgResp_selected_centered(:,:,:,iunit), dirAxis] = dealRespByCue(avgResp_selected(:,:,:,iunit), 1, 0, ...
            prefDir_selected(iunit), param.cardinalDir);
    end

    for icue = 1:2
        ax2(icue)=subplot(3,2,icue+4);
        imagesc(winSamps, dirAxis, squeeze(nanmean(avgResp_selected_centered(:,:,icue,:),4)));
        hline;vline;
        if icue==1
            ylabel('centered direction');
        end
    end
    linkcaxes(ax2);
   mcolorbar;
        


 screen2png(['avgResp_' animal '_' tgtModality ]);
    close all;

end
