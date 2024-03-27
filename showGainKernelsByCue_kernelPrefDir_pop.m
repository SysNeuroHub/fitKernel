[saveServer, rootFolder] = getReady();
load(fullfile(saveServer,'param20230405.mat'),'param');

animal = 'hugo';%'ollie';%'andy';% 'andy' '
tgtModality = 'all';%'eyeSpeed';

dataDir = '/mnt/syncitium/Daisuke/cuesaccade_data/figPSTH_pop20231026hugo/';


%% load gain info <> avg tgt resp [HACK]
load(fullfile(dataDir,'fitPSTH_pop20231026hugo.mat'),...
    'gainInfo_pop','id_pop','kernel_pop','tlags');

%preferred direction & amplitude of each kernel modality
tgtRange = [0.05 0.15; 0.03 0.25; -0.1 0.1];
prefDirOption = 0;
[kernelPrefDir, kernelAmp] = getKernelPrefDirAmp(kernel_pop, tlags, tgtRange, param.cardinalDir, prefDirOption);


%% load kernel info from tgt units
gainKernelName = fullfile(dataDir,'gainkernel_splt_all.mat');
if exist(gainKernelName, 'file')
    load(gainKernelName);
else
    id_all = cell(1);
    kernel_splt_all = cell(1,10);
    prefDir_all = [];
    latency_bhv_all = cell(1);
    latency_neuro_all = cell(1);
    latency_r_all = [];
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
            saveName = fullfile(saveFolder, [saveSuffix '.mat']);

            if exist(saveName_splt,'file') && exist(saveName, 'file')
                try
                    thisKernelInfo = load(saveName_splt, 'kernelInfo');
                    kernel_splt_all(iii,1:10) = thisKernelInfo.kernelInfo.kernel;
                    %prefDir_kernel(iii)

                    thisid = [animal '/' year '/' datech];

                    %% load avg resp to tgt
                    [~, thisIDidx] = intersect(id_pop, thisid);
                    if ~isempty(thisIDidx)
                        avgResp_all(:,:,:,iii) = squeeze(gainInfo_pop(thisIDidx).avgTonsetByCue(:,1,:,:));
                        prefDir_all(iii) = gainInfo_pop(thisIDidx).prefDir;
                    end

                    %% load latency
                    load(saveName,'latency_bhv','latency_neuro','latency_r');
                    latency_bhv_all{iii} = latency_bhv;
                    latency_neuro_all{iii} = latency_neuro;
                    latency_r_all(iii) = latency_r;


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
    save(gainKernelName, 'kernel_splt_all','avgResp_all','prefDir_all','id_all','winSamps',...
        'latency_bhv_all','latency_neuro_all','latency_r_all');
end

%% exclude NG units
load(fullfile(dataDir,'pickUnitsByClass.mat'),"funcClass",'nUnits');

%% tgt units
highAmp = true;
if highAmp
    suffix = '_highAmp';
    param.ampTh = .5;
else
    suffix = '';
end
for itgtModality = 1:3
    switch itgtModality
        case 1
            tgtModality = 'vision';
        case 2
            tgtModality = 'eyeSpeed';
        case 3
            tgtModality = 'eyePosition';
    end
    for iprefDir =1%:8

        prefDir_q = quantizeDir(kernelPrefDir(:,itgtModality), param.cardinalDir);
        if ~highAmp
            tgtID = id_pop(prefDir_q== param.cardinalDir(iprefDir));
        else
            tgtID = id_pop(prefDir_q== param.cardinalDir(iprefDir) & ...
                kernelAmp(:,itgtModality) > param.ampTh);
        end

        tgtID = intersect(funcClass.id_all, tgtID);%exclude NG units

        [~,tgtIDidx] = intersect(id_all, tgtID);

        %kernel_pop: {units kernelType}
        %kernel_splt_selected = kernel_splt_all(tgtIDidx,:)';
        %tlags = thisKernelInfo.kernelInfo.tlags;

        avgResp_selected = avgResp_all(:,:,:,tgtIDidx);
        prefDir_selected = prefDir_all(tgtIDidx);



        %% show individual resp
        mResp = squeeze(mean(avgResp_selected,4));
        seResp = squeeze(ste(avgResp_selected,4));

        dir0 = 1;
        dir180 = 5;

        figure('position',[0 0 1800 1200]);
        ax3(1)=subplot(321);
        boundedline(winSamps, squeeze(mResp(dir180,:,1)), squeeze(seResp(dir180,:,1)),'k','transparency',.5);
        boundedline(winSamps, squeeze(mResp(dir0,:,1)), squeeze(seResp(dir0,:,1)),'transparency',.5);
        title(['wo cue, n=' num2str(numel(tgtIDidx))]); grid on;

        ax3(2)=subplot(322);
        boundedline(winSamps, squeeze(mResp(dir180,:,2)), squeeze(seResp(dir180,:,2)),'k','transparency',.5);
        boundedline(winSamps, squeeze(mResp(dir0,:,2)), squeeze(seResp(dir0,:,2)),'transparency',.5);
        title('w cue'); %legend('180deg','0deg','location','northwest');
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

        screen2png(['avgResp_' animal '_' tgtModality 'kernel_prefDir' num2str(param.cardinalDir(iprefDir)) ...
            suffix]);
        close all;

        %% latency
        latency_neuro_selected = latency_neuro_all(tgtIDidx);
        latency_neuro_selected = cat(1, latency_neuro_selected{:});
        latency_bhv_selected = latency_bhv_all(tgtIDidx);
        latency_bhv_selected = cat(1, latency_bhv_selected{:});
        latency_r_selected = latency_r_all(tgtIDidx);

        figure('visible','on')
        subplot(121);
        hist3([latency_bhv_selected latency_neuro_selected],'CDataMode','auto','FaceColor','interp',...
            'linestyle','none','edges',{0:.02:.5 0:.02:.5});
        xlabel('behavioural latency'); ylabel('neural latency');
        view(2);
        axis square tight;
        set(gca,'tickdir','out');

        subplot(122);
        histogram(latency_r_selected,'BinEdges',-0.7:0.1:0.7); axis square tight;        
        vline(0); set(gca,'tickdir','out');
        xlabel('latency correlation'); ylabel('# units');
        screen2png(['latency_selected_' tgtModality num2str(param.cardinalDir(iprefDir)) ...
            suffix '.png']);
        close all;

    end
end


