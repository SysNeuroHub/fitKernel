%return exemplery cells 2024 March

[saveServer, rootFolder] = getReady();
load(fullfile(saveServer,'param20230405_copy.mat'),'param');

animal = 'hugo';%'ollie';%'andy';% 'andy' '
tgtModality = 'all';%'eyeSpeed';

dataDir = '/mnt/syncitium/Daisuke/cuesaccade_data OBS/figPSTH_pop20231026hugo/';


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


        id_tgtModality{itgtModality, iprefDir} = id_pop(tgtIDidx);
    

    end
end

save(fullfile(dataDir,'exampleCells.mat'), 'id_tgtModality');
