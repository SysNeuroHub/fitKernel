[saveServer, rootFolder] = getReady();
load(fullfile(saveServer,'param20230405.mat'),'param');

animal = 'hugo';%'ollie';%'andy';% 'andy' '
tgtModality = 'all';%'eyeSpeed';



%% load kernel info from tgt units
id_all = cell(1);
kernel_all = cell(1,10);
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
                %kernel_all = [kernel_all thisKernelInfo.kernelInfo];
                kernel_all(iii,1:10) = thisKernelInfo.kernelInfo.kernel;
                thisid = [animal '/' year '/' datech];
                %id_all(1,iii) = [id_all {thisid}];
                id_all{iii} = thisid;
                iii = iii+1;
            catch err
                disp(err);
                continue;
            end
        end
    end
end


%% tgt units
load('pickUnitsByClass_test.mat',"funcClass",'nUnits');
switch tgtModality
    case 'vision'
        tgtID = funcClass.id_v;
    case 'eyeSpeed'
        tgtID = funcClass.id_es;
    case 'eyePosition'
        tgtID = funcClass.id_ep;
    case 'integrator'
        tgtID = funcClass.id_i;
    case 'all'
        tgtID = funcClass.id_all;
end

%thisid = [animal '/' year '/' datech];
 
%tgtUnits_s = cellfun(@(y)(replace(y,filesep,'_')), tgtUnits);

[~,tgtIDidx] = intersect(id_all, tgtID);


for icue = 1:2
    switch icue
        case  1
            suffix = 'wo cue';
        case 2
            suffix = 'w cue';
    end
    cidx = 5*(icue-1)+(1:5);
    kernel_pop = kernel_all(tgtIDidx,cidx)';
    tlags = thisKernelInfo.kernelInfo.tlags(cidx);


    %% show average kernel before centering
    %kernel_pop: {units kernelType}
    [f, kernel_avg] = showKernel3(kernel_pop, tlags, param.cardinalDir, 0);
    screen2png( ['avgKernel_' animal '_' tgtModality '_' suffix]);


    %% centerring by preferred direction
    tgtRange = [0.05 0.15; 0.03 0.25; -0.1 0.1];
    [f, kernel_centered_avg] = showKernel3(kernel_pop, tlags, param.cardinalDir, 1, tgtRange);
    screen2png(['avgKernel_centered_' animal '_' tgtModality '_' suffix]);
end

