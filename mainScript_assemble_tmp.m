% assembly a matfile compiling data across units

%%TODO
%load data
%align number of rows of population data

%% get ready
limSuffix = '';%'_woSuccess';%_woSuccess';
load('/media/daisuke/cuesaccade_data/param20240625.mat');

[saveServer, rootFolder] = getReady();
saveSuffix_p = ['fitPSTH_pop' limSuffix];

%% recorded data
animal = 'hugo';
dataType = 0;%0: each channel, 1: all channels per day
tWin = [0 0.5];%[s]


for yy = 1:3
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
 thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
         regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','03March','18','*_ch20')))));

   nData = length(channels);


   p_hm_pop = cell(nData,1);
   avgAmp_hm_pop = cell(nData,1);
   latency_r_cue_pop = cell(nData,1);
    for idata = 1:length(channels)
        datech = [months{idata} '/' dates{idata} '/' num2str(channels{idata})];
        thisid = [animal '/' year '/' datech];
        disp(thisid);

        saveSuffix = [animal replace(datech,'/','_') '_linear_rReg' limSuffix];

        thisDate = [months{idata} '_' dates{idata}];

        saveFolder = fullfile(saveServer, year,animal);%17/6/23
        saveName = fullfile(saveFolder, [saveSuffix '.mat']);

        if exist(saveName, 'file')

            %datech_pop{idata} = datech;
            try

                %result of mainScript.m and mainScript_latency.m
                S = load(saveName, 'kernelInfo',...
                    'avgAmp_hm', 'p_hm'); %PSTH_f','predicted_all', 'predicted', ...
                    % 'kernelInfo','t_r','cellclassInfo','mFiringRate','t_cat','mdiffCueFOnset','stddiffCueFOnset'); %param
               
                if isfield(S,'kernelInfo')
                   
                    %id_pop{numel(id_pop)+1} = thisid;
                    id_pop{idata} = thisid;

                    eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
                    %eyeData = load(eyeName,'catEvTimes'); %'eyeData_rmotl_cat','startSaccNoTask'
                    %ddData = load(loadNames{idata},'dd'); %slow to read

                    avgAmp_hm_pop{idata} = S.avgAmp_hm;
                    p_hm_pop{idata} = S.p_hm;
                    %output bhv latency - c-t interval rank correlation
                    %why some units have no field?

                    errorIDs{idata} = 0;
                    S = []; eyeData = []; ddData = [];
                end
            catch err
                errorIDs{idata} = 1;
                disp(err.message);
                disp(err.stack);
                    S = []; eyeData = []; ddData = [];
                    %continue;
            end
        end
    end
    %dataByYear = dataByYear(~isnan(dataByYear));
    %alldata = [alldata dataByYear(:)];

%% make a giant structure
assembly.p_hm_pop = p_hm_pop;
assembly.avgAmp_hm_pop = avgAmp_hm_pop;

%% 
% save('fitPSTH_pop20220202','avgPupilResp_pop', '-append');
%save(fullfile(saveServer,[saveSuffix_p animal '.mat']),'assembly');
save(fullfile(saveFolder, 'assembly_tmp.mat'),'assembly','param');
assembly = [];
end