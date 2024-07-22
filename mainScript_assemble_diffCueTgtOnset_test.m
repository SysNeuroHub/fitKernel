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


   mampRespByCue_pop = cell(nData,1);
   p_visRespByCue_pop = cell(nData,1);
   p_cueModulationByCue_pop = cell(nData,1);
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
                S = load(saveName, 'kernelInfo','stats_stratifiedByCue','latency_r_cue'); %PSTH_f','predicted_all', 'predicted', ...
                    % 'kernelInfo','t_r','cellclassInfo','mFiringRate','t_cat','mdiffCueFOnset','stddiffCueFOnset'); %param
               
                if isfield(S,'kernelInfo')
                   
                    %id_pop{numel(id_pop)+1} = thisid;
                    id_pop{idata} = thisid;

                    eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
                    %eyeData = load(eyeName,'catEvTimes'); %'eyeData_rmotl_cat','startSaccNoTask'
                    %ddData = load(loadNames{idata},'dd'); %slow to read


                                   
                    %% explained variance for target response averaged across trials
                    %                 [expval_avgtgt(1,1), corr_avgtgt(1,1)] = getExpVal_avgtgt(S.PSTH_f, S.predicted_all, ...
                    %                     eyeData.catEvTimes, S.t_r, [0 0.5], param.cardinalDir, dd);
                    %                 [expval_avgtgt(2:6,1), corr_avgtgt(2:6,1)] = getExpVal_avgtgt(S.PSTH_f, S.predicted, ...
                    %                     eyeData.catEvTimes, S.t_r, [0 0.5], param.cardinalDir, dd);
                    %                 expval_avgtgt_pop = [expval_avgtgt_pop expval_avgtgt];
                    %                 corr_avgtgt_pop = [corr_avgtgt_pop corr_avgtgt];

                    % %% compute time-resolved explained variance
                    % [choiceOutcome] = getChoiceOutcome(dd);
                    % 
                    % expval_trig = []; corr_trig = [];
                    % for ievtype = 1:4
                    %     % 1: success trial
                    %     % 2: failed quiescent trial
                    %     % 3: failed wrong saccade direction
                    %     % 4: saccade outside trials
                    % 
                    %     if ievtype <= 3
                    %         theseEvents = find(choiceOutcome==ievtype);
                    %         eventTimes = eyeData.catEvTimes.tOnset(theseEvents);
                    %     elseif ievtype == 4
                    %         eventTimes = startSaccNoTask;
                    %     end
                    %     [expval_trig(:,:,ievtype), corr_trig(:,:,ievtype), winSamps] = ...
                    %         getExpVal_trig(S.PSTH_f, [S.predicted_all S.predicted], S.t_r, eventTimes, [-0.5 0.5]);
                    % 
                    % %ax_expval_trig(ievtype) = subplot(4,1,ievtype);
                    % %plot(winSamps, squeeze(expval_trig(:,:,ievtype))');
                    % end
                    % expval_trig_pop = cat(4, expval_trig_pop, expval_trig);
                    % corr_trig_pop = cat(4, corr_trig_pop, expval_trig);

                    % % %% number of target trials
                    % %     ntargetTrials_pop{idata} = sum(~isnan(eyeData.catEvTimes.tOnset));
                    % % 
                    % % %% number of total trials
                    % % ntotTrials_pop{idata} = numel(eyeData.catEvTimes.tOnset);
                    % % 
                    % % %% compute gain
                    % % figTWin = [-0.5 0.5];
                    % % onlySuccess = 0;
                    % % respWin = [0.05 0.35]; %[s]
                    % % y_r = cat(2,S.PSTH_f,S.predicted_all);
                    % % gainInfo = getGainInfo(S.t_r, y_r, param.cardinalDir, eyeData.catEvTimes, ...
                    % %     ddData.dd, figTWin, onlySuccess, respWin);
                    % % gainInfo_pop{idata} = gainInfo;
                    % % 
                    % % 
                    % % %% latency
                    % % S = load(saveName, 'latency_bhv','latency_neuro','latency_r','stats_stratified');
                    % % latency_bhv_pop{idata} = S.latency_bhv;
                    % % latency_neuro_pop{idata} = S.latency_neuro;
                    % % latency_r_pop{idata} = S.latency_r;
                    % % stats_stratified_pop{idata} = S.stats_stratified;
                    latency_r_cue_pop{idata} = S.latency_r_cue;

                    %stats_stratifiedByCue_pop{idata} = S.stats_stratifiedByCue;
                    % if ~isempty(S.stats_stratifiedByCue) && isfield(S.stats_stratifiedByCue, 'multcomp_anova')
                    %     stats_stratifiedByCue_pop{idata} = min(S.stats_stratifiedByCue.multcomp_anova{1}(:,6));
                    % else
                    %     stats_stratifiedByCue_pop{idata} = NaN;
                    % end
                    mampRespByCue_pop{idata} = S.stats_stratifiedByCue.mampResp;
                    p_visRespByCue_pop{idata} = S.stats_stratifiedByCue.p_visResp;
                    p_cueModulationByCue_pop{idata} = S.stats_stratifiedByCue.p_cueModulation;
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
% assembly.mFiringRate_pop = mFiringRate_pop; %needed for inclusionCriteria
% assembly.kernelInfo_pop = kernelInfo_pop;
% assembly.corrcoef_pred_spk_pop = corrcoef_pred_spk_pop;
% assembly.cellclassInfo_pop = cellclassInfo_pop;
% assembly.expval_ind_pop = expval_ind_pop;
% assembly.expval_tgt_pop = expval_tgt_pop;
% assembly.corr_avgtgt_pop = corr_avgtgt_pop;
% assembly.expval_avgtgt_pop = expval_avgtgt_pop;
% assembly.ntargetTrials_pop = ntargetTrials_pop;  %needed for inclusionCriteria
% assembly.ntotTrials_pop = ntotTrials_pop;
% assembly.id_pop = id_pop;
% assembly.Rsqadj_pop = Rsqadj_pop; %needed for pickUnitsByClass
% assembly.Rsqadj_tgt_pop = Rsqadj_tgt_pop; 
% assembly.gainInfo_pop = gainInfo_pop;
% assembly.expval_trig_pop = expval_trig_pop;
% assembly.corr_trig_pop = corr_trig_pop;
% assembly.latency_bhv_pop = latency_bhv_pop;
% assembly.latency_neuro_pop = latency_neuro_pop;
% assembly.latency_r_pop = latency_r_pop;
% assembly.stats_stratified_pop = stats_stratified_pop;
% assembly.expval_tgt_rel_pop = expval_tgt_rel_pop;
% assembly.corr_tgt_rel_pop = corr_tgt_rel_pop;
% assembly.corr_tgt_pop = corr_tgt_pop;
% assembly.stats_stratified_pop = stats_stratifiedByCue_pop;
 assembly.latency_r_cue_pop = latency_r_cue_pop;
assembly.mampRespByCue_pop = mampRespByCue_pop;
assembly.p_visRespByCue_pop = p_visRespByCue_pop;
assembly.p_cueModulationByCue_pop = p_cueModulationByCue_pop;
%% 
% save('fitPSTH_pop20220202','avgPupilResp_pop', '-append');
%save(fullfile(saveServer,[saveSuffix_p animal '.mat']),'assembly');
save(fullfile(saveFolder, 'assembly_tmp.mat'),'assembly','param');
assembly = [];
end