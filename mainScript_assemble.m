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
 % thisdata = find(1-cellfun(@isempty, regexp(loadNames, ...
 %         regexptranslate('wildcard',fullfile(rootFolder, year, 'cuesaccade_data','08August','02','*_ch21')))));

   nData = length(channels);

mFiringRate_pop = cell(nData,1); %needed for inclusionCriteria
kernelInfo_pop = cell(nData,1);
%kernel_pop = cell(nData,1);
%expval_pop = cell(nData,1); %needed for inclusionCriteria
%corrcoef_pop = cell(nData,1);
corrcoef_pred_spk_pop = cell(nData,1);
%PtonsetResp_pop = cell(nData,1); %needed for inclusionCriteria
%PsaccResp_pop = cell(nData,1);
cellclassInfo_pop = cell(nData,1);
expval_ind_pop = cell(nData,1);
expval_tgt_pop = cell(nData,1);
corr_avgtgt_pop = cell(nData,1);
expval_avgtgt_pop = cell(nData,1);
ntargetTrials_pop = cell(nData,1);  %needed for inclusionCriteria
ntotTrials_pop = cell(nData,1);
id_pop = cell(nData,1);
%Rsqadj_pop = cell(nData,1); %needed for pickUnitsByClass
%Rsqadj_tgt_pop = cell(nData,1); 
gainInfo_pop = cell(nData,1);
expval_trig_pop = cell(nData,1);
corr_trig_pop = cell(nData,1);
latency_r_pop = cell(nData,1);
latency_neuro_pop = cell(nData,1);
stats_stratified_pop = cell(nData,1);
corr_tgt_pop = cell(nData,1);
corr_tgt_rel_pop = cell(nData,1);
expval_tgt_rel_pop = cell(nData,1);
mdiffCueFOnset_pop = cell(nData,1);
stddiffCueFOnset_pop = cell(nData,1);
diffCueFOnset_pop = cell(nData,1);
errorIDs= cell(1);
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
                S = load(saveName, 'PSTH_f','predicted_all', 'predicted', ...
                    'kernelInfo','t_r','cellclassInfo','mFiringRate','t_cat'); %param
               
                if isfield(S,'kernelInfo')
                    
                    
                    % % %retrieve just once
                    % % tlags = S.kernelInfo.tlags;
                    % % % param = S.param;
                    % % 
                    % % mFiringRate_pop{idata} = S.mFiringRate;
                    % % %kernel_pop{idata} = S.kernelInfo.kernel;
                    % % %expval_pop{idata} = S.kernelInfo.expval;
                    % % %corrcoef_pop{idata} = S.kernelInfo.corrcoef;
                    % % kernelInfo_pop{idata} = S.kernelInfo;
                    % % 
                    % % R = corrcoef(S.PSTH_f, S.predicted_all);
                    % % corrcoef_pred_spk_pop{idata} = R(1,2);
                    % %id_pop{numel(id_pop)+1} = thisid;
                    id_pop{idata} = thisid;

                    eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
                    eyeData = load(eyeName,'catEvTimes', 'onsets_cat'); %'eyeData_rmotl_cat','startSaccNoTask'
                    % % ddData = load(loadNames{idata},'dd'); %slow to read

                    %% Rsq adj of subjset of variables
                    %%predictorInfoName = fullfile(saveFolder,['predictorInfo_' animal thisDate '.mat']);
                    % % predData = load(predictorInfoName, 'predictorInfo');
                    nPredictorNames = numel(param.predictorNames);
                    % % nsub=4;
                    % % Rsqadj = zeros(nsub,1);
                    % % for jj = 1:nsub
                    % %     switch jj
                    % %         case 1 %full model
                    % %             tgtGroups = 1:nPredictorNames;
                    % %         case 2 %omit vision %added 26/10/2023
                    % %             tgtGroups = setxor(1:nPredictorNames, 1);
                    % %         case 3 %omit eye speed
                    % %             tgtGroups = setxor(1:nPredictorNames, 2);
                    % %         case 4 %omit eye position
                    % %             tgtGroups = setxor(1:nPredictorNames, 3);
                    % %     end
                    % %     [Rsqadjusted,rr,r0] = fitSubset(S.PSTH_f, predData.predictorInfo, ...
                    % %         tgtGroups, param);%, idxTgtOnsets);
                    % % 
                    % %     Rsqadj(jj) = Rsqadjusted;
                    % % end
                    % % Rsqadj_pop{idata} = Rsqadj;
                    % % predData = [];

                    % % %% Rsq_adjusted computed from pre-target
                    % % % periods - NG. can be < 0
                    % % nPredictors_all = sum(cellfun(@numel, S.kernelInfo.kernel));
                    % % Rsqadj_tgt = nan(1, 4);
                    % % [Rsqadj_tgt(1,1)] = getRsqadj_tgt(S.PSTH_f, S.predicted_all, ...
                    % %     eyeData.catEvTimes, S.t_r, tWin, nPredictors_all);
                    % % for jj = 1:3
                    % %     nPredictors = numel(S.kernelInfo.kernel{jj});
                    % %     [Rsqadj_tgt(1,jj+1)] = getRsqadj_tgt(S.PSTH_f, S.predicted(:,jj), ...
                    % %         eyeData.catEvTimes, S.t_r, tWin, nPredictors);
                    % % end
                    % % Rsqadj_tgt_pop{idata} = Rsqadj_tgt;

                    %% response to target & saccade
                    %PtonsetResp_pop{idata} = S.cellclassInfo.PtonsetResp;                  
                    %PsaccResp_pop = [PsaccResp_pop S.cellclassInfo.PsaccResp];
                    % % cellclassInfo_pop{idata} = S.cellclassInfo;

                    %% explained variance per kernel
                    expval_ind = nan(size(S.predicted,2)+1,1);
                    expval_ind(1,1) = getExpVal(S.PSTH_f, S.predicted_all);
                    % expval_ind(1,1) = getExpVal(S.PSTH_f-mean(S.PSTH_f), ...
                    %     S.predicted_all-mean(S.predicted_all));
                    for ivar = 1:nPredictorNames
                        expval_ind(ivar+1,1) = getExpVal(S.PSTH_f, S.predicted(:,ivar));
                        % expval_ind(ivar+1,1) = getExpVal(S.PSTH_f-mean(S.PSTH_f), ...
                        %     S.predicted(:,ivar)-mean(S.predicted(:,ivar)));
                    end
                    expval_ind_pop{idata} = expval_ind;

                    %% explained variance for target response
                    expval_tgt = zeros(nPredictorNames, 1);
                    corr_tgt = zeros(nPredictorNames, 1);
                    [expval_tgt(1,1), corr_tgt(1,1)] = ...
                        getExpVal_tgt(S.PSTH_f, S.predicted_all, eyeData.catEvTimes, S.t_r, tWin);
                    [expval_tgt(2:nPredictorNames+1,1), corr_tgt(2:nPredictorNames+1,1)] = ...
                        getExpVal_tgt(S.PSTH_f, S.predicted, eyeData.catEvTimes, S.t_r, tWin);

                    expval_tgt_pop{idata} = expval_tgt;
                   expval_tgt_rel_pop{idata} = 100*expval_tgt(2:4)./expval_tgt(1);
                   corr_tgt_pop{idata} = corr_tgt;
                   corr_tgt_rel_pop{idata} = 100*corr_tgt(2:4)./corr_tgt(1);

                   %% cue-tgt interval
                   [diffCueFOnset, mdiffCueFOnset,stddiffCueFOnset] = ...
                        getDiffCueTgtOnset(eyeData.onsets_cat, eyeData.catEvTimes); %3/6/24
                    diffCueFOnset_pop{idata} = diffCueFOnset;
                    mdiffCueFOnset_pop{idata} = mdiffCueFOnset;
                    stddiffCueFOnset_pop{idata} = stddiffCueFOnset;


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
assembly.expval_tgt_rel_pop = expval_tgt_rel_pop;
assembly.corr_tgt_rel_pop = corr_tgt_rel_pop;
assembly.corr_tgt_pop = corr_tgt_pop;
assembly.mdiffCueFOnset_pop = mdiffCueFOnset_pop;
assembly.stddiffCueFOnset_pop = stddiffCueFOnset_pop;

%% 
% save('fitPSTH_pop20220202','avgPupilResp_pop', '-append');
%save(fullfile(saveServer,[saveSuffix_p animal '.mat']),'assembly');
save(fullfile(saveFolder, 'assembly_tmp.mat'),'assembly','param');
assembly = [];
end