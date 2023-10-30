[saveServer, rootFolder] = getReady();

load(fullfile(saveServer,'param20230405.mat'),'param');


theseIDs = {'hugo/2021/09September/01/25',... %vision driven
    'hugo/2022/03March/10/20',... %eye speed driven
    'hugo/2022/07July/29/19'}; %eye position driven

Rsqadj = nan(numel(theseIDs),4);
%SSE = nan(numel(theseIDs),3);
%SSR = nan(numel(theseIDs),3);
%SST = nan(numel(theseIDs),3);
rr = [];
for ii = 1:numel(theseIDs)

    aaa = split(theseIDs{ii},'/');
    animal = aaa{1};
    year = aaa{2};
    months = aaa{3};
    dates = aaa{4};
    channels = aaa{5};

    saveFigFolder = fullfile(saveServer, '20230713',year,animal);
    mkdir(saveFigFolder);

    datech = [months filesep dates filesep channels];
    saveSuffix = [animal replace(datech,filesep,'_') ];%'_cue'];

    saveFolder = fullfile(saveServer, year, animal);%17/6/23
    saveName = fullfile(saveFolder, [saveSuffix '.mat']);

    load(saveName, 'PSTH_f');%,'predicted_all', 'predicted','kernelInfo'...
    %,'t_r','param','t_cat','dd');
    %psthNames = cat(2,{'psth','predicted_all'},param.predictorNames);

    thisDate = [months '_' dates];
    eyeName = fullfile(saveFolder,['eyeCat_' animal thisDate '.mat']);
    load(eyeName,'catEvTimes','startSaccNoTask','saccDirNoTask');

    load(fullfile(saveFolder,['predictorInfo_' animal thisDate '.mat']), ...
        'predictorInfo');

    %loadName = getCuesaccadeName(rootFolder, theseIDs{ii});
    %load(loadName,'dd');

    %y_r = cat(2,PSTH_f,predicted_all, predicted);

    %% fitlm using different groups of variables

    groups = [];
    groups{1} = 1:predictorInfo.npredVars(1);
    for iii=2:predictorInfo.nPredictors
        groups{iii} = groups{iii-1}(end)+(1:predictorInfo.npredVars(iii));
    end

    %%
    lagRanges = []; %lag ranges for all variables
    for igroup = 1:size(param.lagRange,1)
        lagRanges = cat(1, lagRanges, repmat(param.lagRange(igroup,:),[predictorInfo.npredVars(igroup) 1]));
    end
    y = PSTH_f;
    [X, tlags, groups_wtlag, groupNames_wtlag] = getPredictorsDelayed(predictorInfo.t_r,...
        predictorInfo.predictors_r, lagRanges, ...
        predictorInfo.npredVars, ...
        param.predictorNames);

    selected_groups{ii} = stepwise_regression_group(X,y,groups_wtlag, groupNames_wtlag);


    for jj = 1:4
        switch jj
            case 1 %full model
                tgtGroups = 1:5;
            case 2 %omit eye speed
                tgtGroups = setxor(1:5, 2);
            case 3 %omit eye position
                tgtGroups = setxor(1:5, 3);
            case 4 %omit vision %25/10/2023
                tgtGroups = setxor(1:5, 1);
        end

        [Rsqadjusted,rr,r0] = fitSubset(PSTH_f, predictorInfo, tgtGroups, param);

        %% stats
        Rsqadj(ii,jj) = Rsqadjusted;%mdl.Rsquared.Adjusted;
        % SSE(ii,jj) = mdl.SSE;
        % SSR(ii,jj) = mdl.SSR;
        % SST(ii,jj) = mdl.SST;
    end
end
% for kk = 1:5
%     subplot(2,3,kk)
%     switch kk
%         case 1
%             thisimage = Rsq; tname = 'Rsquare';
%         case 2
%             thisimage = Rsqadj; tname = 'Rsquare adjusted';
%         case 3
%             thisimage = SSE;  tname = 'SSE';
%         case 4
%             thisimage = SSR;  tname = 'SSR';
%         case 5
%             thisimage = SST;  tname = 'SST';
%     end
%     imagesc(thisimage); colorbar;
%     set(gca,'ytick',1:3,'yticklabel',theseIDs)
%     set(gca,'xtick',1:3,'xticklabel',{'full mdl','wo es','wo ep'})
%     title(tname);
% end
% screen2png('modelComparison_test')

%modality w smallest value = driven by that modality
imagesc(Rsqadj(:,2:4)./Rsqadj(:,1));
set(gca,'Yticklabel',{'vision','eyespeed','eyeposition'});
set(gca,'XTickLabel',{'eyespeed','eyeposition', 'vision'});

%need to limit the analysis window to peristim onset
%otherwise the vision kernel's explanatory power is weak