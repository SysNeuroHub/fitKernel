function [X, tlags, groups_wtlag, groupNames_wtlag] = getPredictorsDelayed(timeVec, ...
    predictor, lagRanges, npredVars, predictorNames)
%X = getPredictorsDelayed(timeVec, predictor, lagRanges)
%
% from ridgeXs
%
% INPUTS:
% preditor: nVar x times


nt = length(timeVec);
fs = 1/median(diff(timeVec));

nVar = size(predictor,1);
nVar_lag = size(lagRanges,1);

lags = round(min(lagRanges(:,1))*fs):round(max(lagRanges(:,2))*fs);
localLags = cell(nVar,1);
for iVar = 1:nVar_lag
    localLags{iVar} = round(lagRanges(iVar,1)*fs):round(lagRanges(iVar,2)*fs);
end

groups = [];
groups{1} = 1:npredVars(1);
for iii=2:numel(npredVars)
    groups{iii} = groups{iii-1}(end)+(1:npredVars(iii));
end
    
tlags = lags/fs;
nLags = length(lags);

groups_wtlag = cell(1,numel(groups));
groupNames_wtlag = cell(1,numel(groups));
 X = zeros(nt, nVar*nLags);
    nLags_tmp = nLags;
    for ibhv = 1:nVar
        
        for ig = 1:numel(groups)
            if ~isempty(find(groups{ig} == ibhv))
                igroup = ig; continue;
            end
        end
        
        for iLag = 1:nLags
            %X(:, iLag) = circshift( predictor(:), lags(iLag)-1);
            if ismember(lags(iLag),localLags{ibhv}) || nVar_lag == 1
                X(:, iLag+(ibhv-1)*nLags) = circshift( predictor(ibhv,:), lags(iLag)-1);
            else
                %values outside the lagRange is substituted by 0 ... is
                %this correct?
                X(:, iLag+(ibhv-1)*nLags) = zeros(1,nt); %nans did not work
            end
            
           groups_wtlag{igroup} = [groups_wtlag{igroup} iLag+(ibhv-1)*nLags];
           groupNames_wtlag{igroup} = predictorNames{igroup};
        end
    end