function [model, selectedGroups] = stepwiseGrouped(X, y, groups)

%% NOT FUNCTIONAL
% TODO:
% omit intruction term between variables
% select groups
% incorpolate temporal delay for each variable

    numGroups = numel(groups);
    selectedGroups = cell(1, numGroups);
    
    for i = 1:numGroups
        groupVars = groups{i};
        mdl = stepwiseglm(X(:, groupVars), y, 'linear', 'Criterion', 'AIC');
        inmodel = mdl.Formula.InModel;
        if any(inmodel)
            selectedGroups{i} = groupVars;
        end
    end
    
    selectedGroups = selectedGroups(~cellfun('isempty', selectedGroups));
    
    allVars = cat(2, selectedGroups{:});
    model = stepwiseglm(X(:, allVars), y, 'linear', 'Criterion', 'AIC');
end
