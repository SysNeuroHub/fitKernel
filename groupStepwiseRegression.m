function [selectedGroups, selectedVariables] = groupStepwiseRegression(X, y, groups, groupNames, criterion)
    numGroups = length(groups);
    selectedGroups = {};
    selectedVariables = {};
    
    % Stepwise regression for each group
    for i = 1:numGroups
        currentGroup = groups{i};
        currentGroupName = groupNames{i};
        
        fprintf('Performing stepwise regression for Group %d: %s\n', i, currentGroupName);
        
        % Variables within the current group
        currentX = X(:, currentGroup);
        
        % Initialize the current selected variables within the group
        currentSelectedVars = [];
        
        % Stepwise regression within the group
        while true
            % Remaining variables in the current group
            remainingVars = setdiff(currentGroup, currentSelectedVars);
            
            % Calculate the criterion for adding variables
            testVars = [currentSelectedVars, remainingVars];
            if ~isempty(testVars)
                [~, ~, ~, ~, stats] = regress(y, currentX(:, testVars));
                RSS = stats(1);
                addCriterion = RSS;
            else
                addCriterion = [];
            end
            % Calculate the criterion for removing variables
            testVars = setdiff(currentSelectedVars, currentSelectedVars);
            if ~isempty(testVars)
                [~, ~, ~, ~, stats] = regress(y, currentX(:, testVars));
                RSS = stats(1);
                removeCriterion = RSS;
            else
                removeCriterion = [];
            end
            
            % Determine the variable to add or remove based on the criterion
            if isempty(currentSelectedVars)
                selectedVar = remainingVars;
            else
                if isempty(addCriterion) && isempty(removeCriterion)
                    % No more variables to add or remove
                    break;
                elseif ~isempty(addCriterion)
                    selectedVar = currentSelectedVars;
                elseif ~isempty(removeCriterion)
                    selectedVar = remainingVars;
                else
                    if addCriterion <= removeCriterion
                        selectedVar = remainingVars;
                    else
                        selectedVar = currentSelectedVars;
                    end
                end
            end
            
            % Update the selected variables within the group
            if isempty(currentSelectedVars)
                currentSelectedVars = selectedVar;
            else
                if ismember(selectedVar, currentSelectedVars)
                    % Remove the selected variable
                    currentSelectedVars = setdiff(currentSelectedVars, selectedVar);
                else
                    % Add the selected variable
                    currentSelectedVars = [currentSelectedVars, selectedVar];
                end
            end
        end
        
        % Store the selected group and variables
        selectedGroups{i} = currentGroupName;
        selectedVariables{i} = currentSelectedVars;
    end
end
