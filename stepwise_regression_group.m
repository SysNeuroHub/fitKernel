function [selected_groups, selected_cols] = stepwise_regression_group(X, y, groups, groupNames, initial_cols)
    % Perform stepwise regression at a group level.
    %
    % Arguments:
    % X -- matrix, independent variables (n_samples x n_features)
    % y -- column vector, dependent variable (n_samples x 1)
    % groups -- cell array of vectors, groups of column indices
    % groupNames -- cell array of strings, names of the groups
    % initial_cols -- row vector, column indices to start with (optional)
    % threshold_in -- float, p-value threshold for variable inclusion (default: 0.05)
    % threshold_out -- float, p-value threshold for variable exclusion (default: 0.10)
    %
    % Returns:
    % selected_cols -- row vector, final selected column indices
    %
    % 12/7/2023 created
    % 13/7/2023 added option for fitting algorithm
    
    option = 'ridge'; %fitlm
    ridgeParams = 100;
    useGPU = 0;

    if nargin < 4
        error('Please provide groups and groupNames.');
    end
    
    if nargin < 5
        initial_cols = [];
    end
    
    if nargin < 6
        threshold_Rsqadj = 0;
    end
  
    
    % Create a mapping from column index to group index
    group_mapping = zeros(1, size(X, 2));
    
    for group_idx = 1:numel(groups)
        group = groups{group_idx};
        group_mapping(group) = group_idx;
    end
    
    selected_cols = initial_cols;
    selected_groups = [];
    
    while true
        changed = false;
        
        % Forward step
        %excluded_cols = setdiff(1:size(X, 2), selected_cols);
        excluded_groups = setdiff(1:numel(groups), selected_groups);
        best_Rsqadj = 0;
       
        best_group_idx = 0;
        best_group_cols = [];
        
        %for col = excluded_cols
        for group_idx = excluded_groups 
            %group_idx = group_mapping(col);
            group_cols = groups{group_idx};
            
%             if ~isempty(setdiff(group_cols, selected_cols))
%                 continue;
%             end
            
            X_subset = [X(:, selected_cols), X(:, group_cols)];
            X_subset = [ones(size(X_subset, 1), 1), X_subset];
            
            switch option
                case 'fitlm'
                    model = fitlm(X_subset, y);
                    Rsqadj = model.Rsquared.Adjusted;
                case 'ridge'
                    b = rReg(X_subset, y, ridgeParams, useGPU); %from ridgeXs.m
                    yHat = X_subset*b(2:end)+b(1);
                    nPredictors = numel(b);
                    [Rsqadj] = getRsqadj(y, yHat, nPredictors);
            end

            if Rsqadj > best_Rsqadj
                best_Rsqadj = Rsqadj;
                best_group_idx = group_idx;
                best_group_cols = group_cols;
            end
        end
        
        if best_Rsqadj > threshold_Rsqadj
            threshold_Rsqadj = best_Rsqadj;
            selected_cols = [selected_cols, best_group_cols];
            selected_groups = [selected_groups, best_group_idx];
            changed = true;
        end
        
         %% Backward step
        for group_idx = selected_groups 
            group_cols = groups{group_idx};
            
%             if ~isempty(setdiff(group_cols, selected_cols))
%                 continue;
%             end
           
            X_subset = [X(:, setxor(selected_cols, group_cols))];
            X_subset = [ones(size(X_subset, 1), 1), X_subset];
            
            switch option
                case 'fitlm'
                    model = fitlm(X_subset, y);
                    Rsqadj = model.Rsquared.Adjusted;
                case 'ridge'
                    b = rReg(X_subset, y, ridgeParams, useGPU); %from ridgeXs.m
                    yHat = X_subset*b(2:end)+b(1);
                    nPredictors = numel(b);
                    [Rsqadj] = getRsqadj(y, yHat, nPredictors);
            end
            
            if Rsqadj > best_Rsqadj
                best_Rsqadj = Rsqadj;
                best_group_idx = group_idx;
                best_group_cols = group_cols;
            end
        end
        
        if  best_Rsqadj > threshold_Rsqadj
            threshold_Rsqadj = best_Rsqadj;
            selected_cols = setxor(selected_cols, best_group_cols);
            selected_groups = setxor(selected_groups, best_group_idx);
            changed = true;
        end
        
        if ~changed
            break;
        end
    end
    
    % Display selected groups
    disp('Selected groups:');
    for group_idx = 1:numel(groups)
        group_cols = groups{group_idx};
        group_name = groupNames{group_idx};
        if any(ismember(group_cols, selected_cols))
            disp(['- ' group_name]);
        end
    end
end
