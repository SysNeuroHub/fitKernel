function f = drawRectanglesFromMatrix(matrix, valueRange, f, subplotPos)
    % Get the dimensions of the matrix
    [rows, cols] = size(matrix);

    % Create a new figure
    if isempty(valueRange)
        valueRange = [min(matrix(:)) max(matrix(:))];
    end
    minValue = valueRange(1);
    maxValue = valueRange(2);
    
    if  nargin < 3 || isempty(f)
        f = figure;
    end
    
    if nargin < 4
        subplotPos = [ 1 1 1];
    end
    
     % Create or select the specified subplot
    subplot(subplotPos(1), subplotPos(2), subplotPos(3));

    ax = gca; % Get the current axes
    ax.Units = 'normalized'; % Set the axes units to normalized
    
    % Set the initial aspect ratio of the figure window
    initialAspectRatio = cols / rows;

    % Set the aspect ratio mode to manual
    ax.PlotBoxAspectRatioMode = 'manual';
    ax.DataAspectRatioMode = 'manual';

    % Define the callback function for resizing the figure window
    resizeFcn = @(~,~) resizeFigure(matrix, initialAspectRatio, minValue, maxValue);

    % Set the figure's ResizeFcn to the defined callback function
    set(gcf, 'ResizeFcn', resizeFcn);
    
    % Call the function to initially resize the figure
    resizeFigure(matrix, initialAspectRatio, minValue, maxValue);
end

function resizeFigure(matrix, initialAspectRatio, minValue, maxValue)
    % Get the current figure and axes handles
    fig = gcf;
    ax = gca;

    % Get the size of the current figure window
    figPos = fig.Position;
    figWidth = figPos(3);
    figHeight = figPos(4);

    % Compute the aspect ratio of the figure window
    currentAspectRatio = figWidth / figHeight;

    % Calculate the size of each rectangle based on the aspect ratio
    if currentAspectRatio >= initialAspectRatio
        rectWidth = figWidth / size(matrix, 2);
        rectHeight = rectWidth / currentAspectRatio;
    else
        rectHeight = figHeight / size(matrix, 1);
        rectWidth = rectHeight * currentAspectRatio;
    end

    % Clear the axes
    cla;

    % Iterate over each element in the matrix
    for row = 1:size(matrix, 1)
        for col = 1:size(matrix, 2)
            % Compute the position of the rectangle
            x = (col - 1) * rectWidth;
            y = figHeight - (row * rectHeight);

            % Get the value at the current position
            value = matrix(row, col);

            if value >= minValue %skip values below minvalue
                % Map the value to a color based on the minimum and maximum values
                color = (value - minValue) / (maxValue - minValue);
                color = [color, color, color];
               
                % Draw the rectangle
                rectangle('Position', [x, y, rectWidth, rectHeight], ...
                    'FaceColor', color, 'EdgeColor','none');
            end
            
        end
    end

    % Set the aspect ratio and axis limits
    ax.DataAspectRatio = [1, 1, 1];
    ax.PlotBoxAspectRatio = [currentAspectRatio, 1, 1];
    ax.XLim = [0, figWidth];
    ax.YLim = [0, figHeight];

    % Turn off the axis labels
    axis off;

    % Create a colorbar
    colorbar('delete');
    [~,c] = mcolorbar(gca);
    c.Ticks = [0 1];
    c.TickLabels = [minValue maxValue];
    colormap(c,'gray');
end
