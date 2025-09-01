function adjustLegendPatchWidth(lgd, factor)
% adjustLegendPatchWidth Adjusts histogram/bar patch width in legend (R2023a+)
%
%   adjustLegendPatchWidth(LGD, FACTOR) rescales the width of the
%   patch icons in the legend by FACTOR (0.5 = half width, 2 = double).
%
%   LGD is the legend handle returned by legend(...).

    if nargin < 2
        factor = 0.5; % default shrink
    end

    % Access legend entries (undocumented, but stable since R2018a)
    entries = lgd.EntryContainer.NodeChildren;

    for k = 1:numel(entries)
        icon = entries(k).Icon;
        if isempty(icon) || ~isvalid(icon)
            continue
        end

        % Histogram/bar entries typically have one child Patch
        iconChildren = icon.Children;
        for c = 1:numel(iconChildren)
            if isa(iconChildren(c), 'matlab.graphics.primitive.Patch')
                xd = iconChildren(c).XData;
                yd = iconChildren(c).YData;

                if numel(xd) == 5 && numel(yd) == 5
                    % Standard rectangle [x1 x2 x2 x1 x1]
                    xMid = mean([min(xd), max(xd)]);
                    halfWidth = (max(xd) - min(xd)) * factor / 2;
                    newX = [xMid-halfWidth, xMid+halfWidth, ...
                            xMid+halfWidth, xMid-halfWidth, xMid-halfWidth];
                    iconChildren(c).XData = newX;
                end
            end
        end
    end
end
