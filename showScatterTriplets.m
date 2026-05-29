function [fig, pval, rho] = showScatterTriplets(values, names, valueRange, selectedIDs, corrType, animalid, show3D)
%fig = showScatterTriplets(values, names, valueRange, selectedIDs)
% INPUTS:
% values: [3 x points]
% names: {3x1}
%
%fig = showScatterTriplets(values, names, valueRange, selectedIDs, 'circular')
% computes circular correlation. For this values should be degree not radian
%
% TODO: replace circ_corrcc with non-parametric one (p.145 of Fisher's book)

nUnits = size(values,2);
msize = 8;
if nargin < 7
    show3D = false;
end
if nargin < 6
    animalid = ones(1, nUnits);
end
if nargin < 5
    corrType = 'linear';
end
if nargin < 4
    selectedIDs = [];
end
if nargin < 3
    valueRange = [];
end

units_selectedIDs = zeros(1, nUnits);
units_selectedIDs(selectedIDs) = 1;

if show3D
    figPosition = [0 0 1200 223];
    nCols = 4;
else
    figPosition = [0 0 900 223];%[0 0 600 200];
    nCols = 3;
end
fig = figure('position', figPosition);

for aa = 1:numel(unique(animalid))
    switch aa
        case 1
            asymbol = 'o'; acolor = 'k';%[.5 .5 .5];%'b';
        case 2
            asymbol = 'diamond'; acolor =  'k';%[.5 .5 .5];%'r';
    end

    for ii = 1:3
        switch ii
            case 1
                v = [2 3]; xcolor = 'g'; ycolor='b';
            case 2
                v = [2 1]; xcolor = 'g'; ycolor='r';
            case 3
                v = [1 3]; xcolor = 'r'; ycolor = 'b';
        end

        subplot(1,nCols,ii);
        xvalues = values(v(1), :);
        yvalues = values(v(2), :);

        if strcmp(corrType, 'linear')
            [rho, pval] = corr(xvalues(animalid == aa)', yvalues(animalid == aa)', 'type','Spearman');
        elseif strcmp(corrType, 'circular')
            [rho, pval] = circ_corrcc(xvalues(animalid == aa)'*pi/180, yvalues(animalid == aa)'*pi/180);
        end

        if ~isempty(valueRange)
            xvalues(xvalues(animalid == aa)<valueRange(1))=nan;%valueRange(1);
            xvalues(xvalues(animalid == aa)>valueRange(2))=nan;%valueRange(2);
            yvalues(yvalues(animalid == aa)<valueRange(1))=nan;%valueRange(1);
            yvalues(yvalues(animalid == aa)>valueRange(2))=nan;%valueRange(2);
        end
        %plot(xvalues, yvalues, 'k.'); hold on;
        s0 = scatter(xvalues(animalid==aa), yvalues(animalid==aa), msize, asymbol);
        s0.MarkerEdgeColor = acolor;
        hold on;
        if sum(units_selectedIDs.*animalid == aa) > 0
            %c = autumn(numel(selectedIDs));
            c = 1-summer(sum(units_selectedIDs.*animalid == aa));
            s = scatter(xvalues(units_selectedIDs.*animalid == aa), yvalues(units_selectedIDs.*animalid == aa), ...
                msize, c, 'filled', asymbol);%,'linewidth',2);
            %s.MarkerEdgeColor = 'k';
        end
        title(['rho:' num2str(rho) ', pval:' num2str(pval)])
        xlabel(names{v(1)}); ylabel(names{v(2)});
        axis equal square;
        if ~isempty(valueRange)
            xlim(valueRange);ylim(valueRange);
        end

        ax = gca;
        ax.XColor = xcolor;
        ax.YColor = ycolor;

        set(ax,'tickdir','out');
        if aa==2
            squareplot(ax, valueRange);
        end
        if aa==2 && ii ==1
            legend('M1','','M2');
        end
    end

    if show3D
        subplot(1,4,4);
        zvalues = values(3, :);

        s0 = scatter3(xvalues(animalid==aa), yvalues(animalid==aa), zvalues(animalid==aa), msize, asymbol);
        s0.MarkerEdgeColor = acolor;
        hold on;
        if sum(units_selectedIDs.*animalid == aa) > 0
            %c = autumn(numel(selectedIDs));
            c = 1-summer(numel(selectedIDs));
            % s = scatter(xvalues(units_selectedIDs.*animalid == aa), yvalues(units_selectedIDs.*animalid == aa), ...
            %     msize, c, 'filled', asymbol);%,'linewidth',2);
            s = scatter3(xvalues(units_selectedIDs.*animalid == aa), yvalues(units_selectedIDs.*animalid == aa), ...
                zvalues(units_selectedIDs.*animalid == aa), msize, c, 'filled', asymbol);%,'linewidth',2);
        end
        %        title(['rho:' num2str(rho) ', pval:' num2str(pval)])
        xlabel(names{1}); ylabel(names{2});zlabel(names{3});
        axis equal square;
        if ~isempty(valueRange)
            xlim(valueRange);ylim(valueRange);zlim(valueRange);
        end

        ax = gca;
        set(ax,'tickdir','out');
    end
end

end

function ax = squareplot(ax, xyrange)
%ax = squareplot(ax)
%makes a plot figure square and add a diagonal line
%ax = square(ax, xyrange)
%lets use specified range to plot
% currently does NOT work well after marginplot
% currently only tested with figure plot plot command
% 2017/10/12 created
%2018/2/19 added 2nd input.
% TODO: if image, minimum = minimum + 0.5;

if nargin<1
    ax=gca;
else
    axes(ax);
end

if nargin < 2 || isempty(xyrange)
    axis tight;
    xyrange = [ax.XLim ax.YLim];
end

hold on;
h=line([min(xyrange) max(xyrange)],[min(xyrange) max(xyrange)],'color','k','linestyle','--');
h.HandleVisibility = 'off'; %20/6/22

uistack(h,'bottom');
xlim([min(xyrange) max(xyrange)]);
ylim([min(xyrange) max(xyrange)]);

axis square;
set(ax,'YTick',ax.XTick);
set(gca,'tickdir','out');
end

