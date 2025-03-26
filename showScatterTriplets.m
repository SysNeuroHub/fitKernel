function [fig, pval, rho] = showScatterTriplets(values, names, valueRange, selectedIDs, corrType, animalid)
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
figPosition = [0 0 670 223];%[0 0 600 200];
msize = 8;

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

        subplot(1,3,ii);
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
            c = 1-summer(numel(selectedIDs));
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
            squareplots(ax, valueRange);
        end
        if aa==2 && ii ==1
            legend('M1','','M2');
        end
    end
end

