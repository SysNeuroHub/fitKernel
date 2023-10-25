function [fig, pval, rho] = showScatterTriplets(values, names, valueRange, selectedIDs, corrType)
%fig = showScatterTriplets(values, names, valueRange, selectedIDs)
%
%fig = showScatterTriplets(values, names, valueRange, selectedIDs, 'circular')
% computes circular correlation. For this values should be degree not radian
%
% TODO: replace circ_corrcc with non-parametric one (p.145 of Fisher's book)
if nargin < 5
    corrType = 'linear';
end
if nargin < 4
    selectedIDs = [];
end
if nargin < 3
    valueRange = [];
end

fig = figure('position',[0 0 1500 500]);

for ii = 1:3
    switch ii
        case 1
            v = [2 3];
        case 2
            v = [2 1];
        case 3
            v = [1 3];
    end
    
    subplot(1,3,ii);
    xvalues = values(v(1),:);
    yvalues = values(v(2),:);
    if strcmp(corrType, 'linear')
        [rho, pval] = corr(xvalues', yvalues', 'type','Spearman');
    elseif strcmp(corrType, 'circular')
        [rho, pval] = circ_corrcc(xvalues'*pi/180, yvalues'*pi/180);
    end
    if ~isempty(valueRange)
        xvalues(xvalues<valueRange(1))=valueRange(1);
        xvalues(xvalues>valueRange(2))=valueRange(2);
        yvalues(yvalues<valueRange(1))=valueRange(1);
        yvalues(yvalues>valueRange(2))=valueRange(2);
    end
    plot(xvalues, yvalues,'k.'); hold on;
    if ~isempty(selectedIDs)
        c = autumn(numel(selectedIDs));
        scatter(xvalues(selectedIDs), yvalues(selectedIDs),10,c,'filled');
    end
    title(['rho:' num2str(rho) ', pval:' num2str(pval)])
    xlabel(names{v(1)}); ylabel(names{v(2)});
    axis equal square;
    if ~isempty(valueRange)
        xlim(valueRange);ylim(valueRange);
    end
    set(gca,'tickdir','out');
end

