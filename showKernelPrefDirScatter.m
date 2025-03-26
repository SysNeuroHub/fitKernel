function [f, pval] = showKernelPrefDirScatter(kernel_pop, tlags, tgtRange, param, animalid_pop)
%[f, pval] = showKernelPrefDirs(kernel_pop, tlags, tgtRange, param, animalid_pop)

nUnits = size(kernel_pop,2);
figPosition = [0 0 600 400];

if nargin < 5
    animalid_pop = ones(1, nUnits);
end

[prefDir, amp] = getKernelPrefDirAmp(kernel_pop, tlags, tgtRange, param.cardinalDir);
%tuned = amp>param.ampTh;
%tuned = 1:numel(okunits);

f = figure('position',figPosition);
for aa = 1:numel(unique(animalid_pop))
    if aa==1
        asymbol = 'o'; histColor = [.5 .5 .5];
    elseif aa==2
        asymbol = 'square'; histColor = [.1 .1 .1];
    end

    for ii = 1:3
        switch ii
            case 1
                v = [1 2];
            case 2
                v = [1 3];
            case 3
                v = [2 3];
        end
        %doubleTuned = find(tuned(:,v(1))+tuned(:,v(2))==2);
        
        %% scatter plot
        % subplot(3,2,2*ii-1);
         subplot(2,3,ii);
        plot(prefDir(animalid_pop==aa,v(1)), prefDir(animalid_pop==aa,v(2)), asymbol);hold on
        %plot(prefDir(doubleTuned,v(1)), prefDir(doubleTuned,v(2)), 'b.');
        squareplot;
        %[rho, pval] = circ_corrcc(prefDir(doubleTuned,v(1))*pi/180, prefDir(doubleTuned,v(2))*pi/180);
        [rho, pval] = circ_corrcc(prefDir(animalid_pop==aa,v(1))*pi/180, prefDir(animalid_pop==aa,v(2))*pi/180);
        title(['rho:' num2str(rho) ', pval:' num2str(pval)]);
        xlabel(param.predictorNames{v(1)});
        ylabel(param.predictorNames{v(2)});
        set(gca,'tickdir','out');

        %% histogram
        %subplot(3,2,2*ii)
        subplot(2,3,ii+3);
        histogram(prefDir(animalid_pop==aa,v(1)) - prefDir(animalid_pop==aa,v(2)),-180:5:180,...
            'FaceColor', histColor); hold on;
        %histogram(prefDir(doubleTuned,v(1)) - prefDir(doubleTuned,v(2)), -180:5:180,  'facecolor', 'b');
        axis tight square; box off;
        ylim([0 35]);
        if ii==1
            ylabel('# Units');
        end
       if aa == numel(unique(animalid_pop))
            vline(0);
        end
        xlabel([param.predictorNames{v(1)} '-' param.predictorNames{v(2)}]);
        %pval = circ_medtest(prefDir(doubleTuned,v(1)) - prefDir(doubleTuned,v(2)),0);
        %pval = circ_medtest(prefDir(:,v(1)) - prefDir(:,v(2)),0);
        %title(['pval:' num2str(pval)]);
        set(gca,'tickdir','out','xtick',[-90 0 90]);

    end
end
legend('M1','M2');