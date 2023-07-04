
tgtRange = [0.05 0.15; 0.03 0.25; -0.1 0.1];

fitPars = [];
for col = 1:3
    
    allmatrix = reshape(zeros(size(kernel_pop{col,1})),[],1);
    for row = 1:nrow
        allmatrix(:,row) = reshape(kernel_pop{col,row},1,[]);
    end
    orisize = size(kernel_pop{col,1});
    allmatrix = reshape(allmatrix, orisize(1), orisize(2),[]);
    
    tgtTimes = intersect(find(tlags{col}(:,1)>tgtRange(col,1)), ...
        find(tlags{col}(:,1)<tgtRange(col,2)));
    
    for idata = 1:size(allmatrix,3)
        resp = mean(allmatrix(tgtTimes,:,idata),1);
        
        [fitPars(:,idata,col), fitErr] ...
            = fitoriWrapped(param.cardinalDir, resp,...
            [], [nan nan 0 min(resp) nan],'',20, []);
        %fitPars: [Dp, Rp, Rn, Ro, sigma]
        %		Dp is the preferred direction (bet 0 and 360)
        %		Rp is the response to the preferred direction;
        %		Rn is the response to the opposite direction;
        %		Ro is the background response (useful only in some cases)

    end
end

tuning = abs(squeeze(fitPars(2,:,:)./fitPars(4,:,:)));%looks good!
for col=1:3
    tuned(:,col) = (tuning(:,col)>2);
end
allTuned = find(sum(tuned,2)==3);

%% very few are tuned for all the three modalities(?)
loglog(abs(tuning(:,1)), abs(tuning(:,2)), '.');

prefDir = squeeze(fitPars(1,:,:));


for ii = 1:3
    switch ii
        case 1
            v = [1 2];
        case 2
            v = [1 3];
        case 3
            v = [2 3];
    end
    
    subplot(1,3,ii);
    
    plot(prefDir(:,v(1)), prefDir(:,v(2)), '.','color',[.7 .7 .7]);hold on
    plot(prefDir(allTuned,v(1)), prefDir(allTuned,v(2)), 'b.');
    axis equal square;
    [rho, pval] = circ_corrcc(prefDir(allTuned,v(1))*pi/180, prefDir(allTuned,v(2))*pi/180);
    title(['rho:' num2str(rho) ', pval:' num2str(pval)]);
    xlabel(param.predictorNames{v(1)});
    ylabel(param.predictorNames{v(2)});
end
savePaperFigure(gcf,['tuning_pop_' animal]);

save(['fitPSTH_pop20230619' animal],'fitPars','-append');




