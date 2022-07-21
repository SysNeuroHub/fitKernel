function [f, psth_snippet, pdiam_snippet, dist_snippet, taxis_snippet] = ...
    pupilFigure(dd, eyeData_rmblk_tr, psth_tr, evName, tWindow)

%f = pupilFigure(dd, eyeData_rmblk_tr, psth_tr, evName)

% evName = 'cOnset';
% evName = 'tOnset';
% evName = 'fOnset';
% evName = 'fixationdt';
% evName = 'srt';

nTypes = size(psth_tr,2);

fs_eye = eyeData_rmblk_tr(1).fs(1);
taxis_snippet = tWindow(1):1/fs_eye:tWindow(2);

theseTr = intersect(find(dd.successTrials), find(~isnan(dd.(evName))));
eventTimes = dd.(evName);

% pdiam_prctile = getPupilDiameter(eyeData_rmblk_cat);
            
pdiam_snippet = zeros(length(taxis_snippet), length(theseTr));
dist_snippet = zeros(length(taxis_snippet), length(theseTr));
psth_snippet = zeros(length(taxis_snippet), length(theseTr), nTypes);
for itr = 1:length(theseTr)
    thisTr = theseTr(itr);
    taxis_signal = eyeData_rmblk_tr(thisTr).t - eventTimes(thisTr);
    
    pdiam_snippet(:,itr) = interp1(taxis_signal, sqrt(eyeData_rmblk_tr(thisTr).parea), taxis_snippet);
    
    dist_snippet(:,itr) =  interp1(taxis_signal, ...
        sqrt(eyeData_rmblk_tr(thisTr).x.^2+eyeData_rmblk_tr(thisTr).y.^2), taxis_snippet);
    
    for itype = 1:nTypes
        psth_snippet(:,itr, itype) = interp1(taxis_signal, psth_tr{thisTr, itype}, taxis_snippet);
    end
end

f = figure('position',[0 0 1000 1000]);
subplot(411);
plot(taxis_snippet,pdiam_snippet, 'color',[.5 .5 .5]);
hold on
errorbar(taxis_snippet, nanmedian(pdiam_snippet,2), ...
    1/sqrt(length(theseTr))*nanstd(pdiam_snippet,[],2), 'linewidth',2);
ylabel('pdiam');
grid on
ylim(prctile(pdiam_snippet(:),[20 80]));
title([num2str(length(theseTr)) ' trials']);

subplot(412);
plot(taxis_snippet,dist_snippet, 'color',[.5 .5 .5]);
hold on
errorbar(taxis_snippet, nanmedian(dist_snippet,2), ...
    1/sqrt(length(theseTr))*nanstd(dist_snippet,[],2), 'linewidth',2);
ylabel('distance');
grid on
axis tight

subplot(413);

plot(taxis_snippet, squeeze(psth_snippet(:,:,1)), 'color',[.5 .5 .5]);
hold on
errorbar(taxis_snippet, nanmedian(psth_snippet(:,:,1),2), ...
    1/sqrt(length(theseTr))*nanstd(psth_snippet(:,:,1),[],2), 'linewidth',2);
c = squeeze(nanmedian(psth_snippet(:,:,1),2));
ylim([.9*min(c) 1.1*max(max(c),min(c)+1)]);
ylabel('observed psth');

subplot(414);
for itype = 1+1:nTypes
    errorbar(taxis_snippet, nanmedian(psth_snippet(:,:,itype),2), ...
        1/sqrt(length(theseTr))*nanstd(psth_snippet(:,:,itype),[],2), 'linewidth',2);
    hold on;
end
ylabel('predicted psth');
grid on
axis tight

xlabel(['time from ' evName ' [s]']);

%screen2png(['peri' evName '_predicted']);
end