

 [okunits, mfiringRateOK, expvalOK, ntargetTrOK, ptonsetRespOK] ...
    = inclusionCriteria(mFiringRate_pop, expval_pop, ntargetTrials_pop, PtonsetResp_pop, param);


expvalTrace = trace2Event(expvalOK); %RED
ntargetTrTrace = trace2Event(ntargetTrOK); %GREEN
ptonsetRespTrace = trace2Event(ptonsetRespOK); %BLUE

figure('position',[0 0 700 1400]);
ax(1)=subplot(311);
plot(ntotTrials_pop);hold on
plot(ntargetTrials_pop);
legend('nTotTrials','ntargetTrials');
hline(param.ntargetTrTh);
%vbox(expvalTrace(:,1),expvalTrace(:,2),gca,[1 .5 .5 .5]);
 vbox(ntargetTrTrace(:,1),ntargetTrTrace(:,2),gca,[.5 1 .5 .5]);
% vbox(ptonsetRespTrace(:,1),ptonsetRespTrace(:,2),gca,[.5 .5 1 .5]);
title(animal);
axis padded

ax(2)=subplot(312);
semilogy(PsaccResp_pop);hold on
semilogy(PtonsetResp_pop);
legend('psaccresp','ptonsetresp');
hline(param.ptonsetRespTh);
%vbox(expvalTrace(:,1),expvalTrace(:,2),gca,[1 .5 .5 .5]);
%vbox(ntargetTrTrace(:,1),ntargetTrTrace(:,2),gca,[.5 1 .5 .5]);
vbox(ptonsetRespTrace(:,1),ptonsetRespTrace(:,2),gca,[.5 .5 1 .5]);
ylim([1e-10 1e0])
axis padded

ax(3)=subplot(313);
yyaxis left
plot(expval_pop);ylim([0 20]);
hline(param.expvalTh);
yyaxis right
plot(corrcoef_pop);hold on
plot(corrcoef_pred_spk_pop);
legend('expval','corr','corr pred spk');
vbox(expvalTrace(:,1),expvalTrace(:,2),gca,[1 .5 .5 .5]);
%vbox(ntargetTrTrace(:,1),ntargetTrTrace(:,2),gca,[.5 1 .5 .5]);
%vbox(ptonsetRespTrace(:,1),ptonsetRespTrace(:,2),gca,[.5 .5 1 .5]);
xlabel('unit index');
axis padded 

linkaxes(ax,'x');
xlim([1 numel(expval_pop)])
screen2png(['validUnits_' animal]);