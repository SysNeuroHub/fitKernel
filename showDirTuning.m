function f = showDirTuning(t_r, y_r, catEvTimes, cardinalDir, param)
%NOT YET IMPLEMENTED
dirs = unique(dd.targetloc);

%         theseTr = find(dd.successTrials);
%         [mDir, winSamps_sacc, singleDirResp, sortedDirLabels, uniqueDirLabels] ...
%             = eventLockedAvg(cat(1,PSTH_f',predicted_all), ...
%             predictorInfo.t_r, catEvTimes.tOnset(theseTr), pi/180*dd.targetloc(theseTr), param.figTWin);
%         cvar = circ_var(sortedDirLabels, singleDirResp);
%         cmean = circ_mean(sortedDirLabels, singleDirResp);

[mDir, sdDir, cmean, cvar, nTrials_dtune] = dirTuning(psth_tr, eyeData_rmotl_tr, dd, param.respWin);
[mDir_pred, sdDir_pred, cmean_pred, cvar_pred ] = dirTuning(psth_predicted_all_tr, eyeData_rmotl_tr, dd, param.respWin);
seDir = sdDir./sqrt(nTrials_dtune);
seDir_pred = sdDir_pred./sqrt(nTrials_dtune);

f = figure('position',[0 0 1000 400]);
errorbar(dirs/180*pi, mDir, seDir, 'b');hold on
errorbar(dirs/180*pi, mDir_pred, seDir_pred, 'g');
legend('observed','fitted','location','southeast');
tname = sprintf('direction tuning %d-%d[ms]\nobserved cvar: %.2f, predicted cvar: %.2f',...
    1e3*param.respWin(1), 1e3*param.respWin(2),cvar, cvar_pred);
title(tname);

%screen2png(fullfile(saveFigFolder,['dirTuning_' saveSuffix]));