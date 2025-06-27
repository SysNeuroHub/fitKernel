function  [startSaccNoTask_spkOkUCue, endSaccNoTask_spkOkUCue, saccDirNoTask_spkOkUCue, dirIndexNoTask_spkOkUCue] = ...
    getSaccNoTask(t_cat, catEvTimes, eyeData_rmotl_cat, blinks, outliers, t_tr, onsets_cat, spkOkTrials, param)
% [startSaccNoTask_spkOkUCue, endSaccNoTask_spkOkUCue, saccDirNoTask_spkOkUCue, dirIndexNoTask_spkOkUCue] = ...
    % getSaccNoTask(t_cat, catEvTimes, eyeData_rmotl_cat, blinks, outliers, t_tr, onsets_cat, spkOkTrials, param)

tOnset = catEvTimes.tOnset;
cOnset = catEvTimes.cOnset; %choice onset not cue %% WHY THIS CONDITION??
validEvents = intersect(find(~isnan(tOnset)), find(~isnan(cOnset)));
tOnset = tOnset(validEvents);
cOnset = cOnset(validEvents);

tcOnset_trace = event2Trace(t_cat, [tOnset; cOnset], 2*0.5);
excEventT_cat = (tcOnset_trace + blinks + outliers > 0); %28/1/22
[startSaccNoTask, endSaccNoTask] = selectSaccades(catEvTimes.saccadeStartTimes, ...
    catEvTimes.saccadeEndTimes, t_cat, excEventT_cat);%param.minSaccInterval); %SLOW
saccNoTaskTrace = event2Trace(t_cat, [startSaccNoTask endSaccNoTask]);
spkOkUCueTrace_tmp = getIncludeTrace(t_cat, t_cat, t_tr, onsets_cat, spkOkTrials);
saccNoTask_spkOkUCue_Trace = saccNoTaskTrace.*spkOkUCueTrace_tmp;

[~, startSaccNoTask_spkOkUCue, endSaccNoTask_spkOkUCue] = trace2Event(saccNoTask_spkOkUCue_Trace, t_cat);
[saccDirNoTask_spkOkUCue, dirIndexNoTask_spkOkUCue] = getSaccDir(startSaccNoTask_spkOkUCue, endSaccNoTask_spkOkUCue, ...
    eyeData_rmotl_cat, param.cardinalDir);