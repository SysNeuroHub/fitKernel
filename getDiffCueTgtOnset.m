function [diffCueTgtOnset, mdiffCueTgtOnset, stddiffCueTgtOnset] = ...
    getDiffCueTgtOnset(onsets_cat, catEvTimes)
% computes delay between Target onset and cue Onset of each trial
% mdiffCueFOnset, stddiffCueFOnset: mean and std difference within cued trials
%created from getTgtLatencyCorr
% inf: without cue
% nan: non-completed?

%[latency_neuro, validEvents,~,~,fneuro] = getTgtNeuroLatency(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, Thresh, param, dd);
validEvents = find(~isnan(catEvTimes.tOnset) .* ~isnan(catEvTimes.fOnset)); %with or without cue

wCue = ~isinf(onsets_cat.cueOnset(validEvents));
woCue = isinf(onsets_cat.cueOnset(validEvents));
%success =  ~isnan(  dd.cOnset(validEvents));
%failure =     isnan(  dd.cOnset(validEvents)); %timeout?

diffCueTgtOnset = onsets_cat.tOnset - onsets_cat.cueOnset;

mdiffCueTgtOnset = mean(diffCueTgtOnset(validEvents(wCue)));
stddiffCueTgtOnset = std(diffCueTgtOnset(validEvents(wCue)));