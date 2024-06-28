function [latency_neuro, latency_bhv, latcorr, f, fneuro] = getTgtLatencyCorr(PSTH_f, t_r, onsets_cat, ...
    catEvTimes, tWin_t, param, dd, validEvents )
% [latency_neuro, latency_bhv, r, f, fneuro] = getTgtLatencyCorr(PSTH_f, t_r, onsets_cat, ...
% catEvTimes, tWin_t, Thresh, param, dd)
% computes neural and behavioural latencies with or without cue

corrOption = 'Spearman'; %'Pearson'

%% neural latency
[latency_neuro, ~,~,fneuro] = getTgtNeuroLatency(PSTH_f,  t_r,  onsets_cat,  ...
    catEvTimes, tWin_t, param.threshParam, param, dd, validEvents);

%% behavioural latency
option = 1;
latency_bhv = getTgtBhvLatency(onsets_cat, dd, validEvents, option);

wCue = ~isinf(onsets_cat.cueOnset(validEvents));
woCue = isinf(onsets_cat.cueOnset(validEvents));

[choiceOutcome] = getChoiceOutcome(dd);
hesitant = getChoiceOutcome_hesitant(onsets_cat, dd, validEvents);

success = (choiceOutcome(validEvents) == 1) .* ~hesitant;
fail = choiceOutcome(validEvents) >= 2;


theseEvents = validEvents;%(success);%better to use all conditions??
onsetTimes = catEvTimes.tOnset; 
onsetTimes = onsetTimes(theseEvents);
tgtDir = getTgtDir(dd.targetloc, param.cardinalDir);
tgtDir = tgtDir(theseEvents);

%% select trials with preferred direction (showTonsetResp)
[prefDir, prefDirTrials] = getPrefDir(PSTH_f, t_r, onsetTimes, tgtDir, param);
tgt_pref = zeros(numel(validEvents),1);
tgt_pref(prefDirTrials) = 1;

%% select trials with most frequent direction
[prevDir, prevDirTrials] = getPrevDir(tgtDir, param);
tgt_prev = zeros(numel(validEvents),1);
tgt_prev(prevDirTrials) = 1;

validLatency = ~isnan(latency_neuro);

%correlation beteween behaviour and neuro using only successful trials
try
    [latcorr.r, latcorr.p] = corr(latency_bhv(logical(validLatency.*~isnan(latency_bhv))), ...
        latency_neuro(logical(validLatency.*~isnan(latency_bhv))), 'type', corrOption);
catch err
    latcorr.r = nan;
    latcorr.p = nan;
end
try
    [latcorr.r_success, latcorr.p_success] = corr(latency_bhv(logical(success.*validLatency)), ...
        latency_neuro(logical(success.*validLatency)), 'type', corrOption);
catch err
    latcorr.r_success = nan;
    latcorr.p_success = nan;
end
try
    [latcorr.r_pref, latcorr.p_pref] = corr(latency_bhv(logical(validLatency.*~isnan(latency_bhv).*tgt_pref)), ...
        latency_neuro(logical(validLatency.*~isnan(latency_bhv).*tgt_pref)), 'type', corrOption);
catch err
    latcorr.r_pref = nan;
    latcorr.p_pref = nan;
end
try
    [latcorr.r_success_pref, latcorr.p_success_pref] = corr(latency_bhv(logical(success.*validLatency.*tgt_pref)), ...
        latency_neuro(logical(success.*validLatency.*tgt_pref)), 'type', corrOption);
catch err
    latcorr.r_success_pref  = nan;
    latcorr.p_success_pref = nan;
end
try
    [latcorr.r_success_prev, latcorr.p_success_prev] = corr(latency_bhv(logical(success.*validLatency.*tgt_prev)), ...
        latency_neuro(logical(success.*validLatency.*tgt_prev)), 'type', corrOption);
catch err
    latcorr.r_success_prev  = nan;
    latcorr.p_success_prev = nan;
end

if nargout > 3
    f = figure;
    plot(latency_bhv(logical(wCue.*success)), latency_neuro(logical(wCue.*success)),'ro','DisplayName','success w cue'); hold on;
    plot(latency_bhv(logical(woCue.*success)), latency_neuro(logical(woCue.*success)),'bo','DisplayName','succes wo cue'); hold on;
    plot(latency_bhv(logical(wCue.*success.*tgt_pref)), latency_neuro(logical(wCue.*success.*tgt_pref)),'o','MarkerFaceColor','r','DisplayName','success w cue prefDir'); hold on;
    plot(latency_bhv(logical(woCue.*success.*tgt_pref)), latency_neuro(logical(woCue.*success.*tgt_pref)),'o','MarkerFaceColor','b','DisplayName','succes wo cue prefDir'); hold on;
    plot(latency_bhv(logical(wCue.*fail)), latency_neuro(logical(wCue.*fail)),'rx','DisplayName','failure w cue'); hold on;
    plot(latency_bhv(logical(woCue.*fail)), latency_neuro(logical(woCue.*fail)),'bx','DisplayName','failure wo cue'); hold on;
    xlabel('behavioural latency');
    ylabel('neural latency');
    legend('location','eastoutside');
    tname = sprintf(['all dir: %.2f (p=%.2f), \n' ...
        'success all dir: %.2f (p=%.2f) \n' ...
        'success pref dir: %.2f (p=%.2f)\n' ...
         'success prev dir: %.2f (p=%.2f)'], ...
        latcorr.r, latcorr.p, ...
        latcorr.r_success, latcorr.p_success, ...
        latcorr.r_success_pref, latcorr.p_success_pref,...
        latcorr.r_success_prev, latcorr.p_success_prev);
    title(tname); %title(['r=' num2str(r) ', r(success)=' num2str(r_success)]);
    squareplot;
    axis padded; set(gca,'TickDir','out');
end