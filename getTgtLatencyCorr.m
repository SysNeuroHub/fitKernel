function [latency_neuro, latency_bhv, latcorr, f, fneuro] = getTgtLatencyCorr(PSTH_f, t_r, onsets_cat, ...
    catEvTimes, tWin_t, Thresh, param, dd )
% [latency_neuro, latency_bhv, r, f, fneuro] = getTgtLatencyCorr(PSTH_f, t_r, onsets_cat, ...
% catEvTimes, tWin_t, Thresh, param, dd)
% computes neural and behavioural latencies with or without cue

%validEvents = intersect(find(~isnan(onset)), find(dd.cueOn==icue-1));
%< this condition only includes all trials irrespective of the trial outcome
validEvents = find(~isnan(catEvTimes.tOnset) .* ~isnan(catEvTimes.fOnset)); %with or without cue

%% neural latency
[latency_neuro, ~,~,fneuro] = getTgtNeuroLatency(PSTH_f,  t_r,  onsets_cat,  ...
    catEvTimes, tWin_t, Thresh, param, dd, validEvents);

%% behavioural latency
option = 1;
latency_bhv = getTgtBhvLatency(onsets_cat, dd, validEvents, option);

wCue = ~isinf(onsets_cat.cueOnset(validEvents));
woCue = isinf(onsets_cat.cueOnset(validEvents));

[choiceOutcome] = getChoiceOutcome(dd);
hesitant = getChoiceOutcome_hesitant(onsets_cat, dd, validEvents);

success = (choiceOutcome(validEvents) == 1) .* ~hesitant;
fail = choiceOutcome(validEvents) >= 2;


%% select trials with preferred direction (showTonsetResp)
theseEvents = validEvents;%(success);%better to use all conditions??
onsetTimes = catEvTimes.tOnset; 
onsetTimes = onsetTimes(theseEvents);
tgtDir = getTgtDir(dd.targetloc, param.cardinalDir);
tgtDir = tgtDir(theseEvents);
[prefDir, prefDirTrials] = getPrefDir(PSTH_f, t_r, onsetTimes, tgtDir, param);
tgt_in = zeros(numel(validEvents),1);
tgt_in(prefDirTrials) = 1;


validLatency = ~isnan(latency_neuro);

%correlation beteween behaviour and neuro using only successful trials
latcorr.r = corr(latency_bhv(logical(validLatency.*~isnan(latency_bhv))), latency_neuro(logical(validLatency.*~isnan(latency_bhv))));
latcorr.r_success = corr(latency_bhv(logical(success.*validLatency)), latency_neuro(logical(success.*validLatency)));
latcorr.r_pref = corr(latency_bhv(logical(validLatency.*~isnan(latency_bhv).*tgt_in)), latency_neuro(logical(validLatency.*~isnan(latency_bhv).*tgt_in)));
latcorr.r_success_pref = corr(latency_bhv(logical(success.*validLatency.*tgt_in)), latency_neuro(logical(success.*validLatency.*tgt_in)));

if nargout > 3
    f = figure;
    plot(latency_bhv(logical(wCue.*success)), latency_neuro(logical(wCue.*success)),'ro','DisplayName','success w cue'); hold on;
    plot(latency_bhv(logical(woCue.*success)), latency_neuro(logical(woCue.*success)),'bo','DisplayName','succes wo cue'); hold on;
    plot(latency_bhv(logical(wCue.*success.*tgt_in)), latency_neuro(logical(wCue.*success.*tgt_in)),'o','MarkerFaceColor','r','DisplayName','success w cue prefDir'); hold on;
    plot(latency_bhv(logical(woCue.*success.*tgt_in)), latency_neuro(logical(woCue.*success.*tgt_in)),'o','MarkerFaceColor','b','DisplayName','succes wo cue prefDir'); hold on;
    plot(latency_bhv(logical(wCue.*fail)), latency_neuro(logical(wCue.*fail)),'rx','DisplayName','failure w cue'); hold on;
    plot(latency_bhv(logical(woCue.*fail)), latency_neuro(logical(woCue.*fail)),'bx','DisplayName','failure wo cue'); hold on;
    xlabel('behavioural latency');
    ylabel('neural latency');
    legend('location','eastoutside');
    tname = sprintf('all dir: %.2f, success all dir: %.2f \n success pref dir: %.2f', latcorr.r, latcorr.r_success,  latcorr.r_success_pref);
    title(tname); %title(['r=' num2str(r) ', r(success)=' num2str(r_success)]);
    squareplot;
    axis padded; set(gca,'TickDir','out');
end