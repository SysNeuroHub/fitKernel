function [latency_neuro, latency_bhv, r, f, fneuro] = getTgtLatencyCorr(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, Thresh, param, dd)
% [latency_neuro, latency_bhv, r, f, fneuro] = getTgtLatencyCorr(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, Thresh, param, dd)
% computes neural and behavioural latencies with or without cue

[latency_neuro, validEvents,~,~,fneuro] = getTgtNeuroLatency(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, Thresh, param, dd);

%% behavioural latency
latency_bhv = onsets_cat.cOnset(validEvents) - onsets_cat.tOnset(validEvents);

wCue = ~isinf(onsets_cat.cueOnset(validEvents));
woCue = isinf(onsets_cat.cueOnset(validEvents));
success =  ~isnan(  dd.cOnset(validEvents));
failure =     isnan(  dd.cOnset(validEvents)); %timeout?

validLatency = ~isnan(latency_neuro);

%correlation beteween behaviour and neuro using only successful trials
r=corr(latency_bhv(logical(success.*validLatency)), latency_neuro(logical(success.*validLatency)));

if nargout > 3
    f = figure;
    plot(latency_bhv(logical(wCue.*success)), latency_neuro(logical(wCue.*success)),'ro','DisplayName','success w cue'); hold on;
    plot(latency_bhv(logical(woCue.*success)), latency_neuro(logical(woCue.*success)),'bo','DisplayName','succes wo cue'); hold on;
    plot(latency_bhv(logical(wCue.*failure)), latency_neuro(logical(wCue.*failure)),'rx','DisplayName','failure w cue'); hold on;
    plot(latency_bhv(logical(woCue.*failure)), latency_neuro(logical(woCue.*failure)),'bx','DisplayName','failure wo cue'); hold on;
    xlabel('behavioural latency');
    ylabel('neural latency');
    title(['r=' num2str(r)]);
    squareplot;
    axis padded; set(gca,'TickDir','out');
end