function [latency_neuro, latency_bhv, latcorr, latcorr_cue, f, fneuro] = ...
    getTgtLatencyCorr(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, param, dd, validEvents )
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

%% cue-tgt interval
diffCueFOnset = getDiffCueTgtOnset(onsets_cat, catEvTimes); 
diffCueFOnset = abs(diffCueFOnset(validEvents));

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
    theseTrials = logical(validLatency.*~isnan(latency_bhv));
    latcorr.trials = theseTrials;
    [latcorr.r, latcorr.p] = corr(latency_bhv(theseTrials), latency_neuro(theseTrials), 'type', corrOption);
catch err
    latcorr.trials = [];
    latcorr.r = nan;
    latcorr.p = nan;
end
try
    theseTrials = logical(success.*validLatency);
    latcorr.trials_success = theseTrials;
    [latcorr.r_success, latcorr.p_success] = corr(latency_bhv(theseTrials), latency_neuro(theseTrials), 'type', corrOption);
catch err
    latcorr.trials_succss = [];
    latcorr.r_success = nan;
    latcorr.p_success = nan;
end
try
    theseTrials = logical(~isnan(latency_bhv).*validLatency.*tgt_pref);
    latcorr.trials_pref = theseTrials;
    [latcorr.r_pref, latcorr.p_pref] = corr(latency_bhv(theseTrials), latency_neuro(theseTrials), 'type', corrOption);
catch err
    latcorr.trials_pref = [];
    latcorr.r_pref = nan;
    latcorr.p_pref = nan;
end
try
    theseTrials = logical(success.*validLatency.*tgt_pref);
    latcorr.trials_success_pref = theseTrials;
    [latcorr.r_success_pref, latcorr.p_success_pref] = corr(latency_bhv(theseTrials), latency_neuro(theseTrials), 'type', corrOption);
catch err
    latcorr.trials_success_pref = [];
    latcorr.r_success_pref  = nan;
    latcorr.p_success_pref = nan;
end
try
    theseTrials = logical(success.*validLatency.*(tgtDir==0).*(dd.cuedLoc(validEvents)==1));
    latcorr.trials_success_prev = theseTrials;
    [latcorr.r_success_prev, latcorr.p_success_prev] = corr(latency_bhv(theseTrials), latency_neuro(theseTrials), 'type', corrOption);
catch err
    latcorr.trials_success_prev = [];
    latcorr.r_success_prev  = nan;
    latcorr.p_success_prev = nan;
end

%% correlation between bhv latency and cue-tgt interval
try
    [latcorr_cue.r_success_prev, latcorr_cue.p_success_prev] = corr(latency_bhv(logical(success.*validLatency.*tgt_prev)), ...
        diffCueFOnset(logical(success.*validLatency.*tgt_prev)), 'type', corrOption);
catch err
    latcorr_cue.r_success_prev  = nan;
    latcorr_cue.p_success_prev = nan;
end


if nargout > 3
    f = figure('position',[ 420           5        1405         849]);

    ax(1) = subplot(2,6,[1:3 7:9]); %bhv v neuro
    plot(latency_bhv(logical(wCue.*success)), latency_neuro(logical(wCue.*success)),'ro','DisplayName','success w cue'); hold on;
    plot(latency_bhv(logical(woCue.*success)), latency_neuro(logical(woCue.*success)),'bo','DisplayName','succes wo cue'); hold on;
    plot(latency_bhv(logical(wCue.*success.*tgt_pref)), latency_neuro(logical(wCue.*success.*tgt_pref)),'o','MarkerFaceColor','r','DisplayName','success w cue prefDir'); hold on;
    plot(latency_bhv(logical(woCue.*success.*tgt_pref)), latency_neuro(logical(woCue.*success.*tgt_pref)),'o','MarkerFaceColor','b','DisplayName','succes wo cue prefDir'); hold on;
    plot(latency_bhv(logical(wCue.*fail)), latency_neuro(logical(wCue.*fail)),'rx','DisplayName','failure w cue'); hold on;
    plot(latency_bhv(logical(woCue.*fail)), latency_neuro(logical(woCue.*fail)),'bx','DisplayName','failure wo cue'); hold on;
    xlabel('behavioural latency');
    ylabel('neural latency');
    legend('location','southoutside');
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
    bhvrange = get(gca,'XLim');

    ax(2) = subplot(2,6,4); %cue v bhv
    plot(zeros(numel(find(woCue.*success)),1), latency_bhv(logical(woCue.*success)),'bo','DisplayName','success wo cue'); hold on;
    plot(zeros(numel(find(woCue.*success.*tgt_pref)),1),latency_bhv(logical(woCue.*success.*tgt_pref)),'b','MarkerFaceColor','r','DisplayName','success wo cue prefDir'); hold on;
    plot(zeros(numel(find(woCue.*fail),1)), latency_bhv(logical(woCue.*fail)),'bx','DisplayName','failure wo cue'); hold on;
    ylim(bhvrange); ylabel('bhv latency'); grid on;

    ax(3) = subplot(2,6,[5 6]); %cue v bhv
    plot(diffCueFOnset(logical(wCue.*success)), latency_bhv(logical(wCue.*success)),'ro','DisplayName','success w cue'); hold on;
    plot(diffCueFOnset(logical(wCue.*success.*tgt_pref)),latency_bhv(logical(wCue.*success.*tgt_pref)),'o','MarkerFaceColor','r','DisplayName','success w cue prefDir'); hold on;
    plot(diffCueFOnset(logical(wCue.*fail)), latency_bhv(logical(wCue.*fail)),'rx','DisplayName','failure w cue'); hold on;
    ylim(bhvrange); xlim([.48 .72]);grid on;
    xlabel('cue-tgt interval [s]');

    ax(4) = subplot(2,6,10); %cue v bhv
   plot(zeros(numel(find(woCue.*success)),1), latency_neuro(logical(woCue.*success)),'bo','DisplayName','success wo cue'); hold on;
    plot(zeros(numel(find(woCue.*success.*tgt_pref)),1),latency_neuro(logical(woCue.*success.*tgt_pref)),'b','MarkerFaceColor','r','DisplayName','success wo cue prefDir'); hold on;
    plot(zeros(numel(find(woCue.*fail),1)), latency_neuro(logical(woCue.*fail)),'bx','DisplayName','failure wo cue'); hold on;
    ylim(bhvrange); ylabel('neuro latency');grid on;

    ax(5) = subplot(2,6,[11 12]); %cue v bhv
    plot(diffCueFOnset(logical(wCue.*success)), latency_neuro(logical(wCue.*success)),'ro','DisplayName','success w cue'); hold on;
    plot(diffCueFOnset(logical(wCue.*success.*tgt_pref)), latency_neuro(logical(wCue.*success.*tgt_pref)),'o','MarkerFaceColor','r','DisplayName','success w cue prefDir'); hold on;
    plot(diffCueFOnset(logical(wCue.*fail)), latency_neuro(logical(wCue.*fail)),'rx','DisplayName','failure w cue'); hold on;
    ylim(bhvrange);xlim([.48 .72]); xlabel('cue-tgt interval [s]'); grid on;
    
    %linkaxes(ax(:),'y');
end