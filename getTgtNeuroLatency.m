%load('/mnt/syncitium/Daisuke/cuesaccade_data/2022/hugo/eyeCat_hugo09September_06.mat')
%load('/mnt/syncitium/Daisuke/cuesaccade_data/2022/hugo/hugo09September_06_10_linear_rReg.mat')

function [latency_neuro, thresh_neuro, tgtDir, fig] = ...
    getTgtNeuroLatency(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, ThreshParam, param, dd, validEvents)
% [latency_neuro, validEvents, thresh_neuro, tgtDir] = getTgtNeuroLatency(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, Thresh, param, dd)

% TODO: each row by trial type success tgt in / success tgt out / fail tgt
% in / fail tgt out / quiescent tgt in / quiescent tgt out

if isempty(ThreshParam)
    %threshOption = 'uniform'; %uniform threshold across all trials or inidividual threshold of each trial (India)
    ThreshParam.option = 'individual'; %OR  define by onset of behavioural fixation as 'fonset'; turned out to produce too short latency
    ThreshParam.dur = 0.04;%s
    ThreshParam.Thresh = 4;
end


%% target response
%baseline = 'tonset';
%tWin_f = [0  min(onsets_cat.cueOnset-onsets_cat.fOnset)];


%% from showTonsetByCue:
onsetTimes_t = catEvTimes.tOnset(validEvents);
tgtDir = getTgtDir(dd.targetloc(validEvents), param.cardinalDir);

[~,dirIdx]=intersect(param.cardinalDir, unique(tgtDir));
% triggered by target onset
[~, winSamps_t, singleResp_t] ...
    = eventLockedAvg(PSTH_f', t_r, onsetTimes_t, tgtDir, tWin_t);
singleResp_t = reshape(singleResp_t, size(singleResp_t,1), size(singleResp_t,3));
%singleResp: trials x times


%% neural latency
 [latency_neuro, thresh_neuro] = getSingleTrialLatency(singleResp_t, winSamps_t, tWin_t, ThreshParam);
% from     latencyV1MT/secondary/findSingleTrialLatency.m

%% behavioural latency for sorting trials
latency_bhv_srt = getTgtBhvLatency(onsets_cat, dd, validEvents, 1);
latency_bhv_cOn = getTgtBhvLatency(onsets_cat, dd, validEvents, 0);

%% divide trials
[choiceOutcome] = getChoiceOutcome(dd);
success = (choiceOutcome(validEvents) == 1) .* (latency_bhv_cOn - latency_bhv_srt <= 0.1);
hesitant = (latency_bhv_cOn - latency_bhv_srt > 0.1); %1st saccade did NOT lead to choice
fail = choiceOutcome(validEvents) >= 2; %quiescent + wrong

[prefDir, prefDirTrials_c] = getPrefDir(PSTH_f, t_r, onsetTimes_t, tgtDir, param);
prefDirTrials = zeros(numel(onsetTimes_t),1);
prefDirTrials(prefDirTrials_c) = 1;

[prevDir, prevDirTrials_c] = getPrevDir(tgtDir, param);
prevDirTrials = zeros(numel(onsetTimes_t),1);
prevDirTrials(prevDirTrials_c) = 1;


%% visualization
if nargout>3
    nCols = 4;
    fig = figure('position',[1003         219         nCols*588         664]);

    for istimtype = 1:nCols
        if istimtype==1
            stimtrials = prefDirTrials;
        elseif istimtype == 2
            stimtrials = 1-prefDirTrials;
        elseif istimtype==3
                stimtrials = prevDirTrials;
        elseif istimtype == 4
            stimtrials = ones(size(prefDirTrials));
        end

        %theseTrials = [find(success.*stimtrials); find(quiescent.*stimtrials); find(wrong.*stimtrials)];
        theseTrials = [];
        for ievtype = 1:3
            switch ievtype
                case 1
                    theseTrials_c = find(success.*stimtrials);
                case 2
                    theseTrials_c = find(hesitant.*stimtrials);%find(quiescent.*stimtrials);
                case 3
                    theseTrials_c = find(fail.*stimtrials);
            end
            [~, idx] = sort(latency_bhv_srt(theseTrials_c));
            theseTrials = [theseTrials; theseTrials_c(idx)];
        end
        nTrials =  [sum(success.*stimtrials) sum(hesitant.*stimtrials) sum(fail.*stimtrials)];

        axes(istimtype) = subplot(1,nCols,istimtype);
        if sum(theseTrials)==0
            continue;
        end
        imagesc(winSamps_t,1:numel(theseTrials), singleResp_t(theseTrials,:));
        hold on
        plot(latency_neuro(theseTrials), 1:numel(theseTrials),'r.');
        plot(latency_bhv_srt(theseTrials), 1:numel(theseTrials),'w.');
        vline(0);
        hline(cumsum(nTrials)+.5,gca,'-','m');
        if istimtype == 1
            title(['tgt stim on preferred direction ' num2str(prefDir) ]);
            xlabel('time from target onset [s]');
            ylabel('fail/ hesitant / success');
        elseif istimtype == 2
            title('other stimulus directions');
        elseif istimtype == 3
            title(['tgt stim on most freq direction ' num2str(prevDir) ]);
        elseif istimtype == 4
            title('all stimulus directions');
           end
    end
    linkcaxes(axes);mcolorbar(gca, .5);
end