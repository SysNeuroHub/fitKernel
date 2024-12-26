%load('/mnt/syncitium/Daisuke/cuesaccade_data/2022/hugo/eyeCat_hugo09September_06.mat')
%load('/mnt/syncitium/Daisuke/cuesaccade_data/2022/hugo/hugo09September_06_10_linear_rReg.mat')

function [latency_neuro, thresh_neuro, tgtDir, fig] = ...
    getTgtNeuroLatency(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, ThreshParam, param, dd, validEvents)
% [latency_neuro, validEvents, thresh_neuro, tgtDir] = getTgtNeuroLatency(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, Thresh, param, dd)
% create figure of single-trials sorted by behavioural latency

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

%% delay between Target onset and cue Onset of each trial
diffCueFOnset = getDiffCueTgtOnset(onsets_cat, catEvTimes); %3/6/24
diffCueFOnset = diffCueFOnset(validEvents);

%% from showTonsetByCue:
onsetTimes_t = catEvTimes.tOnset(validEvents);
tgtDir = getTgtDir(dd.targetloc(validEvents), param.cardinalDir);
% onsetTimes_t = catEvTimes.tOnset; %13/12/24
% tgtDir = getTgtDir(dd.targetloc, param.cardinalDir); %13/12/24

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
% latency_bhv_srt = getTgtBhvLatency(onsets_cat, dd, [], 1);
% latency_bhv_cOn = getTgtBhvLatency(onsets_cat, dd, [], 0);

%% divide trials
[choiceOutcome] = getChoiceOutcome(dd);
success = (choiceOutcome(validEvents) == 1) .* (latency_bhv_cOn - latency_bhv_srt <= 0.1);
fail = choiceOutcome(validEvents) >= 2; %quiescent + wrong
%success = (choiceOutcome == 1) .* (latency_bhv_cOn - latency_bhv_srt <= 0.1); %13/12/24
hesitant = (latency_bhv_cOn - latency_bhv_srt > 0.1); %1st saccade did NOT lead to choice
% fail = choiceOutcome >= 2; %quiescent + wrong %13/12/24

[prefDir, prefDirTrials_c] = getPrefDir(PSTH_f, t_r, onsetTimes_t, tgtDir, param);
prefDirTrials = zeros(numel(onsetTimes_t),1);
prefDirTrials(prefDirTrials_c) = 1;
clear prefDirTrials_c

[prevDir, prevDirTrials_c] = getPrevDir(tgtDir, param);
prevDirTrials = zeros(numel(onsetTimes_t),1);
prevDirTrials(prevDirTrials_c) = 1;
clear prevDirTrials_c


%% visualization
if nargout>3
    nCols = 1;
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
        theseTrials = []; nTrials = [];
        for ievtype = 1:3
            switch ievtype
                case 1
                    theseTrials_c = find(success.*stimtrials);
                case 2
                    theseTrials_c = find(hesitant.*stimtrials);%find(quiescent.*stimtrials);
                case 3
                    theseTrials_c = find(fail.*stimtrials);
            end

            theseTrials_c = intersect(theseTrials_c, validEvents);

            [~, idx] = sort(latency_bhv_srt(theseTrials_c));
            theseTrials = [theseTrials; theseTrials_c(idx)];
            nTrials(ievtype) = numel(theseTrials_c);


            axes(istimtype) = subplot(3,nCols,ievtype);
            if sum(theseTrials)==0
                continue;
            end
            imagesc(winSamps_t, 1:numel(theseTrials_c), singleResp_t(theseTrials_c(idx),:));
            hold on
            plot(latency_neuro(theseTrials_c(idx)), 1:numel(theseTrials_c),'r.');
            plot(latency_bhv_srt(theseTrials_c(idx)), 1:numel(theseTrials_c),'w.');
            plot(-diffCueFOnset(theseTrials_c(idx)), 1:numel(theseTrials_c),'g.');
            vline(0);
            if ievtype == 1
                title(['tgt stim on preferred direction ' num2str(prefDir) ]);
                ylabel('success');
            elseif ievtype == 2
                ylabel('hesitent');
            elseif ievtype == 3
                ylabel('failure');
                xlabel('Time from target onset [s]')
            end
        end
    end
    linkcaxes(axes);mcolorbar(gca, .5);
end