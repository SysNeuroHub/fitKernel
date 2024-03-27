%load('/mnt/syncitium/Daisuke/cuesaccade_data/2022/hugo/eyeCat_hugo09September_06.mat')
%load('/mnt/syncitium/Daisuke/cuesaccade_data/2022/hugo/hugo09September_06_10_linear_rReg.mat')

function [latency_neuro, validEvents, thresh_neuro, tgtDir, fig] = getTgtNeuroLatency(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, Thresh, param, dd)
% [latency_neuro, validEvents, thresh_neuro, tgtDir] = getTgtNeuroLatency(PSTH_f, t_r, onsets_cat, catEvTimes, tWin_t, Thresh, param, dd)

%% target response
%tWin_t = [-0.5 0.5]; %[-0.5 0.5]
%Thresh = 3; %4
baseline = 'tonset';
%OR  define by onset of behavioural fixation as 'fonset'; turned out to produce too short latency
tWin_f = [0  min(onsets_cat.cueOnset-onsets_cat.fOnset)];

     
%% from showTonsetByCue:

    %validEvents = intersect(find(~isnan(onset)), find(dd.cueOn==icue-1));
    %< this condition only includes all trials irrespective of the trial outcome
    validEvents = find(~isnan(catEvTimes.tOnset) .* ~isnan(catEvTimes.fOnset)); %with or without cue
       
    
    onsetTimes_t = catEvTimes.tOnset(validEvents);
    tgtDir = getTgtDir(dd.targetloc(validEvents), param.cardinalDir);
    
    [~,dirIdx]=intersect(param.cardinalDir, unique(tgtDir));
    % triggered by target onset
    [~, winSamps_t, singleResp_t] ...
        = eventLockedAvg(PSTH_f', t_r, onsetTimes_t, tgtDir, tWin_t);

    % triggered by behavioural fixation
    onsetTimes_f = catEvTimes.fOnset(validEvents);
    [~, winSamps_f, singleResp_f] ...
        = eventLockedAvg(PSTH_f', t_r, onsetTimes_f, tgtDir, tWin_f);

    %singleResp: trials x 1 x times


%% neural latency
% from     latencyV1MT/secondary/findSingleTrialLatency.m

latencyLimit = [0 tWin_t(2)];
x = 1e-2*round(1e2*winSamps_t);
sdf_t = squeeze(singleResp_t); %trial x time
sdf_f = squeeze(singleResp_f);

latency_neuro = nan(size(sdf_t,1),1);
thresh_neuro = nan(size(sdf_t,1),1);
for itr = 1:size(sdf_t,1)

    if strcmp(baseline, 'tonset')
        sdf_pre = sdf_t(itr,x<0);
    elseif strcmp(baseline, 'fonset');
        sdf_pre = sdf_f(itr,:);
    end

    sd = std(sdf_pre(:));
    mu = mean(sdf_pre(:));
    thresh = mu + Thresh*sd; %compute each trail or across all trials?

    j = find(x >= 0,1);
    if j ~= size(sdf_t,2)
        latency = x(find(sdf_t(itr,j:end) >= thresh,1) + j-1); %Calculates latency as first time point after stimulus where sdf is >= threshold
    else
        latency = NaN;
    end
    if latency <= latencyLimit(1)
        latency = NaN;
    end
    if latency >= latencyLimit(2)
        latency = NaN;
    end
    if isempty(latency)
        latency = NaN;
    end

    latency_neuro(itr) = latency;
    thresh_neuro(itr) = thresh;
end

if nargout>3
    fig = figure('position',[1003         219         588         664]);
    imagesc(x,1:size(sdf_t,1), sdf_t);
    hold on
    plot(latency_neuro,1:size(sdf_t,1),'r*');
    hline(0);
    colorbar;
end




