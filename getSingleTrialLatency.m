function [latency_neuro, thresh_neuro] = getSingleTrialLatency(singleResp_t, winSamps_t, tWin_t, threshOption, Thresh)
%[latency_neuro, thresh_neuro] = getSingleTrialLatency(singleResp_t, winSamps_t, tWin_t, threshOption)
% from     latencyV1MT/secondary/findSingleTrialLatency.m
% singleResp_t: trial x time

baseline = 'tonset'; %used to have an option of fonset but turned out to be noisier

if nargin < 5
    Thresh = 4; %s.d.
end
if nargin < 4
    threshOption = 'individual'; %india's choice
end

latencyLimit = [0 tWin_t(2)];
x = 1e-2*round(1e2*winSamps_t);
sdf_t = singleResp_t; %trial x time

%% define threshold
thresh_neuro = nan(size(sdf_t,1),1);
if strcmp(threshOption, 'individual')
    for itr = 1:size(sdf_t,1)

        if strcmp(baseline, 'tonset')
            sdf_pre = sdf_t(itr,x<0);
        end

        sd = std(sdf_pre(:));
        mu = mean(sdf_pre(:));
        thresh = mu + Thresh*sd; %compute each trail or across all trials?
        thresh_neuro(itr) = thresh;
    end
elseif strcmp(threshOption, 'uniform');
    if strcmp(baseline, 'tonset')
        sdf_pre = sdf_t(:,x<0);
    end

    sd = std(sdf_pre(:));
    mu = mean(sdf_pre(:));
    thresh = mu + Thresh*sd; %compute each trail or across all trials?
    thresh_neuro(:) = thresh;
end


%% detect latency
latency_neuro = nan(size(sdf_t,1),1);
for itr = 1:size(sdf_t,1)
    j = find(x >= 0,1);
    if j ~= size(sdf_t,2)
        latency = x(find(sdf_t(itr,j:end) >= thresh_neuro(itr),1) + j-1); %Calculates latency as first time point after stimulus where sdf is >= threshold
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
end