function [latency_neuro, thresh_neuro] = getSingleTrialLatency(singleResp_t, winSamps_t, ...
    tWin_t, threshParam)
%[latency_neuro, thresh_neuro] = getSingleTrialLatency(singleResp_t, winSamps_t, tWin_t, threshOption)
% from     latencyV1MT/secondary/findSingleTrialLatency.m
% singleResp_t: trial x time

baseline = 'tonset'; %used to have an option of fonset but turned out to be noisier

if nargin < 4
    threshParam.dur = 0;
    threshParam.thresh = 4; %s.d.
    threshParam.option = 'individual'; %india's choice
end

assert(strcmp(threshParam.option, 'individual') + strcmp(threshParam.option, 'uniform') )

latencyLimit = [0 tWin_t(2)];
x = 1e-2*round(1e2*winSamps_t);
sdf_t = singleResp_t; %trial x time

dt = median(diff(x));
tidx_dur = round(threshParam.dur/dt)-1;

%% define threshold
thresh_neuro = nan(size(sdf_t,1),1);
if strcmp(threshParam.option, 'individual')
    for itr = 1:size(sdf_t,1)

        if strcmp(baseline, 'tonset')
            sdf_pre = sdf_t(itr,x<0);
        end

        sd = std(sdf_pre(:));
        mu = mean(sdf_pre(:));
        thresh = mu + threshParam.thresh*sd; %compute each trail or across all trials?
        thresh_neuro(itr) = thresh;
    end
elseif strcmp(threshParam.option, 'uniform');
    if strcmp(baseline, 'tonset')
        sdf_pre = sdf_t(:,x<0);
    end

    sd = std(sdf_pre(:));
    mu = mean(sdf_pre(:));
    thresh = mu + threshParam.thresh*sd; %compute each trail or across all trials?
    thresh_neuro(:) = thresh;
end


%% detect latency
latency_neuro = nan(size(sdf_t,1),1);
for itr = 1:size(sdf_t,1)
    j = find(x >= 0,1);
    tidx_a = find(sdf_t(itr,j:end) >= thresh_neuro(itr)) + j-1;

    if j ~= size(sdf_t,2) && ~isempty(tidx_a)
        %Calculates latency as first time point after stimulus where sdf is >= threshold
        
        for xxx = 1:numel(tidx_a)
            tidx = tidx_a(xxx);

            %minimal duration surpassing the threshold
            if isempty(find(tidx+tidx_dur<=numel(x), 1))
                latency = NaN; continue;
            elseif sum(sdf_t(itr, tidx:tidx+tidx_dur) >= thresh_neuro(itr)) == tidx_dur+1
                latency = x(tidx); break;
            else 
                latency = NaN;
            end
        end
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