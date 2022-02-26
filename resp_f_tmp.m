function [avgConsetResp_f, avgConsetPhase] = resp_f_tmp(singleConsetResp, sortedConsetLabels, uniqueConsetLabels, cutoffFreq)

fs_r = 50;%
%cutoffFreq = [7 13];
ftype = 'bandpass';
order = 2;
Wn = cutoffFreq/(fs_r/2);
[b,a] = butter(order, Wn, ftype);

singleConsetResp_f = [];
for itr = 1:size(singleConsetResp,1)
    for imodality = 1:size(singleConsetResp,2)
        singleConsetResp_f(itr,imodality,:) = filtfilt(b,a,double(squeeze(singleConsetResp(itr,imodality,:))));
        singleConsetPhase(itr,imodality,:) =  angle(hilbert(singleConsetResp_f(itr,imodality,:)));       
    end
end
avgConsetResp_f = [];
for idir = 1:length(uniqueConsetLabels)
    theseTr = find(sortedConsetLabels==uniqueConsetLabels(idir));
    avgConsetResp_f(idir,:,:) = mean(singleConsetResp_f(theseTr,:,:));
    avgConsetPhase(idir,:,:) = circ_mean(singleConsetPhase(theseTr,:,:),[],1);
end