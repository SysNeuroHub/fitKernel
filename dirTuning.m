function [mresp, sdresp, cmean, cvar, nTrials] = dirTuning(signal_tr, eyeData_tr, dd, respWin, preWin)
% mresp = dirTuning(signal_tr, eyeData, dd, respWin)
% returns mean response across trials during time window after target onset
% (dd.tOnset)
%
% mresp = dirTuning(signal_tr, eyeData, dd, respWin, preWin)
% computes response as mean activity during respWin - preWin
%
% INPUTS:
% respWin: temporal window for analysis [s]



if nargin < 5
    preWin = [];
end

dirs = unique(dd.targetloc);
mresp = zeros(length(dirs),1);
sdresp = zeros(length(dirs),1);
nTrials = zeros(length(dirs),1);

indDirRad = [];
indResp = [];
for idir = 1:length(dirs)
    
    theseTr = intersect(find(dd.successTrials), find(dd.targetloc==dirs(idir)));
    
    resp_c = [];
    for itr = 1:length(theseTr)
        thisTr = theseTr(itr);
        theseTimes  = intersect(find(eyeData_tr(thisTr).t - dd.tOnset(thisTr) > respWin(1)), ...
            find(eyeData_tr(thisTr).t - dd.tOnset(thisTr) < respWin(2)));
        resp_c(itr,1) = sum(signal_tr{thisTr}(theseTimes))/diff(respWin); %[spikes/s]

        if ~isempty(preWin)
            theseTimes  = intersect(find(eyeData_tr(thisTr).t - dd.tOnset(thisTr) > preWin(1)), ...
                find(eyeData_tr(thisTr).t - dd.tOnset(thisTr) < preWin(2)));
            resp_c(itr) = resp_c(itr) - sum(signal_tr{thisTr}(theseTimes))/diff(preWin); %[spikes/s]
        end
    end
    
    indDirRad = cat(1, indDirRad, dirs(idir)*pi/180*ones(length(theseTr),1));
    indResp = cat(1, indResp, resp_c);
    mresp(idir) = nanmean(resp_c);
    sdresp(idir) = nanstd(resp_c);
    nTrials(idir) = length(theseTr);
end

cvar = circ_var(indDirRad, indResp);
cmean = circ_mean(indDirRad, indResp);
