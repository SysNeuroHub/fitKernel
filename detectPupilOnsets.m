
function [dlStartTimes, dlEndTimes, csStartTimes, csEndTimes] = ...
    detectPupilOnsets(eyeData_rmotl_cat, cutoffFreq, durTh, sizeTh, excludeTimes)
%[dlStartTimes, dlEndTimes, csStartTimes, csEndTimes] = ...
%    detectPupilOnsets(eyeData_rmotl_cat, cutoffFreq, durTh, sizeTh, excludeTimes)
%returns pupil dilation and constriction starts/ends Einstein 2017 JNS

%TODO: make this function faster...

if nargin < 2 || isempty(cutoffFreq)
    cutoffFreq = 1;%[hz]
end
if nargin < 3 || isempty(durTh)
    durTh = 0.5; %[%]
end
if nargin < 4 || isempty(sizeTh)
    sizeTh = 5; %[s]
end
if nargin < 5
    excludeTimes = [];
end

%% convert to pupil diameter (prctile?)
t = eyeData_rmotl_cat.t;
ntotFrames = length(t);
[pdiam_prctile, pdiam] = getPupilDiameter(eyeData_rmotl_cat);
signal = pdiam - mean(pdiam);

%% filter pupil diamter
order = 3;
dt=median(diff(t));
fs = 1/dt;
Wn = cutoffFreq/(fs/2);
[b,a]=butter(order, Wn, 'low');

signal_c = filtfilt(b,a,cat(1,flipud(signal), ...
    signal, flipud(signal)));
signal_f = signal_c(ntotFrames+1:2*ntotFrames);

dsignal_f = [0; diff(signal_f)/dt];

% %% convert to phase
% analytic = hilbert(signal_f);
% sigphase = angle(analytic);

%% detect dilation/constriction
dilationPeriod = find(dsignal_f>0);
constrictionPeriod = find(dsignal_f<0);

dilationTrace_init = zeros(length(t),1);
dilationTrace_init(dilationPeriod) = 1;
constrictionTrace_init = zeros(length(t),1);
constrictionTrace_init(constrictionPeriod) = 1;

%% exclude events based on duration and size
[~,dlStart, dlEnd] = trace2Event(dilationTrace_init);
[~,csStart, csEnd] = trace2Event(constrictionTrace_init);

drSize = pdiam_prctile(dlEnd)-pdiam_prctile(dlStart);
drDuration = dt*(dlEnd - dlStart);
drOK = intersect(find(drSize> sizeTh), find(drDuration > durTh));

csSize = pdiam_prctile(csStart)-pdiam_prctile(csEnd);
csDuration = dt*(csEnd - csStart);
csOK = intersect(find(csSize> sizeTh), find(csDuration > durTh));

dlStartTimes_th = t(dlStart(drOK));
dlEndTimes_th = t(dlEnd(drOK));
csStartTimes_th = t(csStart(csOK));
csEndTimes_th = t(csEnd(csOK));

%% exclude events inside blinks or outliers or saccades
% excludeTimes = [catEvTimes.blinkStartTimes catEvTimes.blinkEndTimes; ...
%     catEvTimes.outlierStartTimes catEvTimes.outlierEndTimes];

drOverlapIdx = detectOverlapEvents(t, [dlStartTimes_th dlEndTimes_th], excludeTimes);
drNOverlapIdx = setdiff(1:length(dlStartTimes_th),drOverlapIdx);
dlStartTimes = dlStartTimes_th(drNOverlapIdx);
dlEndTimes = dlEndTimes_th(drNOverlapIdx);

csOverlapIdx = detectOverlapEvents(t, [csStartTimes_th csEndTimes_th], excludeTimes);
csNOverlapIdx = setdiff(1:length(csStartTimes_th),csOverlapIdx);
csStartTimes = csStartTimes_th(csNOverlapIdx);
csEndTimes = csEndTimes_th(csNOverlapIdx);


%pupilTimes = [dlStartTimes; dlEndTimes; csStartTimes; csEndTimes]

%sanity check
% dilationTrace = event2Trace(t, [dlStartTimes dlEndTimes]);
% constrictionTrace = event2Trace(t, [csStartTimes csEndTimes]);
% excludeTrace = event2Trace(t, excludeTimes);
% 
% subplot(211);
% plot(t, signal, t, signal_f);
% hold on
% subplot(212);
% plot(t,excludeTrace,'g',t, dilationTrace,'r',t, constrictionTrace, 'b');
% axis padded
% linksubaxes('x');

% trigger PSTH by pupil dilation/constriction

