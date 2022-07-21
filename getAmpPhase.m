function [filtered, amp, phase] = getAmpPhase(trace, fs, cutoffFreq)

ftype = 'bandpass';
order = 2;
Wn = cutoffFreq/(fs/2);
[b,a] = butter(order, Wn, ftype);

filtered = filtfilt(b,a,double(trace));
analytic = hilbert(filtered);
amp = abs(analytic);
phase =  angle(analytic);
