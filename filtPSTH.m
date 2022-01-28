function PSTH_f = filtPSTH(PSTH, dt, sigma, causal)
%apply gaussian filter to psth (time x channels)
%dt in s
%sigma in s
%causal: 
% 0 non-causal full gaussian
% 1 causal full gaussian
% 2 causal half gaussian

if nargin < 4
    causal = 0;
end


%sigma = 1000;%ms %10: Lee...Lisburger
sigma_p = sigma/dt;%points
N = round(sigma_p*10); %points
alpha = (N-1)/(2*sigma_p);
w = gausswin(N, alpha);
if causal==2
    N = round(N/2);
    w = flipud(w(1:N));
end
if causal>0
    PSTH_c = filter(w, sum(w),cat(1, flipud(PSTH(1:N)),PSTH,flipud(PSTH(end-N:end))));
else
    PSTH_c = filtfilt(w, sum(w),cat(1, flipud(PSTH(1:N)),PSTH,flipud(PSTH(end-N:end))));
end
PSTH_f = PSTH_c(N+1:end-N-1);
