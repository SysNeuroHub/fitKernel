function [R] = corrPSTH_eye_ms(spk_cat, eyeData_cat, signal_r, dt_r, cutoffFreqs)

if nargin < 5
    cutoffFreqs = [1e-2 1e-1 1e0]; %[Hz] %high-pass cutoff freq for fast pupil kernel
end

visualize = 1;

%fs_eye = 1/median([eyeData_rmblk.dt]);
omitDuration = 5; %omit initial and last segments for fitting[s]


%% regress out the session-wide correlation
t_r = (eyeData_cat.t(1):dt_r:eyeData_cat.t(end))';

PSTH_r = getPSTH(spk_cat, t_r);

theseTimes = intersect(find(t_r>t_r(1)+omitDuration), find(t_r < t_r(end)-omitDuration));
PSTH_r = PSTH_r(theseTimes);
signal_r = signal_r(theseTimes);
t_r = t_r(theseTimes);

PSTH_r = PSTH_r - mode(PSTH_r);

% signal_ori = cat(1,flipud(eyeData_cat.(fdname)), ...
%     eyeData_cat.(fdname), flipud(eyeData_cat.(fdname)));
% signal_r = interp1(eyeData_cat.t, eyeData_cat.(fdname), t_r, ...
%         'linear', median(eyeData_cat.(fdname)));
ntotFrames = length(signal_r);
signal_r = signal_r - mean(signal_r);
signal_r_ext = double(cat(1, flipud(signal_r), signal_r, flipud(signal_r)));

PSTH_r_ext = cat(1, flipud(PSTH_r), PSTH_r, flipud(PSTH_r));

fs_r = 1/dt_r;
% residual = PSTH_r;

if visualize
    figure('position', [0 0 600 1000]);
    subplot(length(cutoffFreqs)+1,1,1);
    yyaxis left; plot(t_r, PSTH_r);
    yyaxis right; plot(t_r, signal_r);
    
    axis tight
end
for ift = 1:length(cutoffFreqs)
    if ift==1
        ftype = 'low';
        cutoffFreq = cutoffFreqs(1);
    else
        ftype = 'bandpass';
        cutoffFreq = cutoffFreqs(ift-1:ift);
    end
    
    order = 2;
    Wn = cutoffFreq/(fs_r/2);
    [b,a]=butter(order, Wn, ftype);
    signal_c = filtfilt(b,a,signal_r_ext);

    signal_f = signal_c(ntotFrames+1:2*ntotFrames);
    
    %residual_f = filtfilt(b,a,residual);
    PSTH_c = filtfilt(b,a,PSTH_r_ext);
    PSTH_f = PSTH_c(ntotFrames+1:2*ntotFrames);
    
    
Rc = corrcoef(PSTH_f, signal_f);
R(ift) = Rc(1,2);

    if visualize
%         subplot(length(cutoffFreqs),2,2*ift-1);
        subplot(length(cutoffFreqs)+1,1,ift+1);
     %   R2=corrcoef(PSTH_f, predicted_slow2(:,ift));
        yyaxis left; plot(t_r, PSTH_f);
        ylabel([num2str(cutoffFreq) ' [Hz]']);
        yyaxis right; plot(t_r, signal_f);
         axis tight
        title(['R(vs linear reg): ' num2str(R(ift))]);% ', R(vs kernel fit): ' num2str(R2(1,2))]);
        
        %subplot(length(cutoffFreqs),2,2*ift);
        %plot(tlags, rr_slow);
        %ylabel('kernel')
        
    end
end
end

