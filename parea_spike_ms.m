function [R] = parea_spike_ms(PSTH_r, pdiam_r, t_r, cutoffFreqs)

if nargin < 4
    cutoffFreqs = [1e-2 1e-1 1e0]; %[Hz] %high-pass cutoff freq for fast pupil kernel
end

visualize = 1;




%trim begin/end
% omitDuration = 5; %omit initial and last segments for fitting[s]
% tidx = intersect(find(t_r>t_r(1)+omitDuration), find(t_r<t_r(end)-omitDuration));
% PSTH_r = PSTH_r(tidx);
% pdiam_r = pdiam_r(tidx);
% t_r = t_r(tidx);

ntotFrames = length(pdiam_r);
signal_r_ext = cat(1, flipud(pdiam_r), pdiam_r, flipud(pdiam_r));
%signal_r_ext = cat(1, zeros(ntotFrames,1), pdiam_r, zeros(ntotFrames,1));

PSTH_r_ext = cat(1, flipud(PSTH_r), PSTH_r, flipud(PSTH_r));

dt_r = median(diff(t_r));
fs_r = 1/dt_r;

if visualize
    figure('position', [0 0 1600 1000]);
    subplot(length(cutoffFreqs)+1,1,1);
    yyaxis left
    plot(t_r, PSTH_r);%,'Color',[0.5 0.5 1]);hold on
    %plot(t_r, sgolayfilt(PSTH_r, 3, 1+1e3), 'b','linewidth',2);
    
    yyaxis right
    plot(t_r, pdiam_r);%, 'Color', [1 .5 .5]);
    %plot(t_r, sgolayfilt(pdiam_r, 3, 1+1e3), 'r','linewidth',2);
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
    [b,a] = butter(order, Wn, ftype);
    
    signal_c = filtfilt(b,a,double(signal_r_ext));%NG!!
    signal_f = signal_c(ntotFrames+1:2*ntotFrames);
    
    PSTH_c = filtfilt(b,a,double(PSTH_r_ext));
    PSTH_f = PSTH_c(ntotFrames+1:2*ntotFrames);
    
    
    Rc=corrcoef(PSTH_f, signal_f);
    R(ift) = Rc(1,2);
    
    if visualize
        subplot(length(cutoffFreqs)+1,1,ift+1);
        yyaxis left
        plot(t_r, PSTH_f);
        
        yyaxis right
        plot(t_r, signal_f);
        
        axis tight
        title(['R(vs linear reg): ' num2str(R(ift))]);% ', R(vs kernel fit): ' num2str(R2(1,2))]);
        ylabel([num2str(cutoffFreq) ' [Hz]']);
        
    end
end



%
% if visualize
%     figure;
%     PSTH_regout_slow = PSTH_r-sum(predicted_slow,2);
%     plot(signal_r, PSTH_r, signal_r, sum(PSTH_regout_slow,2));
%     hold on
%     axis square;
%     legend('before removal of slow trend','after removal of slow trend');
%     xlabel('signal low-pass filtered');
%     ylabel('psth');
%     R=corrcoef(PSTH_regout_slow, signal_r);
% end
