function [eyeData_rmblk, blinks] = removeBlinksEDF(eyeData, meta_cat, marginSize, visualize)
% [eyeData_rmblk, blinks] = removeBlinks(eyeData, marginSize, visualize)


if nargin < 3
    visualize = 0;
end

%% detect and interpolate blinks
if nargin < 3
    marginSize = 50; %frames
end

marginDur = marginSize*median(diff(eyeData.t));%seconds

% beginIdx = [];
% endIdx = [];
% for iblk = 1:length(meta_cat.STARTBLINK)
%     [~, beginIdx_c] = min(abs(eyeData.t - meta_cat.STARTBLINK(iblk)));
%     beginIdx = [beginIdx; beginIdx_c];
% 
%     [~, endIdx_c] = min(abs(eyeData.t - meta_cat.ENDBLINK(iblk)));
%     endIdx = [endIdx; endIdx_c];
% end
% blk_rising = zeros(length(eyeData.parea),1);
% blk_falling = zeros(length(eyeData.parea),1);
% blk_rising(max(1, beginIdx-marginSize))=1;
% blk_falling(min(length(eyeData.parea), endIdx+marginSize))=1;
% blkidx = find(cumsum(blk_rising)-cumsum(blk_falling) >= 1);

blinks = event2Trace(eyeData.t, [meta_cat.STARTBLINK; meta_cat.ENDBLINK], marginDur);
blkidx = find(blinks);
%startBlk = eyeData.t(blk_rising);
%endBlk = eyeData.t(blk_falling);

nblkidx = setxor(1:length(eyeData.parea), blkidx);
x_rmblinks = interp1(eyeData.t(nblkidx), eyeData.x(nblkidx), eyeData.t, 'linear',median(eyeData.x(nblkidx)));
y_rmblinks = interp1(eyeData.t(nblkidx), eyeData.y(nblkidx), eyeData.t, 'linear',median(eyeData.y(nblkidx)));
pwdth_rmblinks = interp1(eyeData.t(nblkidx), eyeData.pwdth(nblkidx), eyeData.t, 'linear', median(eyeData.pwdth(nblkidx)));
phght_rmblinks = interp1(eyeData.t(nblkidx), eyeData.phght(nblkidx), eyeData.t, 'linear',median(eyeData.phght(nblkidx)));

eyeData_rmblk = marmodata.eye(eyeData.t, x_rmblinks, y_rmblinks, pwdth_rmblinks, phght_rmblinks);


%% omit time points with outliers? ... not in this function

%% blink event time course (may include eyelink misdetection)
%0: no blink
%1: blink
% blinks = zeros(length(eyeData.t),1);
% blinks(blkidx) = 1;

% blinkOnIdx = cat(1, blkidx(1), blkidx(find(diff(blkidx)>1)+1));
% blinkOffIdx = cat(1, blkidx(find(diff(blkidx)>1)), blkidx(end));
% 
% blinkOnT = eyeData.t(blinkOnIdx);
% blinkOffT = eyeData.t(blinkOffIdx);

%% confirm if blinks are thoroughly removed
if visualize
    figure;
    plot(eyeData.t, eyeData.parea);
    hold on
    plot(eyeData.t, eyeData_rmblk.parea);
    legend('recorded','blink removed w edf data','location','southeast');
    ylabel('parea');
    xlabel('time [s]');
    axis tight;
    grid on;
end
