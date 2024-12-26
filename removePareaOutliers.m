function [eyeData_rmotl, outliers] = removePareaOutliers(eyeData, marginSize, visualize, th, diffth)
% [eyeData_rmotl, outliers] = removePareaOutliers(eyeData, marginSize, visualize, th, diffth)
% created from removeBlinks

if nargin < 5
    diffth = 5;
end

if nargin < 4
    th = 5;
end

if nargin < 3
    visualize = 0;
end

%% detect and interpolate blinks
if nargin < 2
    marginSize = 40; %frames
end


blkidx = [];
for idir = 1:2
    switch idir
        case 1
            parea_th = median(eyeData.parea(eyeData.parea>0)) - ...
                th *std(eyeData.parea(eyeData.parea>0));
            blkidx_init=find(eyeData.parea < parea_th);
        case 2
            parea_th = median(eyeData.parea(eyeData.parea>0)) + ...
                th *std(eyeData.parea(eyeData.parea>0));
            blkidx_init=find(eyeData.parea > parea_th);
    end        
    
    if ~isempty(blkidx_init)
        beginIdx = [blkidx_init(1); blkidx_init(find(diff(blkidx_init)>1)+1)];
        endIdx = [blkidx_init(find(diff(blkidx_init)>1)); blkidx_init(end)];
        blk_rising = zeros(length(eyeData.parea),1);
        blk_falling = zeros(length(eyeData.parea),1);
        blk_rising(max(1, beginIdx-marginSize))=1;
        blk_falling(min(length(eyeData.parea), endIdx+marginSize))=1;
        blkidx_c = find(cumsum(blk_rising)-cumsum(blk_falling) >= 1);
    else
        blkidx_c = [];
    end
    blkidx = cat(1, blkidx, blkidx_c);
end

nblkidx = setxor(1:length(eyeData.parea), blkidx);
x_rmblinks = interp1(eyeData.t(nblkidx), eyeData.x(nblkidx), eyeData.t, 'linear',median(eyeData.x(nblkidx)));
y_rmblinks = interp1(eyeData.t(nblkidx), eyeData.y(nblkidx), eyeData.t, 'linear',median(eyeData.y(nblkidx)));
pwdth_rmblinks = interp1(eyeData.t(nblkidx), eyeData.pwdth(nblkidx), eyeData.t, 'linear', median(eyeData.pwdth(nblkidx)));
phght_rmblinks = interp1(eyeData.t(nblkidx), eyeData.phght(nblkidx), eyeData.t, 'linear',median(eyeData.phght(nblkidx)));

eyeData_rmotl_c = marmodata.eye(eyeData.t, x_rmblinks, y_rmblinks, pwdth_rmblinks, phght_rmblinks);




%% detect blinks using findpeaks
parea_diffth = diffth*std(diff(eyeData_rmotl_c.parea));


blkidx_diff = [];
for idir = 1:2
    switch idir
        case 1
            
            [~, blkidx_diff_falling] = findpeaks([0; -diff(eyeData_rmotl_c.parea)], ...
                'MinPeakHeight', parea_diffth,'MinPeakDistance',marginSize);
            
            
            blkidx_diff_rising = [];
            for iii = 1:length(blkidx_diff_falling)
                snippet = diff([eyeData_rmotl_c.parea(blkidx_diff_falling(iii)...
                    : min(numel(eyeData_rmotl_c.parea),blkidx_diff_falling(iii)+4*marginSize))]);
                [~,idx] = max(snippet);
                if iii < length(blkidx_diff_falling)
                    blkidx_diff_rising(iii) = min(idx+blkidx_diff_falling(iii), blkidx_diff_falling(iii+1)-1);
                else
                    blkidx_diff_rising(iii) = idx+blkidx_diff_falling(iii);
                end
            end
        case 2
            [~, blkidx_diff_rising] = findpeaks(diff([0; eyeData_rmotl_c.parea]),...
                'MinPeakHeight', parea_diffth,'MinPeakDistance',marginSize);
            
            blkidx_diff_falling = [];
            for iii = 1:length(blkidx_diff_rising)
                upToThisIdx = min(blkidx_diff_rising(iii)+4*marginSize, length(eyeData_rmotl_c.parea));
                snippet = diff([eyeData_rmotl_c.parea(blkidx_diff_rising(iii):upToThisIdx)]);
                [~,idx] = max(snippet);
                if iii < length(blkidx_diff_falling)
                    blkidx_diff_falling(iii) = min(idx+blkidx_diff_rising(iii), blkidx_diff_rising(iii+1)-1);
                else
                    blkidx_diff_falling(iii) = idx+blkidx_diff_rising(iii);
                end
            end
    end
    
    if ~isempty(blkidx_diff_rising)
        blk_diff_rising = zeros(length(eyeData_rmotl_c.parea),1);
        blk_diff_falling = zeros(length(eyeData_rmotl_c.parea),1);
        blk_diff_rising(min(length(eyeData_rmotl_c.parea),blkidx_diff_rising+marginSize))=1;
        blk_diff_falling(max(1,blkidx_diff_falling-marginSize))=1;
        blkidx_diff_c = find(cumsum(blk_diff_falling)-cumsum(blk_diff_rising) == 1);
    else
        blkidx_diff_c = [];
    end
    blkidx_diff = cat(1,blkidx_diff, blkidx_diff_c);
end
blkidx_diff = sort(blkidx_diff);



nblkidx_d = setxor(1:length(eyeData_rmotl_c.parea), blkidx_diff);
x_rmblinks_d = interp1(eyeData_rmotl_c.t(nblkidx_d), eyeData_rmotl_c.x(nblkidx_d), eyeData_rmotl_c.t, 'linear',median(eyeData_rmotl_c.x(nblkidx_d)));
y_rmblinks_d = interp1(eyeData_rmotl_c.t(nblkidx_d), eyeData_rmotl_c.y(nblkidx_d), eyeData_rmotl_c.t, 'linear',median(eyeData_rmotl_c.y(nblkidx_d)));
pwdth_rmblinks_d = interp1(eyeData_rmotl_c.t(nblkidx_d), eyeData_rmotl_c.pwdth(nblkidx_d), eyeData_rmotl_c.t, 'linear', median(eyeData_rmotl_c.pwdth(nblkidx_d)));
phght_rmblinks_d = interp1(eyeData_rmotl_c.t(nblkidx_d), eyeData_rmotl_c.phght(nblkidx_d), eyeData_rmotl_c.t, 'linear',median(eyeData_rmotl_c.phght(nblkidx_d)));

eyeData_rmotl = marmodata.eye(eyeData_rmotl_c.t, x_rmblinks_d, y_rmblinks_d, pwdth_rmblinks_d, phght_rmblinks_d);


%% confirm if blinks are thoroughly removed
if visualize
    figure;
    plot(eyeData.t, eyeData.parea);
    hold on
    plot(eyeData_rmotl_c.t, eyeData_rmotl_c.parea, eyeData_rmotl.t, eyeData_rmotl.parea);
    % scatter(eyeData_rmotl.t(blkidx_margin), parea_th, 2,'g');
    legend('recorded','blink removed w threshold' ,'blink removed w differential','location','southeast');
    ylabel('parea');
    xlabel('time [s]');
    axis tight;
    grid on;
end

%% omit time points with outliers?
% cf Nishimoto 2011: "Any motion-energy signal outliers more than 3.0 standard deviations from the mean area truncated to 3.0 in order to improve stability in the model estimation procedure."
% ngth = 3; %sd
% mparea = median(eyeData_rmotl.parea);
% sdparea = std(eyeData_rmotl.parea);
% ngidx = setxor(find(eyeData_rmotl.parea > mparea+ngth*sdparea), ...
%     find(eyeData_rmotl.parea < mparea-ngth*sdparea));
% goodidx = setxor(1:length(eyeData_rmotl.parea), ngidx);
% plot(eyeData_rmotl.t(goodidx), eyeData_rmotl.parea(goodidx));


%% blink event time course (may include eyelink misdetection)
%0: no blink
%1: blink
outliers = zeros(length(eyeData_rmotl.t),1);
outliers(cat(1,blkidx, blkidx_diff)) = 1;

% blinkOnIdx = cat(1, blkidx(1), blkidx(find(diff(blkidx)>1)+1));
% blinkOffIdx = cat(1, blkidx(find(diff(blkidx)>1)), blkidx(end));
%
% blinkOnT = eyeData_rmotl.t(blinkOnIdx);
% blinkOffT = eyeData_rmotl.t(blinkOffIdx);
