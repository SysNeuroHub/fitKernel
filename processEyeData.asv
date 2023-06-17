function [eyeData_rmotl_cat, catEvTimes, t_tr, onsets_cat,meta_cat,blinks,outliers] ...
    = processEyeData(eyeData, dd, param)
%  [eyeData_rmotl_cat, catEvTimes, t_tr, onsets_cat,meta_cat,blinks,outliers] ...
%     = processEyeData(eyeData, dd, param)
% returns concatenated eye data after removing outliers
%
%created from fitPSTH_test.m on 27/1/22

[eyeData_cat, onsets_cat, meta_cat] = concatenate_eye(eyeData, dd);
t_tr={eyeData.t};
%t_cat = getTCat(t_tr);
t_cat = eyeData_cat.t;

%% detect and interpolate blinks
%removing outliers helps for larger kernels
[eyeData_rmblk_cat, blinks] = removeBlinksEDF(eyeData_cat, meta_cat, param.marginSize,1);
close;

%% remove parea outliers based on diff
[eyeData_rmotl_cat, outliers] = removePareaOutliers(eyeData_rmblk_cat, ...
    param.marginSize, param.pareaTh, param.pareaDiffTh);

%check startsacc and endsacc (to be done in concatenate_eye)
if length(meta_cat.STARTSACC) ~= length(meta_cat.ENDSACC)%41
    nSaccs = min(length(meta_cat.STARTSACC), length(meta_cat.ENDSACC));
    ngIdx=find(meta_cat.ENDSACC(1:nSaccs)-meta_cat.STARTSACC(1:nSaccs)<0);
    
    if isempty(ngIdx)
        meta_cat.STARTSACC = meta_cat.STARTSACC(1:nSaccs);
        meta_cat.ENDSACC = meta_cat.ENDSACC(1:nSaccs);
        %else
        % FILL ME?
    end
end

assert(isempty(find(meta_cat.ENDSACC-meta_cat.STARTSACC<0)));
okSacc = find(meta_cat.ENDSACC-meta_cat.STARTSACC>0);
meta_cat.STARTSACC = meta_cat.STARTSACC(okSacc);
meta_cat.ENDSACC = meta_cat.ENDSACC(okSacc);
%             tic
%             [eyeData_cat, onsets_cat, meta_cat] = concatenate_eye(eyeData, dd);%takes 90s
%             t0=toc
%             t_tr={eyeData.t};
%
%
%             %% detect and interpolate blinks
%             %removing outliers helps for larger kernels
%             tic
%             [eyeData_rmblk_cat, blinks] = removeBlinksEDF(eyeData_cat, meta_cat, param.marginSize,1);
%             close;
%             t1=toc
%
%             %% remove parea outliers based on diff
%             tic
%             [eyeData_rmotl_cat, outliers] = removePareaOutliers(eyeData_rmblk_cat, ...
%                 param.marginSize, param.pareaTh, param.pareaDiffTh);
%             screen2png(['rmOtl_' saveSuffix]);
%             close;
%             t2=toc

%% detect saccades
[startSacc, endSacc] = selectSaccades(meta_cat.STARTSACC, meta_cat.ENDSACC,...
    t_cat, (blinks+outliers>0));    % exculde blink and outlier periods

catEvTimes = onsets_cat;
[~,catEvTimes.blinkStartTimes, catEvTimes.blinkEndTimes] = trace2Event(blinks, t_cat);
[~,catEvTimes.outlierStartTimes, catEvTimes.outlierEndTimes] = trace2Event(outliers, t_cat);
catEvTimes.saccadeStartTimes = startSacc;
catEvTimes.saccadeEndTimes = endSacc;