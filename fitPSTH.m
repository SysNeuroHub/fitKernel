function [predicted_fr, predicted_fr_each, observed, kernelInfo] = fitPSTH(spk_cat, ...
    t_r, predictors_r, npredVars, sigma, lagRange, ridgeParam, snonlin)
% [predicted, observed, kernelInfo] = fitPSTH(spk_cat, ...
%     predictors_r, t_r, sigma, lagRange, ridgeParam)
%
% INPUTS
% t_r: timesatmps to use for fitting and output
% predictors_r: predictor sequences sampled at t_r [nVar x times]
% sigma: temporal smoothg factor for psth
% lagRange: temporal range to obtain kernels [min max] [s]
% ridgeParam: if vector, choose the best by cross-validation (time-consuming)
%
% created from denoisePSTHvis - ok I will make it simple
% preparation of predictor variables are done outside of this function

%WARNING: snonlin=1 for static nonlinearity is not correct.
%estimation of kernel is fine (almost same to the result of Yates 2017 algorighm)
%but intercept is NOT correct

nonLinOutParam = 5;%.5;% %preprocNonLinearOut
if nargin < 7
    snonlin = 0;
end
if nargin < 6
    ridgeParam = [0 1e-1 1 1e2 1e3]; %10
end

normalize = 0; %if 1, apply fitting to normalized PSTH. This is necessary to prevent the dissociation of mean activity between observed and predicted
visualize =0;
omitDuration = 0;%5; %omit initial and last segments for fitting[s]

dt_r = median(diff(t_r));

%prepare PSTH
PSTH_r = getPSTH(spk_cat, t_r);

if normalize %11/1/22
    PSTH_r = norm_std_mean(PSTH_r);
end

if isempty(sigma)
    PSTH_f = PSTH_r;
else
    PSTH_f = filtPSTH(PSTH_r, dt_r, sigma, 2);%causal
end

nanIdx = find(isnan(sum(predictors_r,1)));
PSTH_f(nanIdx) = nan;

%test 7/5/22
%PSTH_f = detrend(PSTH_f);


%% ridge regression of PSTH by behavioral signals
regIdx = intersect(find(t_r>=t_r(1)+omitDuration), find(t_r<=t_r(end)-omitDuration));
timeVec = t_r(regIdx)';
%observed = PSTH_f(regIdx) - predicted_slow(regIdx);
observed = PSTH_f(regIdx);

%static nonlinearity 
%instead of applying log to predictors as in Nishimoto 2011 (compressive nonlinearity for fMRI),
%apply log to observed signals, which is equivalent to exponenential static
%nonlinearity
if snonlin==1
    observed_forfit = log(nonLinOutParam + observed);
else
    observed_forfit = observed;
end

%predictor = log(20 + predictors_r(:,regIdx));%worse fitting
%predictor = exp(predictors_r(:,regIdx)); %worst fitting
%predictor(isinf(predictor)) = max(predictor(~isinf(predictor)));
predictor = predictors_r(:,regIdx);


%cross-validation to determine ridgeparam
if length(ridgeParam)>1
    KFolds = 5;
    mse_cv = zeros(length(ridgeParam),1);
    rr_cv = []; mexpval_cv = []; r0_cv=[];
    for irp = 1:length(ridgeParam)
        [mse_c, rr_cv(:,:,:,irp), r0_cv(irp), expval_c] = ridgeXs_cv(KFolds, timeVec, ...
            predictor, observed_forfit, lagRange, ridgeParam(irp));
        mse_cv(irp) = mean(mse_c);
        mexpval_cv(irp) = mean(expval_c);
    end
    
    [~,thisRp] = min(mse_cv);
    ridgeParam = ridgeParam(thisRp);
    kernel_cv = rr_cv(:,:,:,thisRp);
    
    kernelInfo.intercept_cv = r0_cv;
    kernelInfo.kernel_cv = kernel_cv;
    kernelInfo.mse_cv = mse_cv(thisRp);
    kernelInfo.expval_cv = mexpval_cv(thisRp);
end

[kernel, r0, predicted] = ridgeXs(timeVec, predictor, observed_forfit, ...
    lagRange, ridgeParam);

% mse = mean((observed_forfit - predicted').^2);
% expval = 100*(1 - mse / mean((observed_forfit - mean(observed_forfit)).^2));
% R = corrcoef(observed_forfit, predicted');


%% decompose model response to each filter
predicted_each = zeros(numel(npredVars), length(t_r));
for ivar = 1:numel(npredVars)
    if ivar==1
        theseVarIdx = 1:npredVars(1);
    else
        theseVarIdx = sum(npredVars(1:ivar-1))+1:sum(npredVars(1:ivar));
    end
    if size(lagRange,1)== 1
        thisLagRange = lagRange;
    else
        thisLagRange = [min(lagRange(:,1)) max(lagRange(:,2))];
    end
    
    predicted_each(ivar, :) = predictXs(t_r, predictors_r(theseVarIdx,:), ...
        r0, kernel(:,theseVarIdx), thisLagRange);
end
        
        
if snonlin
    predicted_fr = exp(predicted);%  - nonLinOutParam;%/ exp(nonLinOutParam*sum(kernel));
    predicted_fr_each = exp(predicted_each);% - nonLinOutParam;
    kernelInfo.kernel = kernel;
    %kernelInfo.kernel = exp(kernel) - nonLinOutParam; %NG
else
    predicted_fr = predicted;
    predicted_fr_each = predicted;
    kernelInfo.kernel = kernel;
end

mse = mean((observed - predicted_fr').^2);
expval = 100*(1 - mse / mean((observed - mean(observed)).^2));
R = corrcoef(observed, predicted_fr');

fs = 1/median(diff(single(timeVec)));
%lags = round(lagRange(1)*fs):round(lagRange(2)*fs);
lags = min(round(lagRange(:,1)*fs)):max(round(lagRange(:,2)*fs));
tlags = lags/fs;

% R=corrcoef(cat(2,observed, predicted'));

kernelInfo.intercept = r0;
kernelInfo.fs = fs;
kernelInfo.tlags = tlags;
kernelInfo.ridgeParam = ridgeParam;
kernelInfo.mse = mse;
kernelInfo.expval = expval;
kernelInfo.corrcoef = R(1,2);

if visualize
%     figure;
%     %% figure kernels
%     plot(tlags, squeeze(kernel_cv), 'color', [.5 .5 .5]);
%     hold on;
%     plot(tlags, kernel, 'k','linewidth',2);
%     grid on;
%     title(['corr coef observed vs predicted by position(x) ' num2str(R(1,2))]);
%     xlabel('delay [s]');
%     grid on;
%     marginplots;
    
    %screen2png(['Kernels_filterSigma' num2str(sigma) 'ms']);
    
    
    %% figure traces over time
    figure('position',[0 0 1900 1400]);
    ax3(1)=subplot(211); plot(t_r, predictors_r);
    ylabel('signal');grid on;
    ax3(2)=subplot(212); plot(timeVec, observed, timeVec, predicted_fr);
    grid on;
    legend('observed PSTH (filtered)', ['fitted']);
    ylabel('psth');xlabel('time [s]');
    axis tight;
    linkaxes(ax3(:), 'x');
    %screen2png(['tcourse_resampleSigma' num2str(sigma) 'ms']);
end


end

