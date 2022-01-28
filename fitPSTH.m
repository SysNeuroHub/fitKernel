function [predicted, observed, kernelInfo] = fitPSTH(spk_cat, ...
    t_r, predictors_r, sigma, lagRange, ridgeParam)
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

if nargin < 6
    ridgeParam = [0 1e-1 1 1e2 1e3]; %10
end

normalize = 0; %if 1, apply fitting to normalized PSTH. This is necessary to prevent the dissociation of mean activity between observed and predicted
visualize = 0;
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

%% ridge regression of PSTH by behavioral signals
regIdx = intersect(find(t_r>=t_r(1)+omitDuration), find(t_r<=t_r(end)-omitDuration));
timeVec = t_r(regIdx)';
%observed = PSTH_f(regIdx) - predicted_slow(regIdx);
observed = PSTH_f(regIdx);

predictor = predictors_r(:,regIdx);



%cross-validation to determine ridgeparam
if length(ridgeParam)>1
    KFolds = 5;
    mse_cv = zeros(length(ridgeParam),1);
    rr_cv = []; mexpval_cv = []; r0_cv=[];
    for irp = 1:length(ridgeParam)
        [mse_c, rr_cv(:,:,:,irp), r0_cv(irp), expval_c] = ridgeXs_cv(KFolds, timeVec, ...
            predictor, observed, lagRange, ridgeParam(irp));
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

[kernel, r0, predicted] = ridgeXs(timeVec, predictor, observed, ...
    lagRange, ridgeParam);

mse = mean((observed - predicted').^2);
expval = 100*(1 - mse / mean((observed - mean(observed)).^2));
R = corrcoef(observed, predicted');

    
fs = 1/median(diff(single(timeVec)));
lags = round(lagRange(1)*fs):round(lagRange(2)*fs);
tlags = lags/fs;

% R=corrcoef(cat(2,observed, predicted'));

kernelInfo.kernel = kernel;
kernelInfo.intercept = r0;
kernelInfo.fs = fs;
kernelInfo.tlags = tlags;
kernelInfo.ridgeParam = ridgeParam;
kernelInfo.mse = mse;
kernelInfo.expval = expval;
kernelInfo.corrcoef = R(1,2);

if visualize
    figure;
    %% figure kernels
    plot(tlags, squeeze(kernel_cv), 'color', [.5 .5 .5]);
    hold on;
    plot(tlags, kernel, 'k','linewidth',2);
    grid on;
    title(['corr coef observed vs predicted by position(x) ' num2str(R(1,2))]);
    xlabel('delay [s]');
    grid on;
    marginplots;
    
    %screen2png(['Kernels_filterSigma' num2str(sigma) 'ms']);
    
    
    %% figure traces over time
    figure('position',[0 0 1900 1400]);
    ax3(1)=subplot(211); plot(t_r, signal_hp);
    ylabel('signal');grid on;
    ax3(2)=subplot(212); plot(timeVec, observed, timeVec, predicted);
    grid on;
    legend('observed PSTH (filtered)', ['fitted with' fdname]);
    ylabel('psth');xlabel('time [s]');
    axis tight;
    linkaxes(ax3(:), 'x');
    %screen2png(['tcourse_resampleSigma' num2str(sigma) 'ms']);
end


end

