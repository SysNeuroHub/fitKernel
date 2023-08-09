function [predicted, predicted_each, PSTH_f, kernelInfo] = fitPSTH_cv(spk_cat, ...
    t_r, predictorNames, predictors_r, npredVars, sigma, kernelInterval, ...
    lagRange, ridgeParam, trIdx_r, option)
%[predicted, predicted_each, PSTH_f, kernelInfo] = fitPSTH_cv(spk_cat, ...
%    t_r, predictorNames, predictors_r, npredVars, sigma, kernelInterval, ...
%    lagRange, ridgeParam, trIdx_r)
%created from https://github.com/pillowlab/neuroGLM/blob/master/docs/tutorial.md
%
% currently using regression.xvalidationIdx from obsolte neuroGLM by jyts.
% the path for the obsolte neuroGLM should be below the new neuroGLM
% Future version should stop using this function
%
% fmincon occupies most of the time, GPU computation is not yet implemented
%cf. https://au.mathworks.com/matlabcentral/answers/1605350-fmincon-running-on-gpu
%
% predicted: predicted trace by all kernels
% predicted_each: predicted traces by each kernel
% sum(predicted_each) ~= predicted because static nonlinearity was applied differently
%
% 13/7/23 added option=5 (rReg)

sparse = 1; %whether to use sparse matrix in compileSparseDesignMatrix
useGPU = 0;

if nargin < 11
    option = 4;
end

useSptrain = 0;

unitOfTime = 's';
uniqueID = 1;
KFolds = 5;
kernelInfo.basisType = 'raised cosine';
detrend = 1; %22/7/22
%nBasisFunctions = 20;
%offset = -25; %slide the kernel window back in time


if size(lagRange,1) == 1
    lagRange = repmat(lagRange, [numel(npredVars) 1]);
end

dt_r = median(diff(t_r));

PSTH_r = getPSTH(spk_cat, t_r);
PSTH_f = filtPSTH(PSTH_r, dt_r, sigma, 2, detrend);


expt = buildGLM.initExperiment(unitOfTime, dt_r, uniqueID);

%Registering variables to the experiment
for ivar = 1:numel(npredVars)
    expt = buildGLM.registerContinuous(expt, predictorNames{ivar}, [], npredVars(ivar));
end
if useSptrain
    expt = buildGLM.registerSpikeTrain(expt, 'sptrain', 'Our Neuron'); % Spike train!!!
else
    expt = buildGLM.registerContinuous(expt, 'spfilt', 'filtered PSTH', 1);
end

%Loading the data for each trial
K = t_r(1);
t_r = t_r - K;
spk_cat = spk_cat - K;
for itr = 1:numel(trIdx_r)
    t_cat_tr = t_r(trIdx_r{itr});
    duration = numel(trIdx_r{itr})*dt_r;
    
    trial = buildGLM.newTrial(expt, duration);
    
    if useSptrain
        trial.sptrain = spk_cat(spk_cat>=t_cat_tr(1) & spk_cat<=t_cat_tr(end)) ...
            - t_cat_tr(1); %must start from 0
    else
        trial.spfilt = dt_r*PSTH_f(trIdx_r{itr});
    end
    
    for ivar = 1:numel(npredVars)
        if ivar==1
            idx = 1:npredVars(1);
        else
            idx = sum(npredVars(1:ivar-1))+1:sum(npredVars(1:ivar));
        end
        trial.(predictorNames{ivar}) = predictors_r(idx,trIdx_r{itr})';
    end
    %add the trial object to the experiment object with an associated trial index kTrial:
    expt = buildGLM.addTrial(expt, trial, itr);
end

% creating a design specification object.
dspec = buildGLM.initDesignSpec(expt);


nBasisFunctions = zeros(numel(npredVars),1);
for ivar = 1:numel(npredVars)
    covLabel = predictorNames{ivar};
    stimHandle = basisFactory.rawStim(covLabel);
    
    kernelDur = diff(lagRange(ivar,:));
    nBasisFunctions(ivar) = ceil(kernelDur/kernelInterval);
    bs=basisFactory.makeSmoothTemporalBasis(kernelInfo.basisType, kernelDur, ...
        nBasisFunctions(ivar), expt.binfun);
    
    %offset = lagRange(ivar,1)/kernelInterval;
    offset = ceil(lagRange(ivar,1)/dt_r);
    
    dspec = buildGLM.addCovariate(dspec, covLabel,covLabel, stimHandle, ...
        bs, offset);
end



xvFolds = regression.xvalidationIdx(numel(trIdx_r), KFolds, false, true);

w = [];
predicted = zeros(numel(t_r),1);
predicted_each = zeros(numel(t_r), numel(npredVars));
expval = zeros(1,KFolds);
mse = zeros(1,KFolds);
R = zeros(1,KFolds);
for ifold = 1:KFolds
    disp(['fitPSTH_cv:' num2str(ifold) '/' num2str(KFolds)]);
    dm = buildGLM.compileSparseDesignMatrix(dspec, xvFolds{ifold,1});
    
    %If your design matrix is not very sparse (less than 10% sparse, for example),
    %it's better to conver the design matrix to a full (dense) matrix for speed.
    if sparse
        dm.X = full(dm.X);
    end

    % rho = 1; % ridge parameter (play around with this to see how it effects fitting. Use cross validation on the training set to set it properly)
    % dm.addBiasColumn('right'); % augment with column of ones
    dm = buildGLM.addBiasColumn(dm);
    
    %% Get the spike trains back to regress against
    % identical to dt_r*PSTH_f
    if useSptrain
        y = buildGLM.getBinnedSpikeTrain(expt, 'sptrain', dm.trialIndices);
    else
        y = buildGLM.getResponseVariable(expt, 'spfilt', dm.trialIndices);
    end
    
    %% Doing the actual regression
    switch option
        case 1
            % option1: simple least squares
            dmXg = gpuArray([dm.X]);
            yg = gpuArray(y);
            w(:,ifold) = dmXg' * dmXg \ dmXg' * yg;
            
            %w(:,ifold) = dm.X' * dm.X \ dm.X' * y;
            
        case 2
            % option2: Maximum likelihood estimation using glmfit
            [w(:,ifold), dev, stats] = glmfit(dm.X, y, 'poisson', 'link', 'log');
            
        case 3
            % option3: ridge regression with static nonlinearity
            dspec.model.regressionMode='RIDGEFIXED';
            dspec.model.nlfun = @expfun;
            ndx = []; %glmspike.m
            M = regression.doRegressionPoisson(dm.X, y, dspec, ndx, dt_r, ridgeParam); %from Yate's classy neuroGLM: SLOW
            w(:,ifold) = M.khat;
            
        case 4
            %option4: from neuroGLM/tutorial.m
            wInit = dm.X \ y;
            
            %% Use matRegress for Poisson regression
            % it requires `fminunc` from MATLAB's optimization toolbox
            %addpath('C:\Users\dshi0006\git\neuroGLM\matRegress')
            
            fnlin = @nlfuns.exp; % inverse link function (a.k.a. nonlinearity)
            lfunc = @(w)(glms.neglog.poisson(w, dm.X, y, fnlin)); % cost/loss function
            
            opts = optimoptions(@fminunc, 'Algorithm', 'trust-region', ...
                'GradObj', 'on', 'Hessian','on');%,'UseParallel',true);
            
            try
                [wml, nlogli, exitflag, ostruct, grad, hessian] = fminunc(lfunc, wInit, opts);
                w(:,ifold) = wml;
            catch err
                w(:,ifold) = wInit;
            end
            
            %wvar = diag(inv(hessian));
        case 5 %ridge regression used in 2022
            w(:,ifold)= rReg(dm.X(:,2:end), y, ridgeParam, useGPU);
            
    end
    
    %% Simulate from model for test data
    % dmTest = buildGLM.compileSparseDesignMatrix(dspec, testTrialIndices);
    % yPred = generatePrediction(w, model, dmTest); %does not exist
    dm_pred = buildGLM.compileSparseDesignMatrix(dspec, xvFolds{ifold,2});
    dm_pred = buildGLM.addBiasColumn(dm_pred);
    
    switch option
        case {1,5}
            yPred_xv = dm_pred.X*w(:,ifold);
            %yPred_xv  = dm_pred.X*w(2:end,ifold)+w(1,ifold);
        case {2,3,4}
            yPred_xv = exp(dm_pred.X*w(:,ifold)); %for option3
    end
    if ~useSptrain
        yPred_xv = yPred_xv/dt_r; %[Hz]
    end
    
    tidx_r_xv=[];
    for itr = xvFolds{ifold,2}
        tidx_r_xv = [tidx_r_xv; trIdx_r{itr}];
    end
    predicted(tidx_r_xv) = yPred_xv;
    
    [expval(ifold), mse(ifold), R(ifold)] = getExpVal(PSTH_f(tidx_r_xv), ...
        predicted(tidx_r_xv)+mean(PSTH_f(tidx_r_xv))-mean(predicted(tidx_r_xv)));
    %adjust mean of predicted for the case of gradual increase of baseline firing ...kind of cheating
    
    
    head=1;%set 0 if intercept is not estimated
    for ivar = 1:numel(npredVars)
        dspec_sub = dspec;
        dspec_sub.covar = dspec.covar(ivar);
        dspec_sub.edim = nBasisFunctions(ivar)*npredVars(ivar);
        dm_pred_sub = buildGLM.compileSparseDesignMatrix(dspec_sub, xvFolds{ifold,2});
        %dm_pred_sub = buildGLM.addBiasColumn(dm_pred_sub);
        
        widx = (1:nBasisFunctions(ivar)*npredVars(ivar))+head;
        head = max(widx);
        switch option
            case {2,3,4}
                yPred_xv_sub = exp(dm_pred_sub.X*w(widx,ifold)+w(1,ifold)); %for option3
            case {1,5}
                yPred_xv_sub = dm_pred_sub.X*w(widx,ifold)+w(1,ifold);
        end
        if ~useSptrain
            yPred_xv_sub = yPred_xv_sub/dt_r; %[Hz]
        end
        predicted_each(tidx_r_xv,ivar) = yPred_xv_sub;
    end
end

% Post regression weight reconstruction
mw = mean(w,2); %avg across folds. cf. cmdExample
ws = buildGLM.combineWeights(dm, mw);
%ws.(covLabel).tr: time axis
%ws.(covLabel).data: kernel
r0 = mw(1);% intercept

% dm = buildGLM.compileSparseDesignMatrix(dspec, 1:numel(trIdx_r));
% dm.X = full(dm.X);
% dm = buildGLM.addBiasColumn(dm);
% yPred = exp(dm.X*mw);

for ivar = 1:numel(npredVars)
    kernelInfo.kernel{ivar} = ws.(predictorNames{ivar}).data;
    kernelInfo.tlags{ivar} = ws.(predictorNames{ivar}).tr;
end
kernelInfo.intercept = r0;
kernelInfo.fs = 1/dt_r;
kernelInfo.ridgeParam = ridgeParam;
kernelInfo.mse = mean(mse);
kernelInfo.expval = mean(expval);
kernelInfo.corrcoef = mean(R);
kernelInfo.cv.kernel = w;
kernelInfo.cv.mse = mse;
kernelInfo.cv.corrcoef = R;
kernelInfo.cv.expval = expval;

