function [predicted, predicted_each, observed, kernelInfo] = (spk_cat, ...
    t_r, predictors_r, npredVars, sigma, lagRange, ridgeParam)
%from https://github.com/pillowlab/neuroGLM/blob/master/docs/tutorial.md

unitOfTime = 's';
binSize = 0.02;
uniqueID = 1;
%predictorInfo.npredVars = [8 8 1 1];
nBasisFunctions = 20;
KFolds = 5;
ridgeParam=.2; %glmspike.m %ridgeparameter

%t_r;%

expParam = [];
expt = buildGLM.initExperiment(unitOfTime, binSize, uniqueID);

%Registering variables to the experiment
expt = buildGLM.registerContinuous(expt, 'vision', 'Target Position', predictorInfo.npredVars(1)); %8dim
expt = buildGLM.registerContinuous(expt, 'eyepos', 'Eye Position', predictorInfo.npredVars(2));
expt = buildGLM.registerContinuous(expt, 'pdiam', 'Pupil Diamter', predictorInfo.npredVars(3));
expt = buildGLM.registerContinuous(expt, 'blink', 'Blink Event', predictorInfo.npredVars(4));
expt = buildGLM.registerSpikeTrain(expt, 'sptrain', 'Our Neuron'); % Spike train!!!


%Loading the data for each trial
[trIdx_r, trIdx] = retrieveTrIdx_r(t_cat, t_r, t_tr);
for itr = 1:numel(t_tr)
    t_cat_tr = t_r(trIdx_r{itr});
    above0 = find(t_cat_tr>=0);
    trIdx_r{itr} = trIdx_r{itr}(above0);
    t_cat_tr = t_cat_tr(above0);
    %duration = t_cat_tr(end) - t_cat_tr(1) + dt_r;
    duration = numel(trIdx_r{itr})*dt_r;
    
    trial = buildGLM.newTrial(expt, duration);
    
    %trial.dotson = rand() * duration; % timing variable WHAT IS THIS FOR?
    trial.sptrain = spk_all_cat(spk_all_cat>=t_cat_tr(1) & spk_all_cat<=t_cat_tr(end)) ...
        - t_cat_tr(1); %must start from 0
    
    trial.vision = predictorInfo.predictors_r(1:8,trIdx_r{itr})';
    trial.eyepos = predictorInfo.predictors_r(9:16,trIdx_r{itr})';
    %T = expt.binfun(trial.duration); % number of bins for this trial
    %trial.eyepos = randn(T, 1); % seems like this field is mandatory?
    trial.pdiam = predictorInfo.predictors_r(17,trIdx_r{itr})';
    trial.blink = predictorInfo.predictors_r(18,trIdx_r{itr})';
    
    %add the trial object to the experiment object with an associated trial index kTrial:
    expt = buildGLM.addTrial(expt, trial, itr);
end

% creating a design specification object.
dspec = buildGLM.initDesignSpec(expt);

kernelDur = diff([-0.5 0.5]);
basisType = 'raised cosine';
bs=basisFactory.makeSmoothTemporalBasis(basisType, kernelDur, nBasisFunctions, ...
    expt.binfun);

% from Yate's neuroGLM:
% basisType = 'nonlinearly scaled cosine';
% nloffset = 1.5;%?
% bs=basisFactory.makeSmoothTemporalBasis(basisType, lagRange, nBasisFunctions, ...
%     expt.binfun, nloffset);
% bs.normalize

%% add covariate without basis functions
%dspec = buildGLM.addCovariateRaw(dspec, 'vision', 'vision effect', bs);
%dspec = buildGLM.addCovariateRaw(dspec, 'eyepos', 'Eye position effect', bs);

%dspec = buildGLM.addCovariateTiming(dspec, 'eyepos', 'eyepos', bs);
%< Type of label [eyepos] is not timing

offset = -25;
covLabel = 'vision';
stimHandle = basisFactory.rawStim(covLabel);
dspec = buildGLM.addCovariate(dspec, covLabel,covLabel, stimHandle, bs, offset);

covLabel = 'eyepos';
stimHandle = basisFactory.rawStim(covLabel);
dspec = buildGLM.addCovariate(dspec, covLabel,covLabel, stimHandle, bs, offset);

covLabel = 'pdiam';
stimHandle = basisFactory.rawStim(covLabel);
dspec = buildGLM.addCovariate(dspec, covLabel,covLabel, stimHandle, bs, offset);

covLabel = 'blink';
stimHandle = basisFactory.rawStim(covLabel);
dspec = buildGLM.addCovariate(dspec, covLabel,covLabel, stimHandle, bs, offset);



xvFolds = regression.xvalidationIdx(numel(t_tr), KFolds, false, true);

w = [];
predicted = zeros(numel(t_r),1);
for ifold = 1:KFolds
    dm = buildGLM.compileSparseDesignMatrix(dspec, xvFolds{ifold,1});
    
    %If your design matrix is not very sparse (less than 10% sparse, for example),
    %it's better to conver the design matrix to a full (dense) matrix for speed.
    dm.X = full(dm.X);
    
    % rho = 1; % ridge parameter (play around with this to see how it effects fitting. Use cross validation on the training set to set it properly)
    % dm.addBiasColumn('right'); % augment with column of ones
    dm = buildGLM.addBiasColumn(dm);
    
    %% Get the spike trains back to regress against
    y = buildGLM.getBinnedSpikeTrain(expt, 'sptrain', dm.trialIndices);
    
    % Doing the actual regression
    % option1: simple least squares
    % w = dm.X' * dm.X \ dm.X' * y;
    % option2: Maximum likelihood estimation using glmfit
    %[w, dev, stats] = glmfit(dm.X, y, 'poisson', 'link', 'log');
    % option3: ridge regression with static nonlinearity
    dspec.model.regressionMode='RIDGEFIXED';
    dspec.model.nlfun = @expfun;
    dt_r = median(diff(t_r));
    ndx = []; %glmspike.m
    M = regression.doRegressionPoisson(dm.X, y, dspec, ndx, dt_r, ridgeParam); %WARNING: SLOW
    w(:,ifold) = M.khat;
    
    %% Simulate from model for test data
    % dmTest = buildGLM.compileSparseDesignMatrix(dspec, testTrialIndices);
    % yPred = generatePrediction(w, model, dmTest); %does not exist
    dm_pred = buildGLM.compileSparseDesignMatrix(dspec, xvFolds{ifold,2});
    dm_pred = buildGLM.addBiasColumn(dm_pred);

    yPred_xv = exp(dm_pred.X*w(:,ifold)); %for option3

    tidx_r_xv=[];
    for itr = xvFolds{ifold,2}
        tidx_r_xv = [tidx_r_xv; trIdx_r{itr}];
    end
    predicted(tidx_r_xv) = yPred_xv;
    
    head=1;%set 0 if intercept is not estimated
    for isub = 1:4
        dspec_sub = dspec;
        dspec_sub.covar = dspec.covar(isub);
        dspec_sub.edim = nBasisFunctions*predictorInfo.npredVars(isub);
        dm_pred_sub = buildGLM.compileSparseDesignMatrix(dspec_sub, xvFolds{ifold,2});
        %dm_pred_sub = buildGLM.addBiasColumn(dm_pred_sub);
        
        widx = (1:nBasisFunctions*predictorInfo.npredVars(isub))+head;
        head = max(widx);
        yPred_xv_sub = exp(dm_pred_sub.X*w(widx,ifold)); %for option3
    
        predicted_each(tidx_r_xv,isub) = yPred_xv_sub;
    end
end


PSTH_r = getPSTH(spk_all_cat, t_r);
PSTH_f = filtPSTH(PSTH_r, dt_r, sigma, 2);%TO BE FIXED

% Post regression weight reconstruction
mw = mean(w,2); %avg across folds. cf. cmdExample
ws = buildGLM.combineWeights(dm, mw);
%ws.(covLabel).tr: time axis
%ws.(covLabel).data: kernel
r0 = mw(1);% intercept

% dm = buildGLM.compileSparseDesignMatrix(dspec, 1:numel(t_tr));
% dm.X = full(dm.X);
% dm = buildGLM.addBiasColumn(dm);   
% yPred = exp(dm.X*mw);

observed = PSTH_f;
mse = mean((observed - predicted).^2);
expval = 100*(1 - mse / mean((observed - mean(observed)).^2));
R = corrcoef(observed, predicted);

kernelInfo.intercept = r0;
kernelInfo.fs = 1/dt_r;
%kernelInfo.tlags = tlags;
kernelInfo.ridgeParam = ridgeParam;
kernelInfo.mse = mse;
kernelInfo.expval = expval;
kernelInfo.corrcoef = R(1,2);
