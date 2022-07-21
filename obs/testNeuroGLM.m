%from https://github.com/pillowlab/neuroGLM/blob/master/docs/tutorial.md

unitOfTime = 's';
binSize = 0.02;
uniqueID = 1;

%t_r;%

expParam = [];
expt = buildGLM.initExperiment(unitOfTime, binSize, uniqueID);

%Registering variables to the experiment
expt = buildGLM.registerContinuous(expt, 'vision', 'Target Position', 8); %8dim
expt = buildGLM.registerContinuous(expt, 'eyepos', 'Eye Position', 8); 
expt = buildGLM.registerContinuous(expt, 'pdiam', 'Pupil Diamter', 1); 
expt = buildGLM.registerContinuous(expt, 'blink', 'Blink Event', 1); 
expt = buildGLM.registerSpikeTrain(expt, 'sptrain', 'Our Neuron'); % Spike train!!!

%Loading the data for each trial
duration = t_r(end);
trial = buildGLM.newTrial(expt, duration);

trial.dotson = rand() * duration; % timing variable WHAT IS THIS FOR?
trial.sptrain = spk_all_cat(spk_all_cat>0&spk_all_cat<t_r(end));

trial.vision = predictorInfo.predictors_r(1:8,t_r>0)';
trial.eyepos = predictorInfo.predictors_r(9:16,t_r>0)';
%T = expt.binfun(trial.duration); % number of bins for this trial
%trial.eyepos = randn(T, 1); % seems like this field is mandatory?
trial.pdiam = predictorInfo.predictors_r(17,t_r>0)';
trial.blink = predictorInfo.predictors_r(18,t_r>0)';

%add the trial object to the experiment object with an associated trial index kTrial:
kTrial = 1;
expt = buildGLM.addTrial(expt, trial, kTrial);

% creating a design specification object.
dspec = buildGLM.initDesignSpec(expt);

%a set of 8 boxcar basis functions to cover 0.3 s 
%bs = basisFactory.makeSmoothTemporalBasis('boxcar', .3, 8, expt.binfun);

nBasisFunctions = 20;
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


trialIndices = 1;
dm = buildGLM.compileSparseDesignMatrix(dspec, trialIndices);

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
rho=.2; %glmspike.m
ndx = []; %glmspike.m
M = regression.doRegressionPoisson(dm.X, y, dspec, ndx, dt_r, rho); %WARNING: SLOW
w = M.khat;

% Post regression weight reconstruction
ws = buildGLM.combineWeights(dm, M.khat);
%ws.(covLabel).tr: time axis
%ws.(covLabel).data: kernel
r0 = w(1);% intercept


%% Simulate from model for test data
% dmTest = buildGLM.compileSparseDesignMatrix(dspec, testTrialIndices);
% yPred = generatePrediction(w, model, dmTest); %does not exist
yPred = exp(dm.X*w); %for option3

