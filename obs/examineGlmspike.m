%% from glmspike

gammaStim = 500;

% from constructor
obj.designOptions.binSize            = 10;

% from setDesignOptions
obj.binSize = obj.designOptions.binSize;
obj.binfun = str2func(['@(t) (t==0) + ceil(t/' num2str(obj.binSize) ')']);

% from buildDesignBoxcar
stimBasis = basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', ...
    [30 gammaStim],10, obj.binfun, 100);
stimBasis.normalize;
                    
% targBasis = basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', ...
%     [30 gammaStim], 8, obj.binfun, 10);
% targBasis.normalize;
% % add Targets
% obj.addCovariateTiming(trial, 'targson', 'targson', 'Targets Onset', targBasis)

% --- add contrast using spatially averaged contrast
if obj.binSize==1
    stimHandle = @(trial) trial.contrast;
else
    stimHandle = @(trial) basisFactory.binDataVector(trial.contrast/obj.binSize, obj.binSize);
end
obj.addCovariate(trial, 'Motion', 'Motion Contrast', stimHandle, stimBasis);
% in neuroGLM.addCovariate, convolve stimulus and basis functions (should reduce dimension)


% from M=fitRidgeFixed(obj, trial, rho)
rho=.2;

y=obj.getBinnedSpikeTrain(trial, obj.fitNeuron, obj.dm.trialIndices);

%   dspec [struct] or [neuroGLM object]
%       .covar - array of covariates
%           .label
%           .edim
%           .desc
%       .edim  - total dimensionality
%       .unitoftime
%       .binSize
%       .model
%           .regressionMode
%           .bilinearMode
%           .bilinearRank
%           .bilinearCovariate
obj.covar %< neuroGLM.addCovariate
obj.edim = sum([obj.covar(:).edim]); %< neuroGLM.addCovariate
obj.unitoftime %property of neuroGLM: 'string' - 's' or 'ms' indicating a global unit of time
obj.binSize %property of neuroGLM:  - duration of each time bin in units of unitOfTime
obj.model.regressionMode='RIDGEFIXED';
obj.model.bilinearMode%optional field. If ON,  Do bilinear optimization (not used in glmspike)
dspec = obj;
M=regression.doRegressionPoisson(obj.dm.X, y, dspec, [], obj.binSize/1e3, rho);

lambda=M.fnlin(obj.dm.X*M.khat)*M.dt;
M.logli=logliPoisson(lambda, y);
M.df=obj.edim;
M.AIC=2*obj.edim - 2*M.logli;

