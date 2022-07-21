%testing doRegressionPoisson.m - to replace with fitPSTH

%% INPUTS
t_r = predictorInfo.t_r;
sigma = param.psth_sigma;
lagRanges = param.lagRange;

dt_r = median(diff(t_r)); %[ms]

PSTH_r = getPSTH(spk_all_cat, t_r);
PSTH_f = filtPSTH(PSTH_r, dt_r, sigma, 2);%causal

%% compile design matrix from ridgeXs. cf. neuroGLM.compiledesignmatrix
predictor = predictorInfo.predictors_r;
nt = size(predictor,2);
nVar = size(predictor,1);
fs = 1/dt_r;
%tlags = lags/fs;
lags = round(min(lagRanges(:,1))*fs):round(max(lagRanges(:,2))*fs);
nLags = length(lags);
nVar_lag = size(lagRanges,1);

localLags = cell(nVar,1);
for iVar = 1:nVar_lag
    localLags{iVar} = round(lagRanges(iVar,1)*fs):round(lagRanges(iVar,2)*fs);
end

    

X = zeros(nt, nVar*nLags);
nLags_tmp = nLags;
for ibhv = 1:nVar
    for iLag = 1:nLags
        %X(:, iLag) = circshift( predictor(:), lags(iLag)-1);
        if ismember(lags(iLag),localLags{ibhv}) || nVar_lag == 1
            X(:, iLag+(ibhv-1)*nLags) = circshift( predictor(ibhv,:), lags(iLag)-1);
        else
            %values outside the lagRange is substituted by 0 ... is
            %this correct?
            X(:, iLag+(ibhv-1)*nLags) = zeros(1,nt); %nans did not work
        end
    end
end
    
%for estimation of the intercept (baseline firing rate)
X = [ones(size(X, 1), 1), X]; %from neuroGLM.addbiascolumn

%% run doRegressionPoisson
%   X [nsamples x ndims] Stimulus
%   Y [nsamples x 1]     Count
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
%           .nlfun: static nonlinearity (default:@expfun)
%           .bilinearMode (optional)
%           .bilinearRank (optional)
%           .bilinearCovariate (optional)
y = PSTH_f;
rho=.2; %glmspike.m
ndx = []; %glmspike.m

dspec = [];
for ibhv = 1:nVar
     newIdx = ibhv;%numel(fieldnames(dspec.idxmap)) + 1;
    %dspec.addCovariate(trial, 'Motion');
    covLabel = 'dum';
    desc = covLabel;
    stimHandle = [];
    dspec.covar(newIdx)=struct('label', covLabel, 'desc', desc, 'stim', stimHandle, ...
        'offset', 0, 'cond', [], 'basis', [], 'edim', [], 'sdim', []);

    %dspec.edim = sum([dspec.covar(:).edim]);
end
dspec.edim = nVar;
dspec.unitoftime = 's'; %property of neuroGLM: 'string' - 's' or 'ms' indicating a global unit of time
dspec.binSize = dt_r;%property of neuroGLM:  - duration of each time bin in units of unitOfTime
dspec.model.regressionMode='RIDGEFIXED';
dspec.model.nlfun = @expfun;


M = regression.doRegressionPoisson(X, y, dspec, ndx, dt_r, rho); %WARNING: SLOW

lambda=M.fnlin(X*M.khat)*M.dt;
M.logli=logliPoisson(lambda, y);
M.df=dspec.edim;
M.AIC=2*dspec.edim - 2*M.logli;


%% retrieve karnel weights cf neuroGLM.combineWeights(obj, w)
r0 = M.khat(1);
rr = reshape(M.khat(2:end),nLags,nVar);


%% predict psth cf. glmspike.predictSpikeRate
Xc = reshape(X(:,2:end), size(X,1), nLags, nVar);
predicted = r0;
for ibhv = 1:nVar
    predicted = predicted + Xc(:,:,ibhv)*rr(:,ibhv);
end
 predicted = exp(predicted')*dt_r; 

mse = mean((y - predicted').^2);
expval = 100*(1 - mse / mean((y - mean(y)).^2));
R = corrcoef(y, predicted');

kernelInfo.kernel = rr;
kernelInfo.tlags = lags/fs;
kernelInfo.expval = expval;
kernelInfo.corrcoef = R(1,2);


%% figure for kernel fitting
figure('position',[0 0 1000 500]);
subplot(1,2,1);
plot(predictorInfo.t_r, PSTH_f, 'color',[.5 .5 .5]);hold on
plot(predictorInfo.t_r, predicted, 'linewidth',2);
%xlim([1935 1966]);
xlim([100 200])
legend('recorded','fitted');
xlabel('time [s]'); ylabel('firing rate [Hz]');

title(['expval: ' num2str(kernelInfo.expval), ', R: ' num2str(kernelInfo.corrcoef)]);

a2=subplot(3,2,2);
thisIm = kernelInfo.kernel(:,1:predictorInfo.npredVars(1))';
%         crange = prctile(abs(thisIm(:)),99);
crange = prctile(thisIm(:),[1 99]);
imagesc(kernelInfo.tlags, param.cardinalDir, thisIm);
caxis([crange]);
set(gca,'ytick',param.cardinalDir);
xlabel('time from targetOnset [s]');
mcolorbar(a2,.5);

a3=subplot(3,2,4);
thisIm = kernelInfo.kernel(:,predictorInfo.npredVars(1)+1:sum(predictorInfo.npredVars(1:2)))';
%crange = prctile(abs(thisIm(:)),99);
crange = prctile(thisIm(:),[1 99]);
imagesc(kernelInfo.tlags,param.cardinalDir, thisIm);
caxis([crange]);
set(gca,'ytick',param.cardinalDir);
xlabel('time from eye movement [s]');
mcolorbar(a3,.5);

a4=subplot(3,2,6);
plot(kernelInfo.tlags, kernelInfo.kernel(:,sum(predictorInfo.npredVars(1:2))+1:end)');
xlabel('time from pupil dilation/blink [s]');
axis tight;

%set(gca,'ytick',predictorInfo.predictors_r);
screen2png(fullfile(saveFigFolder,['kernels_exp' saveSuffix '_yates']));
