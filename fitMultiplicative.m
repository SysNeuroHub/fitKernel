function [predicted_cue, gainInfo] = fitMultiplicative(PSTH_f, predicted_all, ...
    t_r, predictors_r)

omitDuration = 0;%5; %omit initial and last segments for fitting[s]

nVar = 8;%predictorInfo_cue.npredVars;
regIdx = intersect(find(t_r>=t_r(1)+omitDuration), find(t_r<=t_r(end)-omitDuration));
% timeVec = t_r(regIdx)';
%observed = PSTH_f(regIdx) - predicted_slow(regIdx);
observed = log(PSTH_f(regIdx)) - log(predicted_all(regIdx));
% observed(imag(observed)~=0) = 0;%nan?
incIdx = intersect(find(imag(observed)==0),find(~isinf(observed)));
observed = observed(incIdx);
predictor = predictors_r(:,regIdx(incIdx));

ridgeParam = 0;
%b = [ones(numel(observed),1) observed] \ predictor(ibhv,:)';
b = rReg(predictor', observed(:), ridgeParam, 0);%does not match with the result with \
% b = ridge(observed(:), predictor', ridgeParam, 0);

% plot(exp(b(2:end)));

%% check effect of cue (without fitting)
% ibhv=1;
% histogram(observed(predictor(ibhv,:)==0));
% hold on
% histogram(observed(predictor(ibhv,:)==1));
% 
% nanstd(observed(predictor(ibhv,:)==0))
% nanstd(observed(predictor(ibhv,:)==1))


r0 = b(1);
rr = b(2:end);

gainInfo.r0 =r0;
gainInfo.rr = rr;

predicted_log = r0 + log(predicted_all);
for ibhv = 1:nVar
    predicted_log = predicted_log + predictors_r(ibhv,:)'*rr(ibhv);
end
predicted_cue = exp(predicted_log);


%% sanity check
% plot(t_r, PSTH_f, t_r, predicted_all, t_r, predicted_cue);
% [startIdx, endIdx, cueIdx] = trace2box(predictorInfo_cue.predictors_r,[],'up');
% vbox(t_r(startIdx), t_r(endIdx), gca);%, @hsv, cueIdx);
% legend('observed','GLM','GLM+gain');
% 

