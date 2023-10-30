function [AdjR2, TSS, RSS] = getRsqadj(y, yHat, p)
%[AdjR2, TSS, RSS] = getRsqadj(y, yHat, p)
%
%https://au.mathworks.com/help/stats/coefficient-of-determination-r-squared.html
%
% y: observed trace
% yHat: predicted trace
% p: number of predictors

%[n, p] = size(X); % Number of samples and features
n = numel(y); %Number of observations
TSS = sum((y - mean(y)).^2); % Total sum of squares

RSS = sum((y - yHat).^2); % Residual sum of squares
AdjR2 = 1 - (RSS/(n-p-1)) * (n-1)/TSS; % Adjusted R-squared
