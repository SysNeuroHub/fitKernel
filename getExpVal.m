function [expval, mse, corr] = getExpVal(observed, predicted)
%[expval, mse, corr] = getExpVal(observed, predicted)

predicted = predicted - mean(predicted) + mean(observed); %13/6/24

mse = mean((observed - predicted).^2);
expval = 100*(1 - mse / mean((observed - mean(observed)).^2));
R = corrcoef(observed, predicted);
corr = R(1,2);
