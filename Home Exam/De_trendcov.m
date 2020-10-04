function [ ACV ] = De_trendcov(Yn)
% De_trendcov: This fuction estimates the autocovariance of the de-trended
% input signal, where the de-trending method is commonly known as
% X(n) = Y(n)-Y(n-1).
%
% Input:
%  Yn  - Discrite signal vector.
%
% Output: 
%  ACV - Estimated autocovariance of the de-trended signal Xn = Yn - Y(n-1)

% Number of samples
N = length(Yn);

% Implimenting the backstep term Y(n-1)
Yn1(1) = 0;
Yn1(2:N) = Yn(1:(N-1));

% De-trending method X(n) = Y(n)-Y(n-1)
Xn = Yn - Yn1;

plot(Xn)

% Estimating the autocovariance for each lag h
Estimate = zeros(1,N);
for h=1:N    
    Estimate(h) = sum(( Xn((1:(N-h+1))+(h-1)) - mean(Xn) ).*(( Xn((1:(N-h+1)))-mean(Xn) )));  
end

% Normalize and return the value
ACV = Estimate/N;
end