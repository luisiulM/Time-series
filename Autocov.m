function [ ACV ] = Autocov(Xn,Nlag)
% Autocov: This fuction estimates the autocovariance of the input signal.
%
% Input:
%  Xn   - Discrite signal vector.
%  Nlag - Total number of lag points. Note that Nlag must be lower than the
%         numnber of samples N.
%
% Output: 
%  ACV  - Estimated autocovariance of Xn

% Numebr of samples
N = length(Xn);

% Checks if the Nlag variable is empty, if yes make the total number of lag
% points equal to the numner of samples
if nargin < 2 || isempty(Nlag)
    Nlag = N;
end

% Estimating the autocovariance for each lag h
Estimate = zeros(1,Nlag);
for h=1:Nlag    
    Estimate(h) = sum(( Xn((1:(N-h+1))+(h-1)) - mean(Xn) ).*(( Xn((1:(N-h+1)))-mean(Xn) )));  
end

ACV = Estimate/N;
end
