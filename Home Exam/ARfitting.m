function [ prediction ] = ARfitting(Xn,p)
% ARfitting: This function estimates the prediction of a autoregressive
% random process of order p (AR(p)), of an input signal Xn. For examble 
%
% Input:
%  Xn         - Discrite signal vector.
%  p          - The order of the autoregression model.
%
% Output: 
%  prediction - Autoregression model of the signal Xn with the order p 

% Estimating the autocovariance of the Music signal and making a mirror
% sided vector. That is, if the autocovariance function gives us the values
% [a0 a1 a2 ... ap+1], than Autocov vector is equal to 
% [ap+1 ... a2 a1 a0 a1 a2 ... ap+1]
AutocovX = fliplr(Autocov(Xn,p+1));
AutocovX(p+1:(2*p+1)) = Autocov(Xn,p+1);

% Estimating the covariance matrix and vector
P = AutocovX(p+2:(2*p+1));
R = zeros(p,p);
for i = 1:p
    R(i,:) = AutocovX((p+2-i):(2*p+1-i));
end

% Estimating the coefficients by using Yule-Walker equation.
Coeff = inv(R)* transpose(P);

% Now that we have estimated the coefficients, we can use them to 
% estimate a AR(p).

% Number of samples
N = length(Xn);

% Here we estimate the backstep components and save them as rows. That is
% X[n-3] is on row 3 and X[n-p] is on row p
Xp = zeros(p,p);
for j = 1:p
    Xp(j,(j+1):N) = Xn(1:(N-j));
end

% Estimating the AR(p) series in correspondence to the music data
prediction = transpose(Coeff)*Xp;
end

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