function [S_XX,freq] = PSDEstimate(Xn,fs,pad)
% PSDEstimate: Estimates the Power Ppectral Density of the signal Xn
% with a given frequensy fs, and pad it . 
%  
% Input:
%  Xn   - Discrite signal
%  fs   - Sampling frequency [Hz]
%  pad  - Number of padding samples. If no padding is given, padding is
%         sett to zero.
%
% Output: 
%  S_XX - Power spectral density vector
%  freq - Frequency vector with frequencies |f| < f_NQ (Nyquist frequency)

% Number of samples
N = length(Xn);

% Checks if the padding variable is empty, if yes sett pad = 0
if nargin < 3 || isempty(pad)
    pad = 0;
end

% Number of samples plus padding
M = N+pad;

% Discrete time index vector
n = 0:1:M-1;

% Distance between sampled frequencies
Df = fs/M;

% Distance between sample frequencies in time
Dt = 1/fs;

% Nyquist frequency
f_NQ = fs/2;

% Padd Xn with a trail zeros to length M, if padding is greater then 0
if pad > 0
    Xn(N+1:M) = 0;
end

% Estimated power spectral density
Density = zeros(1,M);
for k = 1:M
    Density(k) = Dt/M*abs( sum(Xn.*exp(-1i*2*pi*(k-1)*Df*n*Dt)) ).^2;
end

% check if M is even or odd and shift accordingly. The mod-function denotes
% the rest between M and 2, if the rest is zero then M is a even number.
if mod(M,2) == 0
    S_XX = circshift(Density,M/2);
end
if mod(M,2) == 1
    S_XX = circshift(Density,(M-1)/2);
end

freq = linspace(-f_NQ, f_NQ, M);
end