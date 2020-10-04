function [S_XX,freq] = fftPSDEstimate(Xn,fs,pad)
% fftPSDEstimate: Estimates the Power Ppectral Density of the signal Xn
% with a given frequensy fs, by using preprogrammed fft-functions. 
%  
% Input:
%  Xn   - Discrite signal
%  fs   - Sampling frequency [Hz]
%  pad  - Number of padding samples. If no padding is given, padding is
%         sett to zero
%
% Output: 
%  S_XX - Power spectral density vector
%  freq - Frequency vector with frequencies |f| < f_NQ (Nyquist frequency)

% Number of samples
N = length(Xn);

% Checks if the padding variable is empty, yes pad = 0
if nargin < 3 || isempty(pad)
    pad = 0;
end

% Number of samples plus padding
M = N+pad;

% Distance between sample frequencies in time
Dt = 1/fs;

% Nyquist frequency
f_NQ = fs/2;

% Estimated power spectral density
Density = Dt/M*abs( fft(Xn,M) ).^2;

% Since fft calculates the Fourier transform for positive Fourier
% frequencies, we must use the fftshift-command to be able to plot the
% periodogram estimate in the frequency interval |f| < f_NQ.
S_XX = fftshift(Density);

freq = linspace(-f_NQ, f_NQ, M);
end