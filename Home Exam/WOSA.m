function [S_XX,freq] = WOSA(Xn,window,fs,overlap,Nfft,DataWindow)
% WOSA: This function estimates the Weighted Overlapped Segment Averaging. 
%  
% Input:
%  Xn         - Discrite signal vector.
%  window     - The length of each segment [positive intenger].
%  fs         - Sampling frequency [Hz]
%  overlap    - Number of overlapping samples [positive intenger]. If no
%               overlap is introduced, overlap is sett to zero
%  Nfft       - Length of the FFT, so zero padding would be done when
%               Nfft > window.
%  DataWindow - Data window used to multiply the data set prior to the
%               Fourier transform [vector with the same length as the
%               window]. If no data window is given then the modified
%               periodograms are computed using a Hamming window.
%
% Output: 
%  S_XX - Power spectral density vector
%  freq - Frequency vector with frequencies |f| < f_NQ (Nyquist frequency)

% Checks if the overlap variable is empty, if yes sett overlap = 0
if nargin < 4 || isempty(overlap)
    O = 0;
else
    O = overlap;
end

% Checks if the Nfft variable is empty, if yes sett Nfft = window
if nargin < 5 || isempty(Nfft)
    Nfft = window;
end

% Number of samples
N = length(Xn);
% Window length
M = window;
% Number of segments
L = ceil(N/(M-O));

% Checks if the DataWindow is empty, if yes use Hamming window
if nargin < 6|| isempty(DataWindow)
    m = 0:1:M-1;
    
    % Hamming window with window length M
    Wm = 0.54 - 0.46 * cos(2*pi*m/(M-1));
else
    Wm = DataWindow;
end
% Normalization factor of Wm
U = 1/M * sum((Wm).^2);

% Distance between sample frequencies in time
Dt = 1/fs;
% Nyquist frequency
f_NQ = fs/2;


% We begin by segmenting the data set into K segments each consisting of M
% samples, and save it as a matrix where every row correspond to a segment
Xk = zeros(L,M);
for k = 1:L
    Ssegm = (M-O)*(k-1) + 1;  % Starting point for segment number k 
    Esegm = (M-O)*(k-1) + M;  % Ending point for segment number k
    
    % If the end point is greater then the length of signal Xn, then save
    % the remaining values Xn on Xk and let the rest be zero.
    if Esegm > N
        Esegm = N;
        Xk(k,1:(Esegm-Ssegm+1)) = Xn(Ssegm:Esegm);
    else
        Xk(k,:) = Xn(Ssegm:Esegm);
    end
end

% Computing the modified spectral estimate 
P = zeros(L,Nfft);
for k = 1:L
    P(k,:) = Dt/(M*U)*abs( fft(Wm.*Xk(k,:),Nfft) ).^2; 
end

% Average periodogram
AverageP = 1/L*sum(P);

% Shift
S_XX = fftshift(AverageP);

freq = linspace(-f_NQ, f_NQ, Nfft);
end