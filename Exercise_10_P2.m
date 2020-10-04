%% 2a) A popular data window is the Hann (Hanning) window. Plot the Hann 
% window.

N = 100;
n = 0:1:N-1;

w_n =1/2*(1 - cos(2*pi*n/(N-1))); % Hann (Hanning) window

plot(n, w_n)
ylabel('Hann window w_n')
xlabel('n discrete steps')

%% 2b) Estimate analytically W(f) = DTFT{w_n} of the Hann window, and plot
% |W(f)|^2 (corresponding to Q(f), which we have seen during lectures) on
% a dB-scale.

Dt = 1;          % Difference in time steps
f_NQ = 1/(2*Dt); % Nyquist frequency
freq = linspace(-f_NQ, f_NQ, N);

W = fft(w_n);
W_shift = fftshift(W);
W_mag = abs(W_shift).^2;

dB_W = mag2db(W_mag);

plot(freq, dB_W)
xlabel('Frequencies')
ylabel('dB scalling')

%% 2c) Plot the periodogram and the Hann window periodogram, and comment
% on the differences.

% Parameters
N = 32;                      % Number of samples
n = 0:1:N-1;
Dt = 1;                      % Difference in time steps
f_NQ = 1/(2*Dt);
freq = linspace(-f_NQ, f_NQ, N);
theta = unifrnd(0,2*pi,[1,1]);  % Uniform distributed variable

sigma = 0.5;
w_t = normrnd(0,sigma^2,[1,N]);  % white noise

A = 2;                         % Amplitude
f0 = 0.1;                      % Frequency [s^-1]
%

X_t = zeros(1,N);
for i = 1:N
    X_t(i) = A*cos(2*pi*f0*(i-1) + theta) + w_t(i); % Periodic process
end

w_n =1/2*(1 - cos(2*pi*n/(N-1))); % Hann window
U = 1/N * sum((w_n).^2);

H = zeros(1,N);
for j = 1:N
   H(j) = X_t(j)*w_n(j); 
end

S_XX = Dt/N*abs( fft(X_t) ).^2;
S_XXshift = fftshift(S_XX);
dB_X = mag2db(S_XXshift);          % dB scale

S_Han = Dt/(N*U)*abs( fft(H) ).^2; % Hanning Periodogram estimate
S_Hanshift = fftshift(S_Han);
dB_W = mag2db(S_Hanshift);          % dB scale

subplot(2,1,1)
plot(freq, dB_X)
xlabel('Frequencies  [s^{-1}]' )
ylabel('Hanning Periodogram estimate S_{XX}')

subplot(2,1,2)
plot(freq, dB_W)
xlabel('Frequencies  [s^{-1}]')
ylabel('dB scale of S_{XX}')

% Comment: 

%% 2d) The Hamming window is also well-known. Plot this window.

N = 100;
n = 0:1:N-1;

w_m = 0.54 - 0.46 * cos(2*pi*n/(N-1)); % Hamming window

plot(n, w_m)
ylabel('Hamming window w_n')
xlabel('n discrete steps')

%% 2e) Estimate analytically W(f) = DTFT{w_n} of the Hamming window, and plot
% |W(f)|^2 on a dB{scale. Compare with the result in (b).
% run: 2a), 2b), 2d)
W_m = fft(w_m);
W_mshift = fftshift(W_m);
W_mmag = abs(W_mshift).^2;

dB_Wm = mag2db(W_mmag);

plot(freq, dB_W)
xlabel('Frequencies')
ylabel('dB scale')
hold on
plot(freq, dB_Wm)
legend('Hanning window','Hamming window')

%% 2f) Generate data from the process 
% X_t= cos(2 pi f_0t+?)+0.001 cos(2 pi f_1t+?),for n = 0,1,...,N?1, and 
% N = 128. Here, ? is defined as above, and ?t= 1s. Plot the periodogram 
% and the Hamming window periodogram, and comment on possible differences.

N = 128; % Number of sample
n = 0:1:N-1;
f_0 = 0.1;
f_1 = 0.4;
w_m = 0.54 - 0.46 * cos(2*pi*n/(N-1)); % Hamming window
U = 1/N * sum((w_m).^2);

Dt = 1;  % Difference in time steps
f_NQ = 1/(2*Dt);
freq = linspace(-f_NQ, f_NQ, N);
theta = unifrnd(0,2*pi,[1,1]);  % Uniform distributed variable

for i = 1:N
   X_t(i) = cos(2*pi*f_0*(i-1)+theta) + 0.001*cos(2*pi*f_1*(i-1)+theta); % Periodic process 
end

Ham = zeros(1,N);
for j = 1:N
   Ham(j) = X_t(j)*w_m(j); 
end

S_XX = Dt/N *abs( fft(X_t) ).^2;     % Periodogram estimate
S_XXshift = fftshift(S_XX);
dB_X = mag2db(S_XXshift);

S_Ham = Dt/(N*U)*abs( fft(Ham) ).^2; % Hamming window periodogram estimate
S_Hamshift = fftshift(S_Ham);
dB_H = mag2db(S_Hamshift);

plot(freq, dB_X)
xlabel('Frequencies  [s^{-1}]' )
ylabel('Periodogram estimate S_{XX}')
hold on
plot(freq, dB_H)
xlabel('Frequencies  [s^{-1}]')
legend('Periodogram','Hamming window')

%% 2g) Segment the dataset in non–overlapping blocks of length M = 32, use
% a data window on each segment, and compute the average of the resulting
% estimates. Plot the result of this operation, and compare with the 
% results above.

% Parameters
N = 256;                        % Number of samples
M = 32;                         % Window length
K = N/M;                        % Number of segments

Dt = 1;                         % Difference in time steps
f_NQ = 1/(2*Dt);
freq = linspace(-f_NQ, f_NQ, N);
theta = unifrnd(0,2*pi,[1,1]);  % Uniform distributed variable

sigma = 0.5;
w_t = normrnd(0,sigma^2,[1,N]); % white noise

A = 1/2;                        % Amplitude
f0 = 0.1;                       % Frequency [s^-1]
%
X_t = zeros(1,N);
for i = 1:N
    X_t(i) = A*cos(2*pi*f0*(i-1) + theta) + w_t(i); % Periodic process
end

% We begin by making a data segment number k of length M
X_k = zeros(K,M);

X_k(1,:) = X_t(1:M*1);

for k = 1:K
    X_k(k,:) = X_t(M*(k-1)+1:M*k);
end