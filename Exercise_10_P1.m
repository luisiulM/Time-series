%% 1a) Generate N = 51 samples of a white Gaussian process, with mean value
% equal to zero. Let the variance be sd_w = 1 and let t = 1.

N = 50;          % Number of samples
n = 0:1:N-1;
Dt = 1;          % Difference in time steps
Df = 1/(N*Dt);
f_NQ = 1/(2*Dt); % Nyquist frequency
sd_w = 1;        % variance of white noise 
mu = 0;
X_n = normrnd(mu,sd_w,[1,N]);

%ArMd = arima('Constant',0,'AR',{0} ,'MA' ,{0} ,'Variance',(sd_w)); % white noise model 
%X_n = simulate(ArMd,N); % Simulating white Gaussian process

%% 1b) Plot the result from (a) on a decibel scale (dB) as a function of fre-
% quency (two-sided). Compare with the true power spectral density.

freq = linspace(-f_NQ, f_NQ, N);

S_XX = zeros(1,N);
for k = 1:N
    
   S_XX(k) = Dt/N*abs( sum(X_n.*exp(-1i*2*pi*(k-1)*Df*n*Dt)) )^2;
    
end

dB_S = mag2db(S_XX);

figure;
plot(freq,dB_S)

%% 1c) Let N = 101 and repeat (a) and (b). Comment on the variance prop-
% erties of the estimator.

N2 = 101;
n2 = 0:1:N2-1;
Df2 = 1/(N2*Dt);
X_n2 = normrnd(mu,sd_w,[1,N2]);

freq2 = linspace(-f_NQ, f_NQ, N2);

S_XX2 = zeros(1,N2);
for k = 1:N2
    
   S_XX2(k) = Dt/N2*abs( sum(X_n2.*exp(-1i*2*pi*(k-1)*Df2*n2*Dt)) )^2;
    
end

dB_S2 = mag2db(S_XX2);

figure;
plot(freq2,dB_S2)

% Comment: 

%% 1d) Let N = 401 and repeat (a) and (b). Comment on the result.

N3 = 101;
n3 = 0:1:N3-1;
Df3 = 1/(N3*Dt);
X_n3 = normrnd(mu,sd_w,[1,N3]);

freq3 = linspace(-f_NQ, f_NQ, N3);

S_XX3 = zeros(1,N3);
for k = 1:N3
  
   S_XX3(k) = Dt/N3*abs( sum(X_n3.*exp(-1i*2*pi*(k-1)*Df3*n3*Dt)) )^2;
    
end

dB_S3 = mag2db(S_XX3);

figure;
plot(freq3,dB_S3)

%% 1e) Repeat (a) and (b), but now using fft. Check that the results are 
% the same.

FFT_S = Dt/N*abs( fft(X_n) ).^2;

dB_FS = mag2db(FFT_S);

figure;
plot(freq,dB_FS)
hold on
plot(freq,dB_S)
xlabel('Frequencies')
ylabel('dB')

%% 1f) Repeat (c), but now by using fft.

FFT_S2 = Dt/N2*abs( fft(X_n2) ).^2;

dB_FS2 = mag2db(FFT_S2);

figure;
plot(freq2,dB_FS2)
hold on
plot(freq2,dB_S2)
xlabel('Frequencies')
ylabel('dB')

%% 1g) Repeat (d), but now by using fft.

FFT_S3 = Dt/N3*abs( fft(X_n3) ).^2;

dB_FS3 = mag2db(FFT_S3);

figure;
plot(freq3,dB_FS3)
hold on
plot(freq3,dB_S3)
xlabel('Frequencies')
ylabel('dB')

%% 1h) Consider once more the dataset from (a), i.e. N = 51. Use the zero-
% padding fft(x, 512), and plot the resulting periodogram estimate. Com-
% pare with the results from (a) and (e). Comment on the result.

freq_pad = linspace(-f_NQ, f_NQ, 512);

FFT_pad = Dt/N*abs( fft(X_n,512) ).^2;

dB_FSpad = mag2db(FFT_pad);

figure;
plot(freq,dB_FS)
hold on
plot(freq,dB_S,'o--')
xlabel('Frequencies')
ylabel('dB')
plot(freq_pad,dB_FSpad, '--')
legend('S_{XX}','S_{XX} with FFT','S_{XX} paded with FFT')
%% 1i) Evaluate the periodogram estimate of the process on the Fourier-
% frequencies. Use a valid frequency axis. Plot the result both on a linear
% and dB scale. Comment on the result.

N = 40;                        % Number of samples
Dt = 0.2;                      % Difference in time steps
f_NQ = 1/(2*Dt);
freq = linspace(-f_NQ, f_NQ, N);
theta = unifrnd(0,2*pi,[1,1]);  % Uniform distributed variable

A = 1;                         % Amplitude
f0 = 0.6;                      % Frequency [s^-1]

X_t = zeros(1,N);
for n = 1:N
    X_t(n) = A*cos(2*pi*f0*(n-1) + theta); % Periodic process
end

S_XX = Dt/N*abs( fft(X_t) ).^2; % Periodogram estimate

S_XXshift = fftshift(S_XX);

dB_FS = mag2db(S_XXshift);          % dB scale

plot(freq, S_XXshift)
xlabel('Frequencies  [s^{-1}]' )
ylabel('Periodogram estimate S_{XX}')
hold on
plot(freq, dB_FS)
legend('normal','dB scale of S_{XX}')

%% 1j)

N = [40 , 75 , 130];           % Number of samples
Dt = 0.2;                      % Difference in time steps
f_NQ = 1/(2*Dt);

freq = linspace(-f_NQ, f_NQ, N(1));
freq2 = linspace(-f_NQ, f_NQ, N(2));
freq3 = linspace(-f_NQ, f_NQ, N(3));

theta = unifrnd(0,2*pi,[1,1]);      % Uniform distributed variable
sigma = [2 , 5 , 7];

w_t = normrnd(0,sigma(1),[1,N(1)]);  % white noise
w_t2 = normrnd(0,sigma(2),[1,N(2)]); % white noise
w_t3 = normrnd(0,sigma(3),[1,N(3)]); % white noise

A = 1;                         % Amplitude
f0 = 0.6;                      % Frequency [s^-1]

X_t = zeros(1,N(1));
X_t2 = zeros(1,N(2));
X_t3 = zeros(1,N(3));
for n = 1:N(1)
    X_t(n) = A*cos(2*pi*f0*(n-1) + theta) + w_t(n); % Periodic process
end
for n = 1:N(2)
    X_t2(n) = A*cos(2*pi*f0*(n-1) + theta) + w_t2(n); % Periodic process
end
for n = 1:N(3)
    X_t3(n) = A*cos(2*pi*f0*(n-1) + theta) + w_t3(n); % Periodic process
end

S_XX = Dt/N(1)*abs( fft(X_t) ).^2; % Periodogram estimate
S_XX2 = Dt/N(2)*abs( fft(X_t2) ).^2;
S_XX3 = Dt/N(3)*abs( fft(X_t3) ).^2;

S_XXshift = fftshift(S_XX);
S_XXshift2 = fftshift(S_XX2);
S_XXshift3 = fftshift(S_XX3);

dB_FS = mag2db(S_XXshift);          % dB scale
dB_FS2 = mag2db(S_XXshift2);
dB_FS3 = mag2db(S_XXshift3);

subplot(2,3,1)
plot(freq, S_XXshift)
xlabel('Frequencies  [s^{-1}]' )
ylabel('Periodogram estimate S_{XX}')
subplot(2,3,4)
plot(freq, dB_FS)
xlabel('Frequencies  [s^{-1}]')
ylabel('dB scale of S_{XX}')

subplot(2,3,2)
plot(freq2, S_XXshift2)
xlabel('Frequencies  [s^{-1}]' )
ylabel('Periodogram estimate S_{XX}')
subplot(2,3,5)
plot(freq2, dB_FS2)
xlabel('Frequencies  [s^{-1}]')
ylabel('dB scale of S_{XX}')

subplot(2,3,3)
plot(freq3, S_XXshift3)
xlabel('Frequencies  [s^{-1}]' )
ylabel('Periodogram estimate S_{XX}')
subplot(2,3,6)
plot(freq3, dB_FS3)
xlabel('Frequencies  [s^{-1}]')
ylabel('dB scale of S_{XX}')
