%% 1f)
% Load Data
load('TouchToneMatrixNoise1.mat')
load('TouchToneMatrixNoise2.mat')

% Count the number of rows for each matrix. P.S. Since both matrices have
% the same dimentions, we only need to consider one of them.
R = size(TouchToneMatrixNoise1,1);

% Number of samples. Just like before, the number of samples is the same
% for both matrices, so we only consider one.
N = size(TouchToneMatrixNoise1,2);

% Number of padd samples
pad = 700;

% Sampling frequency
fs = 8000;

% Creating matries where each row corresponds to the power spectral
% density of each corresponding row of the Touch Tone Matrix. 
S_X1 = zeros(R,N+pad);
S_X2 = zeros(R,N+pad);
for r = 1:R
    [S_X1(r,:),f1] = PSDEstimate(TouchToneMatrixNoise1(r,:),fs,pad);
    [S_X2(r,:),f2] = PSDEstimate(TouchToneMatrixNoise2(r,:),fs,pad);
end

% Ploting
for i = 1:R
    figure;
    plot(f1,S_X1(i,:))
    xlabel('Frequencies [Hz]')
    ylabel('PSD Estimate of TouchToneMatrixNoise1 [W/Hz]')
    
    figure;
    plot(f2,S_X2(i,:))
    xlabel('Frequencies [Hz]')
    ylabel('PSD Estimate of TouchToneMatrixNoise2 [W/Hz]')
end

%% 1g)

% Creating matries where each row corresponds to the power spectral
% density of each corresponding row of the Touch Tone Matrix. 
fftS_X1 = zeros(R,N+pad);
fftS_X2 = zeros(R,N+pad);
for r = 1:R
    fftS_X1(r,:) = fftPSDEstimate(TouchToneMatrixNoise1(r,:),fs,pad);
    fftS_X2(r,:) = fftPSDEstimate(TouchToneMatrixNoise2(r,:),fs,pad);
end

% Ploting
for i = 1:R
    figure;
    plot(f1,S_X1(i,:))
    hold on
    plot(f1,fftS_X1(i,:))
    xlabel('Frequencies [Hz]')
    ylabel('PSD Estimate of TouchToneMatrixNoise1 [W/Hz]')
    legend('PSDEstimate','fftPSDEstimate')
    
    figure;
    plot(f2,S_X2(i,:))
    hold on
    plot(f2,fftS_X2(i,:))
    xlabel('Frequencies [Hz]')
    ylabel('PSD Estimate of TouchToneMatrixNoise2 [W/Hz]')
    legend('PSDEstimate','fftPSDEstimate')
end

%% 1h)
% Window length
window = 400;

% Sampling frequency
fs = 8000;

% Number of overlapping samples
overlap = 100;

% Length of the FFT
Nfft = window + 50;

% Creating matries where each row corresponds to the power spectral
% density of each corresponding row of the Touch Tone Matrix. 
PSD1 = zeros(R,Nfft);
PSD2 = zeros(R,Nfft);
for r = 1:R
    [PSD1(r,:),f1] = WOSA(TouchToneMatrixNoise1(r,:),window,fs,overlap,Nfft);
    [PSD2(r,:),f2] = WOSA(TouchToneMatrixNoise2(r,:),window,fs,overlap,Nfft);
end

% Ploting
for i = 1:R
    figure;
    plot(f1,PSD1(i,:))
    xlabel('Frequencies [Hz]')
    ylabel('WOSA Estimate of TouchToneMatrixNoise1 [W/Hz]')
    
    figure;
    plot(f2,PSD2(i,:))
    xlabel('Frequencies [Hz]')
    ylabel('WOSA Estimate of TouchToneMatrixNoise2 [W/Hz]')
end