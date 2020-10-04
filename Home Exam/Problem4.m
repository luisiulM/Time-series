%% 4
load('Riff.mat')

fs = 8192;  % Sampling frequency [Hz]
a = 0.5;    % Echoes strenght
e = 1;      % Number of echoes
tD = 0.15;  % Time delay

[RiffAndEcho,riff_D] = CreateRiffAndEcho(Riff,a,e,tD);

% Estimating the Power spectral densities
[PSD_R,f1] = fftPSDEstimate(riff_D,fs);
[PSD_RE,f2] = fftPSDEstimate(RiffAndEcho,fs);

% Estimating the difference between the PSD
Diff_PSD = PSD_R - PSD_RE;

figure;
plot(f1,PSD_R)
xlabel('Frequencies axis [Hz]')
ylabel('Estimated Power Spectral Density of Riff [W/Hz] ')
xlim([-4100 4100])

figure;
plot(f2,PSD_RE)
xlabel('Frequencies axis [Hz]')
ylabel('Estimated Power Spectral Density of Riff and Echo [W/Hz] ')
xlim([-4100 4100])

% It doesn't matter whether you use f1 or f2 in the frequency axis since
% both are the same.
figure;
plot(f1,Diff_PSD)
xlabel('Frequencies axis [Hz]')
ylabel('Estimated difference between Riff and RiffandEcho [W/Hz] ')
xlim([-4100 4100])