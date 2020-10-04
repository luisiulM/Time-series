%% 3a)
% Load Music.mat Data
load('Music.mat')

% Defining the order of Autoregressive random process
p = [3 5 11 17];

N = length(Music);
P = length(p);

% Estimating the prediction of the autoregressive process of order p, from
% the Music Data.
ARpred1 = zeros(P,N);
for i = 1:P
    ARpred1(i,:) = ARfitting(Music,p(i));
    figure;
    plot(ARpred1(i,:))
    xlim([-1000 N])
end

%% 3b)

% Here we estimate the prediction of the previous prediction and take the
% difference, in order to see the differences.
ARpred2 = zeros(P,N);
for i = 1:P
    ARpred2(i,:) = ARfitting(ARpred1(i,:),p(i));
    DiffARpred = ARpred2(i,:) - ARpred1(i,:);
    figure;
    plot(DiffARpred)
    xlim([-1000 N])
    ylim([-0.5 0.5])
end