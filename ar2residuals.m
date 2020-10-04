% load and plot data

addpath /Users/hfr003/Downloads
data = dlmread('rec.txt');

plot(data)

% fit ar2 model:
ToEstmdl = arima(2,0,0);
EstMdl=estimate(ToEstmdl,data)


%% compute residuals:

[E,V] = infer(EstMdl,data);

%% plot residuals:

plot(E)

%% plot ACF:

autocorr(E)