%% 2d)
% load data
load('Cardiovascular.mat')

% Number of samples.
N = length(Cardiovascular);

n = 0:1:N-1;

% Plotting 
plot(n,Cardiovascular)
xlim([0 N])
xlabel('Weekly samples n')
ylabel('Average weekly cardiovascular mortality Y_n')

%% 2e)

% load data
load('Cardiovascular.mat')

% Number of samples.
N = length(Cardiovascular);

% Creating the lag h vector 
lag = 0:1:N-1;

% Estimating the de-trended autocovariance of the Cardiovascular data
Decov = De_trendcov(Cardiovascular);

% plotting
stem(lag,Decov,'.')
xlim([-10 510])
xlabel('Lag h')
ylabel('Autocovariance of cardiovascular mortality Y_n')