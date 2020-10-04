%% 1 f)

n = 100;                                                   % Number of simulations
ArMd = arima('Constant',0,'AR',{0.9},'Variance',(1/1900)); % AR(1) with phi 0.9


x_t = simulate(ArMd,n);      % Simulate 100 observations from ArMd.
x_t2 = simulate(ArMd,n*10);  % Simulate 1000 observations from ArMd.
x_t3 = simulate(ArMd,n*100); % Simulate 10000 observations from ArMd.

%%
subplot(3,1,1)
autocorr(x_t)
subplot(3,1,2)
autocorr(x_t2)
subplot(3,1,3)
autocorr(x_t3)