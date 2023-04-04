nobs = size(series_YT,1);
ti = 1981:0.25:1981+(nobs-1)/4; 

figure(100)
subplot(3,1,1)
plot(ti, series_YT(:,1));title('Output Gap');
subplot(3,1,2)
plot(ti, series_YT(:,2)); title('Inflation');
subplot(3,1,3)
plot(ti, series_YT(:,3)); title('Nominal Rate');


