%% data
datapath = './data/';

%% Data (year on year)
datafilename = 'data_us';         % actual data

yy1 = csvread(strcat(datapath, [ datafilename '.csv' ] ), 1, 1);

yy = [ yy1(:,1) yy1(:,2) yy1(:,3)]; % dy pi r

ZZ = [100; 100; 100]; % coefficients between observables and endogenous variables of the measurement equation



%%
Tobs = size(yy,1);