
%---------------------------------------------------------------------
% Path Specification
%---------------------------------------------------------------------
priopath = './prior/';
postpath = './posterior/';
datapath = './data/'; %'./data/';
resupath = './results/';

if data_country == 1
%     datafilename = 'data_JP.csv';         % actual data
    datafilename = 'data_jpn_2020.csv';  
    
    obs_ratio = zeros(3,1);
    obs_ratio(1)= 100; % output
    obs_ratio(2)= 100; % inflation
    obs_ratio(3)= 100; % interest rate
    
    data = csvread(strcat(datapath, datafilename), 1, 1);
%      series_YT = [ data(:,1) data(:,2) data(:,4) ];
    series_YT = [ data(:,1) data(:,2) data(:,3) ];
    
else
    datafilename = 'data_us.csv';         % actual data
    
    obs_ratio = zeros(3,1);
    obs_ratio(1)= 100; % output
    obs_ratio(2)= 100; % inflation
    obs_ratio(3)= 100; % interest rate
    
    series_YT = csvread(strcat(datapath, datafilename), 1, 1);
 
end    
    
    
    
%---------------------------------------------------------------------
% loading DATA
%---------------------------------------------------------------------


nobs = size(series_YT, 1);  % number of observations 
yy_m = mean(series_YT(:,1:3), 1);  % mean

YY = series_YT(:,1:3);




