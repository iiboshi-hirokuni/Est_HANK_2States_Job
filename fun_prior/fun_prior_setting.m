
function  Pr = fun_prior_setting(npara)

%%  prior  setting

datapath = './data/';

% **         0 is no prior             
% **         1 is BETA(mean,stdd)
% **         2 is GAMMA(mean,stdd)
% **         3 is NORMAL(mean,stdd)
% **         4 is INVGAMMA(s^2,nu)
% **         5 is uniform(a,b)

%  prior setting file
    datafilename = 'prior_setting_HANK';   

prior_set = csvread(strcat(datapath, [ datafilename '.csv' ] ), 1, 2); %, [1,2,23,7]);

Pr.pmean = prior_set(1:npara,1);
Pr.pstdd = prior_set(1:npara,2);
Pr.pmask = prior_set(1:npara,3);
Pr.pshape = prior_set(1:npara,4);
Pr.upper = prior_set(1:npara,5);
Pr.lower = prior_set(1:npara,6);


%    P.pmean =[1.5, 0.995, 1.2, 0, 0, 1, 6, 0.024, 1,  ];             
%    P.pstdd  = [0.05, 0.01, 0.05 ];  % standard deviation     
%    P.pmask =[0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
%    P.pshape =[3,  3,   3  ];     

