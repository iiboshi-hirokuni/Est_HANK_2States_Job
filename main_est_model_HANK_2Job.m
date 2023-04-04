% 
%  Sequential Monte Carlo Methods with Particle Filter
%  Canonical New Keynesian Model (Rotemberg Pricing)    
%   -Imposes the zero lower bound on the interest rate
%

clear all
clc

 os = 'windows';
% os = 'mac';

switch os
%  addpath(genpath('/Users/hirokuniiiboshi/Dropbox/matlab_code/Phact'));
case 'mac'
  addpath(genpath('/Users/hirokuniiiboshi/Dropbox/matlab_code/DSGE_HANK/estimation'));
  addpath(genpath('/Users/hirokuniiiboshi/Dropbox/matlab_code/DSGE_HANK/estimation/AutoDiff/@myAD'));
case  'windows'
  addpath(genpath('./AutoDiff'));
  addpath(genpath('./AutoDiff/@myAD'));
end  
  addpath(genpath('./function_phact'));
  addpath(genpath('./function_solv'));
  addpath(genpath('./fun_prior'));
  addpath(genpath('./fun_smc'));
  addpath(genpath('./function_est'));
  
% Turn off the warning message from auto diff
% You should read the warning once, but turn it off after reading it.
%% 
   warning('off','AutoDiff:maxmin')
   addpath(genpath('./function_est'));
   addpath(genpath('./data'));
   %% 

tstart = tic;                           % Job timer start

disp('Start SMC^2 ')
ncores  =  4  % number of core of CPU for parallel computing 

data_country = 1  % 1: Japan, 2:US
def_switch   = 1  % 1st deference for GDP = 1, level = 0

%% setting 
nsim       = ncores*20   % # of particles of parameters
nstage     = 2           % # of stages
npara      = 18;         % # of parameters
% cc1      = 0.035 ;     % adjustment coefficient of MH

%% parallel 
 
  delete(gcp('nocreate')) %   
  parpool('local', ncores) 
%  parpool( ncores)


main_smc_2job



%% end



