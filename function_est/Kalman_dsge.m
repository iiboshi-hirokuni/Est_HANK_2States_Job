function [rloglh,obssmooth,DD] = Kalman_dsge(para,yy,nshock,ZZ,T1,TC,TEPS)

% This procedure evaluates the likelihood function of the 
% monetary DSGE model
% retcode = -1 : non existence
%         = 0  : existence and uniqueness
%         = 1  : existence and non-uniqueness

% Parameters
%  h,  sigma_c, sigma_L, beta, phi, tau, Rk, gam_p, 
%  gam_w, psi_p, psi_w, alpha, psi, k_y, g_y,
%  rho_m, mu_pi, mu_y
%  e_c, e_inv, e_q, e_L, e_w, e_z, e_p, e_g, e_m

% npara   = size(para, 1);
% retcode = 0;

% solve the DSGE model
% 
[T1,TC,TEPS,RC] = dsgesolv_SW(para);
%DD = zeros(size(ZZ, 1),1);

 TT = T1;
 RR = TEPS;

nseries  = size(yy, 2);
nstate   = size(T1, 2);

nobs      = size(yy, 1);
loglh     = 0;
loglhzero = -1E8;
obsmean   = zeros(nobs, 2*nstate);
obssmooth =zeros(nobs, 2*nstate);

obsvar    = zeros(nobs, nseries);
nu_save   = zeros(nseries,nobs);
Ft_save   = zeros(nseries,nseries*nobs); 
% Kg_save   = zeros(nseries,nseries*nobs);
at_save   = zeros(nobs,2*nstate);
pt_save   = zeros(2*nstate,2*nstate*nobs);
% shock     = zeros(nobs, nshock);

% create system matrices for state space model

% These matrices are regime independent


z_star_bar = para(14,1);
psi_bar    = para(15,1);

l_bar   = mean(yy(:,5));  %1.0;
r_n_bar = mean(yy(:,8));
pi_bar  = 0.25;

DD = zeros(nseries,1);
DD(1, 1) = z_star_bar;
DD(2, 1) = z_star_bar;
DD(3, 1) = z_star_bar+psi_bar;
DD(4, 1) = z_star_bar;
DD(5, 1) = l_bar;
DD(6, 1) = pi_bar;
DD(7, 1) = -psi_bar;
DD(8, 1) = r_n_bar;
% DD(9, 1) = 0;
% DD(10, 1) = z_star_bar;

% HH = zeros(nseries,nseries);
HH = 0.01*eye(nseries,nseries);
QQ = createcov(para(31:39,1));
VV = zeros(nshock,nseries);

% Check whether covariance matrix QQ is positive definite

if sum(eig(QQ) < 0) > 0
   loglh = loglhzero;
   return;
end

% We can now define the initial mean and variance for the state vector
%
%At = zeros(nstate,1);
At = zeros(nstate*2,1);


Pt = dlyap(TT,RR*QQ*RR');
Pt = [Pt zeros(nstate);
      zeros(nstate) Pt];
Pt = 0*Pt;  

TT = [TT zeros(nstate);
      diag(ones(nstate,1)) zeros(nstate)];

RR = [RR; zeros(nstate,nshock)];

% compute likelihood with Kalman filter

t = 1;
while t <= nobs
   
   At1 = At;
   Pt1 = Pt;
   
   % Forecasting
   alphahat = TT*At1;
   Phat = TT*Pt1*TT' + RR*QQ*RR';
   yhat = ZZ*alphahat + DD;
   nu   = yy(t,:) - yhat';
   
   Ft   = ZZ*Phat*ZZ' + HH + ZZ*RR*VV + (ZZ*RR*VV)';
   Ft   = 0.5*(Ft + Ft');
   
   K_g = TT*Phat*ZZ'*inv(Ft);   % Kalman Gain
  
   at_save(t,:) = alphahat;
   pt_save(:,(t-1)*2*nstate+1:t*2*nstate) = Phat;
   
   loglh = loglh -0.5*size(yy, 2)*log(2*pi)-0.5*log(det(Ft)) ...
           - 0.5*nu*inv(Ft)*nu';
   
   % Updating
   At = alphahat + (Phat*ZZ' + RR*VV)*inv(Ft)*nu';
   Pt = Phat - (Phat*ZZ'+RR*VV)*inv(Ft)*(Phat*ZZ'+RR*VV)';
   
   %  store
   obsmean(t,:) =alphahat';
   obsvar(t,:)  = diag(Ft)';
   nu_save(:,t) = nu';
   Ft_save(:,(t-1)*nseries+1:t*nseries) = Ft;
   Kg_save(:,(t-1)*nseries+1:t*nseries) = K_g;
   
   t = t+1;
end  

% state smoothing by Durbin & Koopman (2001, p73-75)

r_t = zeros(nstate*2,1);
N_t = zeros(2*nstate,2*nstate); 
eta = zeros(nshock,nobs);

t = nobs; 
while t > 0
  nu = nu_save(:,t);
  Ft =  Ft_save(:,(t-1)*nseries+1:t*nseries);
  K_g = Kg_save(:,(t-1)*nseries+1:t*nseries);
  Phat = pt_save(:,(t-1)*2*nstate+1:t*2*nstate);
  L_t = TT- K_g*ZZ; 
  
  r_t = ZZ'*inv(Ft)*nu + L_t'*r_t;
  N_t = ZZ'*inv(Ft)*ZZ + L_t'*N_t*L_t;
  alpha_t = at_save(t,:)' + Phat * r_t;     
  V_t   = Phat - Phat*N_t*Phat; 

  yhat  = alpha_t;
  obssmooth(t,:) = yhat';

  eta(:,t) = QQ*RR'*r_t;
  
  t = t-1;
end

% shock = eta';

rloglh = real(loglh);

%-----------------------------
function [omega] = createcov(para)

npara = max(size(para));
omega = zeros(npara, npara);
for i = 1:npara
  omega(i, i) = para(i)^2;
end
