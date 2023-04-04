function [st_save, smooth1, smooth2,rloglh] = decomp_HANK(g1,impact,T,N,method,yy, ZZ,m)


vtime = linspace(0,T,N);
dt = vtime(2)-vtime(1);

% Preallocation
[nvars,~] = size(g1);
values = zeros(nvars,N);

if (method == 'implicit')
    % if N is small, might be faster just to do backslash instead. To use
    %    backslash method, just uncomment line 51/54 and comment line 50/53
    gg1 = inv(speye(size(g1)) - g1*dt);
    % gg1 = speye(size(g1))-g1*dt;

elseif (method == 'explicit')
    gg1 = speye(size(g1))+g1*dt;

end

%% start Kalman filter

 TT = gg1;
 RR = (dt^(1/2))*gg1*impact; 
 
nseries  = size(yy, 2);
nstate   = size(TT, 2);
nshock   = size(impact,2);

nobs      = size(yy, 1);
loglh     = 0;
loglhzero = -1E8;
obsmean   = zeros(nobs, m*nstate);
obssmooth =zeros(nobs, m*nstate);

obsvar    = zeros(nobs, nseries);
nu_save   = zeros(nseries,nobs);
Ft_save   = zeros(nseries,nseries*nobs); 
% Kg_save   = zeros(nseries,nseries*nobs);
at_save   = zeros(nobs,m*nstate);
pt_save   = zeros(m*nstate,m*nstate*nobs);
% shock     = zeros(nobs, nshock);

% obsmean =zeros(nobs, nseries);
smooth1 =zeros(nobs, nseries);
smooth2 =zeros(nobs, nseries);
% obsvar    = zeros(nobs, nseries);
% nu_save   = zeros(nseries,nobs);
% ft_save   = zeros(nseries,nseries*nobs); 
% kg_save   = zeros(nstate,nseries*nobs);
% at_save   = zeros(nobs,2*nstate);
% pt_save   = zeros(2*nstate,2*nstate*nobs);
% shock     = zeros(nobs, nshock);
st_save = zeros(nobs,nseries*nshock);


% create system matrices for state space model
% These matrices are regime independent

DD = zeros(nseries,1);

% HH = zeros(nseries,nseries);
HH = 0.01*eye(nseries,nseries);
QQ = eye(nshock);
VV = zeros(nshock,nseries);

% Check whether covariance matrix QQ is positive definite

if sum(eig(QQ) < 0) > 0
   loglh = loglhzero;
   return;
end

% We can now define the initial mean and variance for the state vector
At = zeros(nstate*m,1);
a1 = At;

% Pt = dlyap(TT,RR*QQ*RR');
Pt = 0.1*eye(nstate);


if m==2
  Pt = [Pt zeros(nstate);
      zeros(nstate) Pt];
  TT = [TT zeros(nstate);
      diag(ones(nstate,1)) zeros(nstate)];
  RR = [RR; zeros(nstate,nshock)];
end

p1 = Pt; 

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
   pt_save(:,(t-1)*m*nstate+1:t*m*nstate) = Phat;
   
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


rloglh = real(loglh);


% switch smooth_type
%     case 'state_smooth'
% state smoothing by Durbin & Koopman (2001, p73-75)
%                 by Durbin & Koopman (2012, p91)

r_t = zeros(nstate*2,1);
N_t = zeros(2*nstate,2*nstate); 

 for t = nobs:-1:1 

  nu = nu_save(:,t);
  Ft =  Ft_save(:,(t-1)*nseries+1:t*nseries);
  K_g = Kg_save(:,(t-1)*nseries+1:t*nseries);
  Phat = pt_save(:,(t-1)*2*nstate+1:t*2*nstate);
  L_t = TT- K_g*ZZ; 
  
  r_t = ZZ'*inv(Ft)*nu + L_t'*r_t;     % Eq.(4.32)
  N_t = ZZ'*inv(Ft)*ZZ + L_t'*N_t*L_t;
  alpha_t = at_save(t,:)' + Phat * r_t;     
  V_t   = Phat - Phat*N_t*Phat; 
  
   smooth1(t,:) =  (ZZ*alpha_t + DD)';  
  
 end


%     case 'disturbance_smooth' 
% disturbance smoothing by Durbin & Koopman (2001, p??)
%                       by Durbin & Koopman (2012, p96)

r_t = zeros(nstate*m,1);
N_t = zeros(m*nstate,m*nstate); 
eta = zeros(nshock,nobs);

for t = nobs:-1:1 

  nu = nu_save(:,t);
  Ft =  Ft_save(:,(t-1)*nseries+1:t*nseries);
  K_g = Kg_save(:,(t-1)*nseries+1:t*nseries);
  Phat = pt_save(:,(t-1)*2*nstate+1:t*2*nstate);
  L_t = TT- K_g*ZZ; 
  eta(:,t) = QQ*RR'*r_t;
  
  r_t = ZZ'*inv(Ft)*nu + L_t'*r_t;     % Eq.(4.32)
  N_t = ZZ'*inv(Ft)*ZZ + L_t'*N_t*L_t; 
  
end

st1 = a1;  %zeros(80,1);
 
 alpha_t = st1 + p1*r_t;  % initialize
%   alpha_t = st1;

for t = 1:nobs
   smooth2(t,:) = (ZZ * alpha_t(1:m*nstate))' + DD';
  alpha_t = TT * alpha_t + RR * eta(:,t);  % Eq.(4.85)
end


% /*   historical decomposition
% */

nvar    = nseries;

for  index = 1:nshock
%    st1 = zeros(40,1);
   r1 = zeros(11,1);
   r1(index) = r_t(index);  
   
%    st1 = a1 + p1*r_t;  %QQ(index,index)/sum(diag(QQ)) ;  % initialize    
    st1 = a1;                                             % initialize   

 for i =1:nobs-1   
      st_save(i,1+(index-1)*nvar:index*nvar ) = (ZZ*st1)';
    
     shock_tmp = zeros(nshock,1);
     shock_tmp(index) = eta(index,i);     
     st = TT*st1 + RR*shock_tmp; 
     st1 = st;
 end
 
end

% end %( end switch)



