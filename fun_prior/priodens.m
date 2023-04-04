function lnprior = priodens(para, pmean, pstdd, pshape)
% /* This procedure computes a prior density for
% ** the structural parameters of the DSGE models
% ** pshape: 0 is point mass, both para and pstdd are ignored
% **         1 is BETA(mean,stdd)
% **         2 is GAMMA(mean,stdd)
% **         3 is NORMAL(mean,stdd)
% **         4 is INVGAMMA(s^2,nu)
% **         5 is uniform(a,b)
% */
% local a, b, lnprior, prioinfo, nprio, i;
lnprior = 0;
a = 0;
b = 0;

[nprio,cprio] = size(pshape);
prioinfo = [zeros(nprio,2),pshape];

% i = 1;
% do until i > nprio;
for i = 1:nprio
   if prioinfo(i,3) == 1 %/* BETA Prior */
      a = (1-pmean(i))*pmean(i)^2/pstdd(i)^2 - pmean(i);
      b = a*(1/pmean(i) - 1);
      lnprior = lnprior + lpdfbeta(para(i),a,b);   
   elseif prioinfo(i,3) == 2 %/* GAMMA PRIOR */
      b = pstdd(i)^2/pmean(i);
      a = pmean(i)/b;
      lnprior = lnprior + lpdfgam(para(i),a,b);
   elseif prioinfo(i,3) == 3 %/* GAUSSIAN PRIOR */
      a = pmean(i);
      b = pstdd(i);
      lnprior = lnprior + lpdfnor(para(i),a,b);
   elseif prioinfo(i,3) == 4 %/* INVGAMMA PRIOR */
      a = pmean(i);  % variance 
      b = pstdd(i);  % degree of freedom
      lnprior = lnprior + lpdfig(para(i),a,b);
   elseif prioinfo(i,3) == 5 %/* Uniform Prior */
      a = pmean(i);  % lower band
      b = pstdd(i);  % upper band
      lnprior = lnprior + log(1/(b-a));
   end
      prioinfo(i,1) = a;
      prioinfo(i,2) = b;
end

%-----------------------------
function ldens = lpdfbeta(x,a,b)
% /* log BETA PDF
% */
% local ldens;
  ldens = lngam(a+b) - lngam(a) - lngam(b) + (a-1)*log(x) + (b-1)*log(1-x);
  
%-----------------------------
function ldens = lpdfgam(x,a,b)
% /* log GAMMA PDF
% */
% local ldens;
  ldens = -lngam(a) -a*log(b)+ (a-1)*log(x) -x/b ;

%-----------------------------
function ldens = lpdfig(x,a,b)
% /* log INVERSE GAMMA
% */
% local ldens;
   
   ldens = log(2) - gammaln(b/2) + (b/2)*log(b*a^2/2) - ( (b+1)/2 )*log(x^2) - b*a^2/(2*x^2);

%-----------------------------
function ldens = lpdfnor(x,a,b)
% /* log NORMAL PDF
% */
% local ldens;
  ldens = -0.5*log(2*pi) - log(b) - 0.5*(x-a)^2/b^2;

%-----------------------------
function lng = lngam(x)
% /* This procedure computes the log-gamma function
% */
% local k,lng;
%   if x <= 3;
%      lng= log(gamma(x));
%   else
%      k = floor(x-2);
%      lng= sum(log(x-1:-1:k)) + log(gamma(x-k));
%   end
lng = gammaln(x);

%-----------------------------


