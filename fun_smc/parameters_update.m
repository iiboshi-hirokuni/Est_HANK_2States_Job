
% Preferences
coefrra = para(1);
frisch  = para(2);
adjfricshgridfrac = para(3);
meanlabeff = 3; 	%so that at h=1/3 outout will be approximately = 1;
labdisutil = meanlabeff ./ ( (0.75 .^(-coefrra)) * ((1.0/3.0).^(1.0/frisch))); %guess labdisutil so that at average wages and average consumption hours =1/3 (sets C/Y = 0.75);
maxhours = 1;

% Production
ceselast = para(4);		 % elasticity of substitution / demand
priceadjust = para(5);

% Policy parameters
taylor_inflation = para(6);	%taylor rule coefficient on inflation
taylor_outputgap = para(7);		%taylor rule coefficient on output
labtax = para(8);			% marginal tax rate
govbondtarget = para(9);		%multiple of quarterly GDP
lumptransferpc = para(10);	 %6% of quarterly GDP in steady state
govbcrule_fixnomB = para(11);


% Aggregate shocks
ssigma_MP = para(12)/10;
ttheta_MP = para(13);
ssigma_FP = para(14)/10;
ttheta_FP = para(15);
ssigma_TFP = para(16)/10;
ttheta_TFP = para(17);

% transition prob
%ymarkov_combined = [-0.5    0.5;  0.035    -0.035]; 
ymarkov_combined(2,1)=para(18);
ymarkov_combined(2,2)=-1*para(18);

