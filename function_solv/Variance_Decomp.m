
%--------------------------
% Variance Decompositions
%--------------------------

% 
% AA = eval(AAj);
% 
% sig_chol = eye(nshock);  % nshock = number of shocks: 3 or 4 
nshock =3;
nvar   = 3;
nirf   = 40;

nhorizon = 40;

%  name of stractural shocks
titlestr = {'MP shock: ','Gov shock: ', 'TFP shock: ',  'Price shock: '};

%  name of endogenous variables
ystr = {'Y_t', '\pi_t', 'r_t' };

 yyirf  = zeros(nirf,nvar,nshock);
 dyyirf  = zeros(nirf,nvar,nshock);


%     for sh_ind = 1:nshock
%         impact = sig_chol(:,sh_ind);
%        
%         
%         % impact of shock on  endogenous variables
%         ss = T0*impact;   
%         
%         s  = [ss;zeros(size(ss,1),1)];
%         dyyirf(1,:,sh_ind ) = AA*s;
%         yyirf(1,:,sh_ind ) = dyyirf(1,:,sh_ind );
%         
%         for t = 2:nhorizon %  # of horizons
%             
%             % dynamics of endogenous variables after shocks 
%             ss1 = T1*ss;  % ss1 = s_t+1, ss = s_t
%             s = [ss1;ss];
%             
%             % change from endogenous variables to observable variables
%             dyyirf(t,:,sh_ind ) = (AA*s)';
%             
%             % change from growth data to level data for output
%             yyirf(t,1,sh_ind) = yyirf(t-1,1,sh_ind ) + dyyirf(t,1,sh_ind ); 
%             % % percentage of inflation and interest rate (no change)
%             yyirf(t,2:3,sh_ind) = dyyirf(t,2:3,sh_ind );            
%             
%             ss = ss1; % proceed to t= t+1
%         end        
%        
%     end

T = 40;
N = 40;
dt = T/N;
yyirf= zeros(N,3,3);

%% MP shocks
for i= 1:3
impact1=impact(:,i);

vAggregateShock	= zeros(1,N);
vAggregateShock(:,1) = 1/sqrt(dt);
trans_mat = inv_state_red*from_spline;
[simulated,vTime] = simulate(G1,impact1,T,N,vAggregateShock,'implicit',trans_mat,[n_v,n_v+n_g-2:n_v+n_g+6]);

[simulated_gg0,~] = simulate(G1,impact1,T,N,vAggregateShock,'implicit',trans_mat,[n_v+1:n_v+n_g-2]);
simulated_gg = reshape(simulated_gg0(:,2),I,J);
simulated_gg(I,:)=zeros(2,1);

yyirf(:,2,i) = simulated(1,:)';
% monetary_shock = simulated(2,:)';
% FP_shock = simulated(3,:)';
% consumption = (simulated(2+5,:)')/vars_SS(n_v+n_g+3);
yyirf(:,1,i) = (simulated(2+6,:)')/vars_SS(n_v+n_g+4);
% lab_sup = (simulated(2+4,:)')/vars_SS(n_v+n_g+2);
% wage = simulated(2+3,:)'/vars_SS(n_v+n_g+1);
yyirf(:,3,i) = simulated(2+8,:)';
% bond     = simulated(2+7,:)'/vars_SS(n_v+n_g+5);
end

  
    ratio = zeros(nshock,1);
    
         disp(' ')
         disp('Variance Decomposition' )
         disp('===========================================')      
         
    for i = 1:3 % endogenous variables
        disp([ 'variable: ', char(ystr(i)) ])
        disp([ '    ', char(titlestr(1)),char(titlestr(2)),char(titlestr(3) ) ])
        disp('--------------------------------------------')   
        
       for t = 1:nhorizon
        
        total_var =  yyirf(t,i,1)^2+ yyirf(t,i,2)^2+ yyirf(t,i,3)^2;
        
        if (t ==1)||(t ==5)||(t ==20)||(t ==40)
            
         for  sh_ind = 1:3
             ratio(sh_ind)= round(yyirf(t,i,sh_ind)^2/total_var*1000)/10;
         end   
            
            disp([ num2str(t) '-th period,    ' ] );            
            disp(['    ' num2str(ratio(1)), ' [ % ],     ',...
                 num2str(ratio(2)), ' [ % ],       ',...
                 num2str(ratio(3)), ' [ % ]' ] );
        end
            
       end 
       disp('===========================================')        
    end    
    
    