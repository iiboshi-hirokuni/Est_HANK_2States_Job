function  [ likenew pf] = step_solve_particle(P,O,tstart,yy,ZZ,HH,...
                            nparticles, zlbflag, policy_flag, npar, stock_shock,stock_state,mex)   


     S = steadystate(P);
     G = grids_even(O,P);
        
    [ pf, converged] = solve_model(P,S,G,O,zlbflag,policy_flag,tstart,mex);
         
   % structural shocks
     QQ = eye(3);
     QQ(1,1)= P.sigma_b^2; % preference
     QQ(2,2)= P.sigma_a^2; % TFP
     QQ(3,3)= P.sigma_r^2; % monetary policy   
     
     
%      HH = 1.0*diag([ P.sigma_y^2  P.sigma_pi^2  P.sigma_R^2 ]); 
   
   % particle filter
  
    if converged == 2
           likenew = -1E6;  
           disp('model is not converged.')
    
    else
%            disp('Start Particle Filter ')    
        
       [ log_lik] = ...
                fun_ParticleFilter(yy,ZZ,HH,QQ, O, G, S, P, pf,...
                      nparticles, zlbflag, policy_flag, npar, stock_shock,stock_state) ;  
  
        if isfinite(log_lik) % if log_lik not =  infinite or NaN
               likenew = log_lik;      
%                 disp( ['log likelihood is ' num2str(log_lik) ]);
        else    
               likenew = -1E6;   
               disp('log_lik is infinite')
        end     
    end    
%      
        
        
        
        