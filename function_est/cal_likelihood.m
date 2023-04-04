function  [rloglh,obsmean, rrho] = cal_likelihood(YY, I,J,n_v,n_g,n_p,n_shocks,para,obs_ratio,def_switch)

   rrho=0;
   obsmean = 0;

   [G1, impact,inv_state_red,from_spline, eu, rrho] = solve_HANK(I,J,n_v,n_g,n_p,n_shocks,para);  
   
   
   if (eu(1) == 1)&&(eu(2)==1) 
     T = 40;
     N = 40;
     dt = T/N;
      vAggregateShock	= zeros(1,N);
      trans_mat = inv_state_red*from_spline;    
        
     ZZ0 = zeros(3,size(trans_mat,1));
     ZZ0(1,n_v+n_g+4) = obs_ratio(1); %output
     ZZ0(2,n_v)       = obs_ratio(2); %inflation
     ZZ0(3,n_v+n_g+6) = obs_ratio(3); %interest rate
     ZZ1 = zeros(3,size(trans_mat,1));
     if def_switch ==1  
       ZZ1(1,n_v+n_g+4) = obs_ratio(1); %output
     end
     ZZ = [ZZ0*trans_mat -ZZ1*trans_mat] ; 
   
     m = 2;
     [rloglh,obsmean] = eval_HANK(G1,impact,T,N,'implicit',YY, ZZ,m);
   else
     rloglh = -1E6; 
   end    
     
   
   