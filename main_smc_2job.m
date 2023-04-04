% 
%  Sequential Monte Carlo Methods 
%

%% load data


   load_data;
%    dataplot;

%% setting of  estimation
set_itr_disp = 1;
cc1        =   0.5 ;  % adjustment coefficient of SMC
N_Blocks = 1;     % Number of random Blocks of sampling 


%% setting of  environment of HANK model
  I   = 100;   % number of grids of one ASSET
  J   = 2;       % number of grids of states of JOB
  n_v = I*J + 1;    % number of JUMP variables (value function + inflation)
  n_g = I*J + 2 ;   % number of ENDOGENOUS state variables (distribution + monetary + Fiscal policy)
  n_p = 6;          % number of static relations: bond-market clearing, labor market clearing, consumption, output, total assets
  n_shocks = 3;     % number of SHOCKS, i.e., monetary policy shock, fiscal policy shock, TFP shock.
  nEErrors = n_v;   
  nVars = n_v + n_g + n_p;
  
  
  %%   
    Pr = fun_prior_setting(npara);
    parasim =  sample_pri(Pr,nsim,npara); 
    lik_stock = zeros(nsim,4);
    
    stock_accept_rate = zeros(nstage,nsim);
    psi = 0;
    weight = ones(1,nsim);
    priornew = zeros(nsim,1);
    parasim_stock   =  zeros(nsim, (npara+4),nstage);
    
    cc = cc1*diag(Pr.pstdd(1:npara) );
%     rrho1 = 0;
%     rrho2 = 0;
 
 %% start SMC^2   
for i = 1:nstage 
    disp( ' ' );   
  disp([ ' ' num2str(i), ' th-stage' ]) ;
  
 
 
%% Step 1. Correction


   %% update of psi of likelihood function 
        if i > 1
          psi = (i/nstage)^2-((i-1)/nstage)^2;
        else
          psi = (i/nstage)^2;
        end    
         disp( [ 'psi = ' num2str(psi) ]  )
   parfor k = 1:nsim  
%   for k = 1:nsim    


  if mod(i,set_itr_disp)==0
         disp( ' ' );   
         disp([ ' Step 1: Correction ' num2str(k), ' th-particle of ' num2str(i), ' th-stage ']) ;
         pause(0.05)
  end       
  
  
    % Projection method for solution of DSGE model  
     [check_para] = fun_check_para(parasim(k,:),Pr,npara);
     
     if check_para == 1
          likenew = -1E6;
           rrho1  =  0; 
           if mod(i,set_itr_disp)==0
                   disp([ 'out of bound of parameters' ]);
           end        
     elseif i == 1
         
              [ likenew,~,rrho1 ]  = cal_likelihood(YY, I,J,n_v,n_g,n_p,...
                                n_shocks, parasim(k,:),obs_ratio,def_switch);
 
     else
             dummy =  lik_stock(k,:);
             likenew = dummy(1,2);
             rrho1    = dummy(1,4);
     end        
                        
%       % calculate prior density  
        if i ==1
              priornew(k) =priodens(parasim(k,:), Pr.pmean, Pr.pstdd, Pr.pshape);   
              postnew = psi*likenew + priornew(k) ;
        else
             dummy =  lik_stock(k,:);
             post_old = dummy(1,1);
             postnew  = post_old + psi*likenew;
        end     
            
        % incremental weights
        incwt    =  psi*likenew;
     
       % save lik and posterior
        if mod(i,set_itr_disp)==0
               disp( [ 'post = ' num2str(postnew) ',   like = '  num2str(likenew) ...
                 ]);
             pause(0.05)
        end     
%                '  para ='  num2str(parasim(k,:))    ] );
       lik_stock(k,:)=[postnew likenew incwt rrho1];
     
  end    
% end 
 
   %% save file 
          file_name = ['./output/save_step1_para_'  num2str(nsim) '_'  num2str(i) ];                      
          save(file_name,'parasim','lik_stock');  
            


 %% Step 2. Selection
  disp( ' ' );   
  disp([ ' Step 2: Selection of ' num2str(i), ' th-stage ']);
      [ para_Resamp , lik_Resamp, weight, ESS ] = resample_para(parasim', lik_stock', npara,nsim, weight);
        para_Resamp = para_Resamp';
        lik_Resamp = lik_Resamp';    
        
        post_Resamp = zeros(nsim,4);
     
        
 
 %% Step 3. Mutation 
      para_new =zeros(nsim,npara);
      
       % variances of parameters of Algorithm 10 of Herbst and Schorfheide (2016, Ch5, p113)    
      if i >1
          V=cov(para_Resamp);
          V = V + diag( 1E-7*ones(npara,1));
          cc = cc1*chol(V,'upper');
      end
          
 for  j = 1:nsim
    
       % RW sampling of candidates of parameters    
       para_new(j,:) = para_Resamp(j,:) + randn(1,npara)*cc  ;
       para_new(j,:) = para_new(j,:).*(1-Pr.pmask)'+para_Resamp(j,:).*Pr.pmask';

       % Projection method for solution of DSGE model
      
 end  

  
%  for j = 1:nsim 
 parfor j = 1:nsim  

     % initial setting
     para_old = para_Resamp(j,:) ;
     post_old = lik_Resamp(j,1);
     post_Resamp(j,:) = [ lik_Resamp(j,1) lik_Resamp(j,2) 0 lik_Resamp(j,4)];
     %

      if mod(j,set_itr_disp)==0
             disp( ' ' );  
            disp([ ' Step 3: Mutation ' num2str(j), ' th-particle of ' num2str(i), '/' ...
                     num2str(nstage), ' th-stage ' ]) ;  
            pause(0.05)     
      end           
          
        [check_para] = fun_check_para(para_new(j,:),Pr,npara);
     
        if check_para == 1
               likenew = -1E6;
                 if mod(j,set_itr_disp)==0
                       disp([ 'out of bound of parameters' ]);
                      pause(0.05) 
                 end      
                   post_Resamp(j,:) = [ lik_Resamp(j,1) lik_Resamp(j,2) 0 lik_Resamp(j,4)];
        else
            
            para_select =rand(1,npara);       

            for k=1:N_Blocks
                para_block_select=(para_select<k/N_Blocks)-(para_select<(k-1)/N_Blocks);
                para_candiate = para_block_select.*para_new(j,:)+(1-para_block_select).*para_old;
           
                [ likenew,~, rrho2] =  cal_likelihood(YY, I,J,n_v,n_g,n_p,...
                                   n_shocks, para_candiate, obs_ratio,def_switch );                  
               
                 %  calculate prior density   
                  priornew(j) =priodens(para_candiate , Pr.pmean, Pr.pstdd, Pr.pshape);   
                  postnew =  (i/nstage)^2*likenew + priornew(j); 
        
                 % MH_step
                 r = min(1,exp(postnew-post_old ));   
                  if (rand < r)   
                       weight(j)=1; 
                       stock_accept_rate(i,j) = 1;
                       para_Resamp(j,:) = para_new(j,:) ;   para_old =  para_new(j,:);          
                       post_Resamp(j,:) = [postnew  likenew   1 rrho2  ];
                         if mod(j,set_itr_disp)==0
                            disp( [ 'post = ' num2str(postnew) ',   like = '  num2str(likenew) ...
                                ]);
                        end         
                  else
                        post_Resamp(j,:) = [ lik_Resamp(j,1) lik_Resamp(j,2) 0 lik_Resamp(j,4)];
                        if mod(j,set_itr_disp)==0
                             disp( ['no change: ', 'post = ', num2str(lik_Resamp(j,1)),...
                               ',   like = '  num2str(lik_Resamp(j,2)) ...
                              ]);
                       end       
                  end 
            end    
        end       
end   
         R =   mean(stock_accept_rate(i,:))*100;
         % scaling factor of Algorithm 10 of Herbst and Schorfheide (2016, Ch5, p113)        
         x=6; 
         cc1 = cc1*(0.95+0.1*exp(x*(R-0.25))/(1+exp(x*(R-0.25))) );

         disp( ' ' ); 
         disp( [  num2str(i), ' th-stage_s  Accept Rate = ' num2str(R) ' %' ]);   
         disp( [  num2str(i), ' th-stage_s  scaling factor = ' num2str(cc1)  ]); 
         pause(0.05)
         
       
         parasim = para_Resamp;
         lik_stock = post_Resamp;
         
         %=====================================
          parasim_stock(:,1:npara,i)                     =  parasim;
          parasim_stock(:,npara+1:npara+4 ,i)   =  lik_stock;          
         %=====================================

        %% save file 
          file_name = ['./output/save_temp_para_' num2str(nsim) '_'  num2str(i) ];                      
          save(file_name,'para_Resamp','post_Resamp', 'stock_accept_rate',"parasim_stock");
          
        %% plot graph
          plot_dist([parasim,lik_stock], npara,Pr,i)
          plot_trace([parasim,lik_stock], npara,Pr,i)
          
        
          
%%  end of Iterations           
end 

%%
if data_country ==1
    country = {'JP'};
else
    country = {'US'};
end    

  file_name = ['./output/save_para_HANK_' char(country) '_' num2str(nsim) '_' num2str(nstage) ];
                      
     save(file_name,'para_Resamp','post_Resamp', 'stock_accept_rate'...
           );  %,'stock_shock','stock_state');
  %%     
   plot_dist([parasim,lik_stock], npara,Pr,i)
   plot_trace([parasim,lik_stock], npara,Pr,i)
   
  stats_sample_para_2job  
  
%% time of computing
 dec = 10^1;
T = toc(tstart);
hh = floor(T/3600);
mm = floor(T/60)-60*hh;
ss = round((T-60*mm-3600*hh)*dec)/dec;

hh = num2str(hh);
mm = num2str(mm);
ss = num2str(ss);

display(['Total Computing time: ' hh 'h' mm 'm' ss 's']);      
         
                  
 