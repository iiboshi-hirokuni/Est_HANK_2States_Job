function [ re_xt,  re_lik, weight_new, ESS  ] = resample_para(xt, lik_xt, npara,nsim, weight)

%% calculate ESS and rho

   ESS = nsim/(sum(weight.^2)/nsim);
   
   if ESS < nsim/2
       rho = 1;
       w1 = ones(1,nsim);
   else
       rho = 0;
       w1 = weight;
   end
            

prob_xt=exp(lik_xt(3,:)-max(lik_xt(3,:))); 

% disp(['prob_xt = ' num2str(prob_xt) ]);
save('prob','prob_xt' );

% disp(lik_xt(1,:));
% disp(prob_xt);

w = (prob_xt.*w1)/sum(prob_xt.*w1);

% disp(['weight = ' num2str(w) ]);

weight1 = nsim*w;

re_xt_1 = zeros(npara+3+2,nsim);



 for j=1:1:nsim    
          re_xt_1(:,j) = mnrnd(1, w) * [ xt;  lik_xt; weight1  ]'  ;    
 end  

re_xt = re_xt_1(1:npara,:);  %% resample pf para

re_lik = re_xt_1(npara+1:npara+3+1,:);

weight_new = re_xt_1(npara+4+1,:);

  disp([ 'ESS =  '  num2str(ESS) ]); 
  disp([ 'Num/2 =  '  num2str(nsim/2) ]);

% total = 0;
% 
% re_xt = xt;
% re_lik = lik_xt;
% 
% for i = 1:nparticles
%     
%     n = round(w(i)*nsim);
%     
%     if (n>0)&&(total<nsim)
%         
%         if (total + n) > nparticles
%             total_last = nparticles;
%         else
%             total_last = total + n;
%         end
%         
%         for j = total+1:total_last
%             re_xt(:,j)   = xt(:,i);
%             re_lik(:,j) = lik_xt(:,j);
%         end
%     end
%     
%     total = total + n;
% end