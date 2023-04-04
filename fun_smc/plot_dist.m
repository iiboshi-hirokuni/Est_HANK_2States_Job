% 
% 
%  Plot distributions of parameters
% 
% 
% 
% 

function plot_dist(parasim, npara, P,ii)

% parasim = [ para_Resamp  post_Resamp stock_accept_rate' ] ;

both_plot = 1; % on =1, off =0 

para_names_p = {'coefrra','frisch','adjfricshgridfrac','ceselast',...
                    'priceadjust','taylor_inflation','taylor_outputgap',...
                    'labtax','govbondtarget','lumptransferpc','govbcrule_fixnomB',...
                    'ssigma_MP','ttheta_MP','ssigma_FP','ttheta_FP','ssigma_TFP',...
                    'ttheta_TFP', 'TP_22',...         
                    'post', 'lik',...
                    'accept_rate(1)', 'accept_rate(2)','accept_rate(3)',...
                    'accept_rate(4)', 'accept_rate(5)','accept_rate(6)',...
                    'accept_rate(7)', 'accept_rate(8)','accept_rate(9)',...
                    'accept_rate(10)', 'accept_rate(11)','accept_rate(12)'...
                     };

%% draw graph both of Prior and posterior densities
  k=0;
   
 j=1 ;
   
%   figure('Name',[num2str(ii*100)],'Position',[25,50,600,600] )    
    figure([1000])     
  
 for i =1:1:npara
     
  if P.pmask(i)==0   
    
  if mod(i,5*4)==1
     k = k +1 ;
 
%       figure('Position',[250,250,1000,250], 'Name', 'Fitting of Macro Factors','NumberTitle','off', 'FileName','Fig3001')  
  end
 
 subplot(4,4,j) % 4(çs)Å~4(óÒ)
 
%      [density0,x0]  = ksdensity(parasim_pri(1:end,i));
     [density,x]  = ksdensity(parasim(1:end,i));
     
%            plot(x,density,'LineStyle','-','Color','b',...
%           'LineWidth',2.5);
%            y=linspace(min(density),max(density),size(x,2));
%    
           hist(parasim(1:end,i)) ; 
           hold on
                data = parasim(1:end,i);
                t_x = [min(data) mean(data) mean(data) mean(data) max(data) ];
                t_y = [0  0 100  0  0 ];
                 plot(t_x, t_y, 'LineWidth',2,'Color','red');
           hold off  
    
      
     title( char(para_names_p(i)),'FontSize',12) ; 
%     if i == 1 && mode_flag==1
%      legend('density','max', 'second', 'third');
%     end  
     j= j+1; 
  end
   
 end
 
 subplot(4,4,j) % 4(çs)Å~4(óÒ)
 
%      [density0,x0]  = ksdensity(parasim_pri(1:end,i));
     [density,x]  = ksdensity(parasim(1:end,npara+4));
     
%            plot(x,density,'LineStyle','-','Color','b',...
%           'LineWidth',2.5);
%            y=linspace(min(density),max(density),size(x,2)); 
%            
             hist(parasim(1:end,npara+4)) ;
       hold on
                data = parasim(1:end,npara+4);
                t_x = [min(data) mean(data) mean(data) mean(data) max(data) ];
                t_y = [0  0 100  0  0 ];
                 plot(t_x, t_y, 'LineWidth',2,'Color','red');
           hold off 
%      title( char(para_names_p(npara+2)),'FontSize',12) ; 

           
     title( 'Dis Fac','FontSize',12) ; 
 
 
 
 subplot(4,4,j+1) % 4(çs)Å~4(óÒ)
 
%      [density0,x0]  = ksdensity(parasim_pri(1:end,i));
     [density,x]  = ksdensity(parasim(1:end,npara+1));
     
%            plot(x,density,'LineStyle','-','Color','b',...
%           'LineWidth',2.5);
%            y=linspace(min(density),max(density),size(x,2));        
%      title( char(para_names_p(npara+1)),'FontSize',12) ; 
%      
       hist(parasim(1:end,npara+1)) ;
       hold on
                data = parasim(1:end,npara+1);
                t_x = [min(data) mean(data) mean(data) mean(data) max(data) ];
                t_y = [0  0 100  0  0 ];
                 plot(t_x, t_y, 'LineWidth',2,'Color','red');
           hold off 
     title( char(para_names_p(npara+1)),'FontSize',12) ; 

 
 subplot(4,4,j+2) % 4(çs)Å~4(óÒ)
 
%      [density0,x0]  = ksdensity(parasim_pri(1:end,i));
     [density,x]  = ksdensity(parasim(1:end,npara+2));
     
           plot(x,density,'LineStyle','-','Color','b',...
          'LineWidth',2.5);
           y=linspace(min(density),max(density),size(x,2));   
           
             hist(parasim(1:end,npara+2)) ;
       hold on
                data = parasim(1:end,npara+2);
                t_x = [min(data) mean(data) mean(data) mean(data) max(data) ];
                t_y = [0  0 100  0  0 ];
                 plot(t_x, t_y, 'LineWidth',2,'Color','red');
           hold off 
     title( char(para_names_p(npara+2)),'FontSize',12) ; 

%      title( char(para_names_p(npara+2)),'FontSize',12) ; 
 
 pause(0.5)




