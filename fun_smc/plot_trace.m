% 
% 
%  Plot traces of sampling parameters
% 
% 
% 

function plot_trace(parasim, npara, P,ii)


para_names_p = {'coefrra','frisch','adjfricshgridfrac','ceselast',...
                    'priceadjust','taylor_inflation','taylor_outputgap',...
                    'labtax','govbondtarget','lumptransferpc','govbcrule_fixnomB',...
                    'ssigma_MP','ttheta_MP','ssigma_FP','ttheta_FP','ssigma_TFP',...
                    'ttheta_TFP','TP_22', ...         
                    'post', 'lik',...
                    'accept_rate(1)', 'accept_rate(2)','accept_rate(3)',...
                    'accept_rate(4)', 'accept_rate(5)','accept_rate(6)',...
                    'accept_rate(7)', 'accept_rate(8)','accept_rate(9)',...
                    'accept_rate(10)', 'accept_rate(11)','accept_rate(12)'...
                     };
 
 
%% draw graph both of Prior and posterior densities
    k=0;
    j = 1;
     figure([100])
    % figure('Name',[num2str(ii*1000)],'Position',[125,150,600,600])  
 for i =1:1:npara
     
  if P.pmask(i)==0 
    
  if mod(i,5*4)==1
     k = k +1 ;

     
%       figure('Position',[250,250,1000,250], 'Name', 'Fitting of Macro Factors','NumberTitle','off', 'FileName','Fig3001')  
  end
 
 subplot(4,4,j) % 4(çs)Å~4(óÒ)
           plot(parasim(:,i)','LineStyle','-','Color','b',...
          'LineWidth',1.0);      
     title( char(para_names_p(i)),'FontSize',12) ;      
      j= j+1; 
  end 
 end
 
 subplot(4,4,j) % 4(çs)Å~4(óÒ)
           plot(parasim(:,npara+4)','LineStyle','-','Color','b',...
          'LineWidth',1.0);      
     title( 'rrho','FontSize',12) ;      
 
 
 subplot(4,4,j+1) % 4(çs)Å~4(óÒ)
           plot(parasim(:,npara+1)','LineStyle','-','Color','b',...
          'LineWidth',1.0);      
     title( char(para_names_p(npara+1)),'FontSize',12) ;      
     
 subplot(4,4,j+2) % 4(çs)Å~4(óÒ)
           plot(parasim(:,npara+2)','LineStyle','-','Color','b',...
          'LineWidth',1.0);      
     title( char(para_names_p(npara+2)),'FontSize',12) ;      
 
 
 pause(0.5)





