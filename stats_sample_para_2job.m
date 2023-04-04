
para_names_p = char('coefrra','frisch','adjfricshgridfrac','ceselast',...
                    'priceadjust','taylor_inflation','taylor_outputgap',...
                    'labtax','govbondtarget','lumptransferpc','govbcrule_fixnomB',...
                    'ssigma_MP','ttheta_MP','ssigma_FP','ttheta_FP','ssigma_TFP',...
                    'ttheta_TFP','TP_labor_22', ...                                         
                    'post', 'lik','','rrho',...
                    'accept_rate(1)', 'accept_rate(2)','accept_rate(3)',...
                    'accept_rate(4)', 'accept_rate(5)','accept_rate(6)',...
                    'accept_rate(7)', 'accept_rate(8)','accept_rate(9)',...
                    'accept_rate(10)', ...
                    'accept_rate(11)', 'accept_rate(12)','accept_rate(13)',...
                    'accept_rate(14)', 'accept_rate(15)','accept_rate(16)',...
                    'accept_rate(17)', 'accept_rate(18)','accept_rate(19)',...
                    'accept_rate(20)', ...
                   'accept_rate(21)', 'accept_rate(22)','accept_rate(23)',...
                    'accept_rate(24)', 'accept_rate(25)','accept_rate(26)',...
                    'accept_rate(27)', 'accept_rate(28)','accept_rate(29)',...
                    'accept_rate(20)', ...
                    'accept_rate(31)', 'accept_rate(32)','accept_rate(33)',...
                    'accept_rate(34)', 'accept_rate(35)','accept_rate(36)',...
                    'accept_rate(37)', 'accept_rate(38)','accept_rate(39)',...
                    'accept_rate(40)', ...
                   'accept_rate(41)', 'accept_rate(42)','accept_rate(43)',...
                    'accept_rate(44)', 'accept_rate(45)','accept_rate(46)',...
                    'accept_rate(47)', 'accept_rate(48)','accept_rate(49)',...
                    'accept_rate(50)', ...
                    'accept_rate(51)', 'accept_rate(52)','accept_rate(53)',...
                    'accept_rate(54)', 'accept_rate(55)','accept_rate(56)',...
                    'accept_rate(57)', 'accept_rate(58)','accept_rate(59)',...
                    'accept_rate(10)' ...
                    );
                 

 parasim = [ para_Resamp post_Resamp(:,1:4) stock_accept_rate' ] ;



%% Output posterior file

   % Calculating of Posterior estimates 

a  =0.90;   % (1-2*rate) percentage of Credibile interval of Parameters  
rate = (1-a)/2;


% calculation of posterior estimates of parameters
%

sort_para=zeros( nsim, npara+4+nstage );


  for i=1:1:(npara+2+nstage)
     sort_para(:,i) = sort(parasim(:,i),1);
  end


%  The Fisrt equation
  para_low = sort_para(ceil(( nsim )*rate),:); 
  para_up  = sort_para(floor(( nsim )*(1-rate)),:);

iBm = min([500, nsim/2]); 


fprintf( [ '\n\n  nstage = ' num2str(nstage) ', particle (para) = ' num2str(nsim) ...
                   ] ); 
% fprintf( [ '\n\n  Measurement error = ' num2str(sqrt(HH(1,1))) ', ' num2str(sqrt(HH(2,2))) ...
%                   ', ' num2str(sqrt(HH( 3,3))) ] );              
              
fprintf('\n--------------------------------------------');
fprintf('------------------------------------');
fprintf('\n\n                        [ESTIMATION RESULT]');
fprintf('\n--------------------------------------------');
fprintf('------------------------------------');
fprintf('\n  Parameter           Mean        Stdev     ');
fprintf('95%%Low     95%%Up    Geweke     Inef.');
fprintf('\n----------------------------------');
fprintf('--------------------------------------------\n');

if nsim > 25
for i = 1:(npara)
 if (Pr.pmask(i)==0) 
     fprintf('%s %10.4f  %10.4f %10.4f %10.4f %9.3f %9.3f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      fGeweke(parasim(:,i), iBm), ...
      ftsvar(parasim(:,i), iBm)/var(parasim(:,i)) ]  );
 end
  
end
      fprintf('----------------------------------');
      fprintf('-------------------------------------------- \n');
for i = (npara+1):(npara+4+nstage)
% fprintf( '%s ',  para_names_p(i)  );
fprintf('%s %10.4f  %10.4f %10.4f %10.4f %9.3f %9.3f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      fGeweke(parasim(:,i), iBm), ...
      ftsvar(parasim(:,i), iBm)/var(parasim(:,i)) ]  );
  if (i == npara)||(i== npara+2)
      fprintf('----------------------------------');
      fprintf('-------------------------------------------- \n');
  end    
end

else
  for i = 1:(npara)
      if (Pr.pmask(i)==0)   
% fprintf( '%s ',  para_names_p(i)  );
fprintf('%s %10.4f  %10.4f %10.4f %10.4f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      ]  );
      end
  end
  fprintf('\n----------------------------------');
fprintf('--------------------------------------------\n');

  for i = (npara+1):(npara+4+nstage)
     % fprintf( '%s ',  para_names_p(i)  );
      fprintf('%s %10.5f  %10.4f %10.5f %10.5f \n' ,...
         para_names_p(i,:), ...   
        [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      ]  );
  end
  
end    


fprintf('-----------------------------------');
fprintf('---------------------------------------------\n');

 
est_date = datestr(date);         
result_name = ['./output/ESTIMATE_HANK_J2_', num2str(nstage), ...
    '-',num2str(nsim), '-',num2str(data_country), '-', est_date,'.txt'];          

fileID = fopen(result_name,'w');
fprintf(fileID, [ '\n\n  nstage = ' num2str(nstage) ', particle (para) = ' num2str(nsim) ...
                    ] ); 

% fprintf(fileID, [ '\n\n  Measurement error = ' num2str(sqrt(HH(1,1))) ', ' num2str(sqrt(HH(2,2))) ...
%                   ', ' num2str(sqrt(HH( 3,3))) ] );                  

fprintf(fileID,'\n\n                        [ESTIMATION RESULT]');
fprintf(fileID,'\n----------------------------------');
fprintf(fileID,'------------------------------------');
fprintf(fileID,'\nParameter         Mean        Stdev     ');
fprintf(fileID,'95%%Low     95%%Up    Geweke     Inef.');
fprintf(fileID,'\n----------------------------------');
fprintf(fileID,'--------------------------------------------\n');

if nsim > 25
for i = 1:(npara)
 if (Pr.pmask(i)==0)   
     fprintf(fileID,'%s %10.4f  %10.4f %9.3f %9.3f %9.3f %9.3f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      fGeweke(parasim(:,i), iBm), ...
      ftsvar(parasim(:,i), iBm)/var(parasim(:,i)) ]  );
 end
end
for i = (npara+1):(npara+3+nstage)
% fprintf( '%s ',  para_names_p(i)  );
fprintf(fileID,'%s %10.4f  %10.4f %9.3f %9.3f %9.3f %9.3f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      fGeweke(parasim(:,i), iBm), ...
      ftsvar(parasim(:,i), iBm)/var(parasim(:,i)) ]  );
end

else
  for i = 1:(npara)
      if (Pr.pmask(i)==0) 
% fprintf( '%s ',  para_names_p(i)  );
fprintf(fileID,'%s %10.4f  %10.4f %9.3f %9.3f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      ]  );
      end
  end
  for i = (npara+1):(npara+2+nstage)
     % fprintf( '%s ',  para_names_p(i)  );
fprintf(fileID,'%s %10.4f  %10.4f %9.3f %9.3f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      ]  );
  end
  
end    

fprintf(fileID,'-----------------------------------');
fprintf(fileID,'-----------------------------------');
fclose(fileID);
