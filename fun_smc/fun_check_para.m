function [check] = fun_check_para(para,Pr,npara)

% npara=20;

Pr.upper;
Pr.lower;

check = 0;

for i =1:npara
   if (para(i) > Pr.upper(i))&&(Pr.pmask(i)==0)
%        para(i)
       check = 1;
       break;
   elseif ( para(i) < Pr.lower(i))&&(Pr.pmask(i)==0)
%        para(i)
       check = 1;
       break;
   end 
end


