% 
% 
%  Plot Figure of  Invariant Distributions (Steady State) of Hetero Agents
%
%

figure('Name','policy_functions');
 subplot(2,2,1);
hold on
 l1= plot(a,(V_SS(:,1)),'b--','LineWidth',2);
 l2= plot(a,(V_SS(:,2)),'r-','LineWidth',2) ; 
 hold off
 title('Value Function','FontSize',14);
  legend([l2 l1], {'High-Pro.-Worker','Low-Pro.-Worker'} ,'FontSize',12);
 xlabel('Asset','FontSize',12);
 xlim([0 15]);
%   ylim([-24600 -24580])
 
subplot(2,2,2);
hold on
 l1= plot(a,(c_SS(:,1)),'b--','LineWidth',2);
 l2= plot(a,(c_SS(:,2)),'r-','LineWidth',2) ; 
 hold off
 title('Policy Function of Consumption','FontSize',14);
%  legend([l1 l2], {'Low-Skill-Labor','High-Skill-Labor'} ,'FontSize',12);
 xlabel('Asset','FontSize',12);
 xlim([0 15]);
 
subplot(2,2,3);
hold on
 l1= plot(a,(s_SS(:,1)),'b--','LineWidth',2);
 l2= plot(a,(s_SS(:,2)),'r-','LineWidth',2) ; 
 hold off
 title('Policy Function of Saving','FontSize',14);
%  legend([l1 l2], {'Low-Skill','High-Skill'} ,'FontSize',12);
 xlabel('Asset','FontSize',12);
 xlim([0 15]);
 
subplot(2,2,4);
hold on
 l1= plot(a,(h_SS(:,1)),'b--','LineWidth',2);
 l2= plot(a,(h_SS(:,2)),'r-','LineWidth',2) ; 
 hold off
 title('Policy Function of Working Hour','FontSize',14);
%  legend([l1 l2], {'Low-Skill','High-Skill'} ,'FontSize',12);
 xlabel('Asset','FontSize',12);
 xlim([0 15]);
 ylim([0.2 0.4])
 