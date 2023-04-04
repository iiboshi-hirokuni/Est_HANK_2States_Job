% 
% 
%  Plot Figure of  Invariant Distributions (Steady State) of Hetero Agents
%
%

 w1 = sum(g_SS(:,1));
 w2 = sum(g_SS(:,2));
%  w3 = sum(g_SS(:,3));
%  w4 = sum(g_SS(:,4));
 
 disp( [ 'ratio of Agent 1=' num2str(w1/(w1+w2)) ]);
 disp( [ 'ratio of Agent 2=' num2str(w2/(w1+w2)) ]);
%  disp( [ 'ratio of Agent 3=' num2str(w3/(w1+w2+w3+w4)) ]);
%  disp( [ 'ratio of Agent 4=' num2str(w4/(w1+w2+w3+w4)) ]);

figure('Name','SS_dist','File','Dist_J2');
subplot(2,1,2);
 hold on
  l1= area(a, g_SS(:,1)+g_SS(:,2)  );
%   l1= area(a,((w1+w2+w3+w4)/w1).*g_SS(:,1));
    set(l1(1),'FaceColor', [0.65 0.65 1.00],'EdgeColor',[0 0.5 1])  % steady state

     y2=0.5/5;     y1 = 0:y2:0.5;
     plot( 6*ones(6,1), y1,'r-','LineWidth', 3 );
  hold off
   title('Stationary Distribution of All Agents ','FontSize',12);
   
   legend('All Agents','Mean','FontSize',12,'Box','off'); 
     xlabel({'Asset'},'FontSize',11);
     ylabel({'Density'},'FontSize',11);
    xlim([0 15]);
   ylim([0 0.25]);

% figure('Name','SS_dist');
subplot(2,1,1);


 hold on
  l1= area(a,((w1+w2)/w1).*g_SS(:,1));
    set(l1(1),'FaceColor', [0.65 0.65 1.00],'EdgeColor',[0 0.5 1])  % steady state
  l2= area(a,((w1+w2)/w2).*g_SS(:,2)) ;
    set(l2(1),'FaceColor',[0.25 1.0 1.00],'EdgeColor',[0.5 1 1 ],'FaceAlpha',0.6)  % steady state
%   l3= area(a,((w1+w2+w3+w4)/w3).*g_SS(:,3)) ;
%     set(l3(1),'FaceColor', [0.25 1.0 1.00],'EdgeColor',[0.5 1 1 ],'FaceAlpha',0.6)  % steady state  
%   l4= area(a,((w1+w2+w3+w4)/w4).*g_SS(:,4)) ;   
%     set(l4(1),'FaceColor', [0.65 1.0 0.65],'EdgeColor',[0.5 1 0.5],'FaceAlpha',0.6)  % steady state
%   l1d= plot(a,((w1+w2+w3+w4)/w1).*(g_SS(:,1)+simulated_gg(:,1)),'b--','LineWidth',2) ;
%   l2d= plot(a,((w1+w2+w3+w4)/w2).*(g_SS(:,2)+simulated_gg(:,2)),'r--','LineWidth',2) ;
%   l3d= plot(a,((w1+w2+w3+w4)/w3).*(g_SS(:,3)+simulated_gg(:,3)),'g--','LineWidth',2) ;
%   l4d= plot(a,((w1+w2+w3+w4)/w4).*(g_SS(:,4)+simulated_gg(:,4)),'g--','LineWidth',2) ;
  
  hold off
   title('Stationary Distribution of Four Job States ','FontSize',12);
   
    legend([l1 l2  ], {['Prod=0.2 ( ' num2str( round(1000*w1/(w1+w2))/10) ' % )' ],...
                       ['Prod=1.0 ( ' num2str( round(1000*w2/(w1+w2))/10) ' % )' ],...
                        }, 'FontSize',12,'Box','off');
%    if j==1
%      legend([l1 l1d l2 l2d l3 l3d l4 l4d], {'Low-Q-Worker (before shock)','Low-Q-Worker (after shock)',...
%                                   'Mid-Q-Worker (before shock)','Mid-Q-Worker (after shock)',...
%                                   'High-Q-Worker (before shock)','High-Q-Worker (after shock)'} ,'FontSize',9);
% %    end
     xlabel({'Asset'},'FontSize',11);
     ylabel({'Density'},'FontSize',11);
   xlim([0 15]);
   ylim([0 0.2]);