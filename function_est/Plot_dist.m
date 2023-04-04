
figure('Name','SS_dist');

shock_name=char('Monetary Policy shock', 'Fiscal Policy shock', 'TFP shock');

for j=1:3
subplot(2,2,j);
  impact1=mm*impact(:,j);

  [simulated_gg0,~] = simulate(G1,impact1,T,N,vAggregateShock,'implicit',trans_mat,[n_v+1:n_v+n_g-2]);
  simulated_gg = reshape(simulated_gg0(:,2),I,J);
  simulated_gg(I,:)=zeros(2,1);   

  w1 = sum(g_SS(:,1));
  w2 = sum(g_SS(:,2));
  
  if j==1
    w1/(w1+w2)
    w2/(w1+w2)
  end
  
 hold on
  l1= area(a,((w1+w2)/w1).*g_SS(:,1));
    set(l1(1),'FaceColor', [0.65 0.65 1.00],'EdgeColor',[0 0.5 0.1])  % steady state
  l2= area(a,((w1+w2)/w2).*g_SS(:,2)) ;
    set(l2(1),'FaceColor', [1.0 0.65 0.65],'EdgeColor',[1 0.5 0.5],'FaceAlpha',0.6)  % steady state
%   l1= plot(a,((w1+w2)/w1).*g_SS(:,1),'b-','LineWidth',2);
%   l2= plot(a,((w1+w2)/w2).*g_SS(:,2),'r-','LineWidth',2) ;
  
  l3= plot(a,((w1+w2)/w1).*(g_SS(:,1)+simulated_gg(:,1)),'b--','LineWidth',2) ;
  l4= plot(a,((w1+w2)/w2).*(g_SS(:,2)+simulated_gg(:,2)),'r--','LineWidth',2) ;
  hold off
   title( [shock_name(j,:)],'FontSize',12);
   if j==1
    legend([l1 l3 l2 l4], {['Prod = 0.2 (before shock) (' num2str(round(1000*w1/(w1+w2))/10) ' % )'  ],...
                           'Prod = 0.2 (after shock)',...
                           ['Prod = 1.0 (before shock) (' num2str(round(1000*w2/(w1+w2))/10) ' % )'  ],...
                           'Prod = 1.0 (after shock)'} ,'FontSize',12);
   end
     xlabel({'Asset'},'FontSize',11);
     ylabel({'Density'},'FontSize',11);
  xlim([0 15]);
  if j==3
      ylim([-0.01 0.2])
  end    

end


