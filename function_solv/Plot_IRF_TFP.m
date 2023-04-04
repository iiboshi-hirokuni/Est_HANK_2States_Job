%% Step 5: Simulate Impulse Response Functions
fprintf('Simulating Model...\n');
t0 = tic;

impact1=impact(:,3);

T = 40;
N = 100;
dt = T/N;
vAggregateShock	= zeros(1,N);
vAggregateShock(:,1) = 1/sqrt(dt);
trans_mat = inv_state_red*from_spline;
[simulated,vTime] = simulate(G1,impact1,T,N,vAggregateShock,'implicit',trans_mat,[n_v,n_v+n_g-2:n_v+n_g+6]);

[simulated_gg0,~] = simulate(G1,impact1,T,N,vAggregateShock,'implicit',trans_mat,[n_v+1:n_v+n_g-2]);
simulated_gg = reshape(simulated_gg0(:,J),I,J);
simulated_gg(I,:)=zeros(J,1);

fprintf('...Done!\n');
fprintf('Time to simulate model: %2.4f seconds\n\n\n',toc(t0));

J2.inflation = simulated(1,:)';
% J2.MP_shock = simulated(2,:)';
J2.TFP_shock = simulated(4,:)';
J2.consumption = (simulated(2+5,:)')/vars_SS(n_v+n_g+3);
J2.Y = (simulated(2+6,:)')/vars_SS(n_v+n_g+4);
J2.lab_sup = (simulated(2+4,:)')/vars_SS(n_v+n_g+2);
J2.wage = simulated(2+3,:)'/vars_SS(n_v+n_g+1);
J2.interest = simulated(2+8,:)'-  simulated(1,:)';
J2.bond     = simulated(2+7,:)';

%% (optional) Step 7: Plot relevant IRFs
line_style = '--';
color = 'blue';
line_style = '--';
color = 'red';
% line_style = ':';
% color = 'black';

figure('Name','TFP Shock');
subplot(3,3,1);
hold on;
plot(vTime,J2.TFP_shock,'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('TFP Shock','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 T]);
grid on;
hold off;

subplot(3,3,2);
hold on;
plot(vTime,J2.interest,'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Interest Rate','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 T]);
grid on;
hold off;

subplot(3,3,3);
hold on;
plot(vTime,J2.inflation,'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Inflation','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 T]);
grid on;
hold off;

subplot(3,3,4);
hold on;
plot(vTime,J2.consumption,'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Consumption','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 T]);
grid on;
hold off;

subplot(3,3,5);
hold on;
plot(vTime,J2.Y,'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('GDP','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 T]);
grid on;
hold off;

subplot(3,3,6);
hold on;
plot(vTime,J2.lab_sup,'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Labor Supply','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 T]);
grid on;
hold off;

subplot(3,3,7);
hold on;
plot(vTime,J2.wage,'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Wage','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 T]);
grid on;
hold off;

subplot(3,3,8);
hold on;
plot(vTime,J2.bond,'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Assets','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 T]);
grid on;
hold off;


%%

if data_country == 1
  save('IRF_TFP_JP_J2.mat','J2');
else
   save('IRF_TFP_US_J2.mat','J2'); 
end    

% plot_G_ss

