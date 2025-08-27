load('h_vec.mat');
load('h_eps_vec.mat');

subplot(1,2,1);

plot(h_vec(100:end,1),h_vec(100:end,2),'LineWidth',1.2);
ylim([1,3]);
xlim([0,90]);
title('$h_{2}(m)$','Interpreter','latex')
xlabel('$m$','Interpreter','latex')
set(gca,'FontSize',12);

subplot(1,2,2);

plot(h_eps_vec(100:end,1),h_eps_vec(100:end,5),'LineWidth',1.2);
ylim([1,1.1]);
xlim([0,90]);
title('$\overline{h}_{2,0.05}(m)$','Interpreter','latex')
xlabel('$m$','Interpreter','latex')
set(gca,'FontSize',12);
x0=10;
y0=10;
width=1000;
height=380;

set(gcf,'position',[x0,y0,width,height])
set(gcf,'color','w');
%export_fig('h2_examples','-eps');