clc;
close all;
clear;
%% α=0.2：0.2：1.0 五条曲线图，论文中一个挤压态误码率的实现
% % define squeezed operator
% r = 1;
% theta = 0:pi/8:2*pi;
% z = r*exp(1j*theta);
% % define displaced operator
% deta_1 = 0.2; % ∆ = 0.2
% eipxi = pi/6; % ξ= pi/6
% alpha_1 = deta_1*exp(1j*eipxi);
% % define 3-PPM modulation
%   K = 3;
% % calculate Pe
% inner_pro_1 = (1/cosh(r))*exp(-(abs(alpha_1)).^2*f(r,theta));
% Pe_1 = 1-((sqrt(1+(K-1).*inner_pro_1)+(K-1).*sqrt(1-inner_pro_1)).^2)/K^2;
% figure(1)
% semilogy(theta,Pe_1,'k--+')
% hold on
% 
% %  ∆ = 0.4
% deta_2 = 0.4; % ∆ = 0.5
% eipxi = pi/6; % ξ= pi/6
% alpha_2 = deta_2*exp(1j*eipxi);
% inner_pro_2 = (1/cosh(r))*exp(-(abs(alpha_2)).^2*f(r,theta));
% Pe_2 = 1-((sqrt(1+(K-1).*inner_pro_2)+(K-1).*sqrt(1-inner_pro_2)).^2)/K^2;
% plot(theta,Pe_2,'k:*')
% hold on
% 
% 
% %  ∆ = 0.6
% deta_3 = 0.6; % ∆ = 0.5
% eipxi = pi/6; % ξ= pi/6
% alpha_3 = deta_3*exp(1j*eipxi);
% inner_pro_3 = (1/cosh(r))*exp(-(abs(alpha_3)).^2*f(r,theta));
% Pe_3 = 1-((sqrt(1+(K-1).*inner_pro_3)+(K-1).*sqrt(1-inner_pro_3)).^2)/K^2;
% plot(theta,Pe_3,'k:p')
% hold on
% 
% %  ∆ = 0.8
% deta_4 = 0.8; % ∆ = 0.5
% eipxi = pi/6; % ξ= pi/6
% alpha_4 = deta_4*exp(1j*eipxi);
% inner_pro_4 = (1/cosh(r))*exp(-(abs(alpha_4)).^2*f(r,theta));
% Pe_4 = 1-((sqrt(1+(K-1).*inner_pro_4)+(K-1).*sqrt(1-inner_pro_4)).^2)/K^2;
% %plot(theta,Pe_4,'k:s')
% hold on
% 
% %  ∆ = 1.0
% deta_5 = 1.0; % ∆ = 0.5
% eipxi = pi/6; % ξ= pi/6
% alpha_5 = deta_5*exp(1j*eipxi);
% inner_pro_5 = (1/cosh(r))*exp(-(abs(alpha_5)).^2*f(r,theta));
% Pe_5 = 1-((sqrt(1+(K-1).*inner_pro_5)+(K-1).*sqrt(1-inner_pro_5)).^2)/K^2;
% %plot(theta,Pe_5,'k-^')
% grid on
% 
% xlim([0 2*pi])
% ylim([0 0.2])
% legend('show')
% legend('α = 0.2','α = 0.4','α = 0.6','α = 0.8','α = 1.0','Location','North')
% % 
% set(gca,'XTick',[0:pi/2:2*pi])
% set(gca,'xtickLabel',{'0','π/2','π','3π/2','2π'})
% xlabel('Squeezed angle θ [rad]')
% ylabel('Error probability Pe')
% 
% %% 不同 r 情况下误码率比较，r=0.1,0.5,1.0,  论文一个挤压态不同挤压程度（r）图的实现
% Ns = 0.02:0.5:5;
% Nr = Ns/log2(K);
% theta = pi;
% r = 0;
% inner_pro = (1/cosh(r))*exp(-(Ns-(sinh(r))^2)*f(r,theta));
% Pe = 1-((sqrt(1+(K-1).*inner_pro)+(K-1).*sqrt(1-inner_pro)).^2)/K^2;
% figure(2)
% semilogy(Ns,Pe,'k:s')
% xlabel('Received photons Ns [photons/symbol]')
% ylabel('Error probability Pe')
% 
% 
% hold on
% 
% Ns = 0.02:0.5:5;
% Nr = Ns/log2(K);
% theta = pi;
% r = 0.1;
% inner_pro = (1/cosh(r))*exp(-(Ns-(sinh(r))^2)*f(r,theta));
% Pe = 1-((sqrt(1+(K-1).*inner_pro)+(K-1).*sqrt(1-inner_pro)).^2)/K^2;
% figure(2)
% semilogy(Ns,Pe,'k--p')
% xlabel('Received photons Ns [photons/symbol]')
% ylabel('Error probability Pe')
% 
% legend('r=0.1,θ=pi')
% 
% hold on
% 
% 
% Ns = 0.28:0.5:8;
% Nr = Ns/log2(K);
% r = 0.3;
% inner_pro = (1/cosh(r))*exp(-(Ns-(sinh(r))^2)*f(r,theta));
% Pe = 1-((sqrt(1+(K-1).*inner_pro)+(K-1).*sqrt(1-inner_pro)).^2)/K^2;
% semilogy(Nr,Pe,'k:*')
% hold on
% 
% 
% Ns = 0.3:0.5:99;
% Nr = Ns/log2(K);
% r = 0.5;
% inner_pro = (1/cosh(r))*exp(-(Ns-(sinh(r))^2)*f(r,theta));
% Pe = 1-((sqrt(1+(K-1).*inner_pro)+(K-1).*sqrt(1-inner_pro)).^2)/K^2;
% semilogy(Nr,Pe,'k--+')
% xlim([0 4])
% ylim([0 0.15])
% legend('coherent(r=0)','r=0.1,   θ=π','r=0.3，θ=π','r=0.5，θ=π')

%% 一个挤压态Winger函数投影，椭圆
% define squeezed operator
r = 0.5;
theta = pi;
z = r*exp(1j*theta);
% define displaced operator
deta = 4; % ∆ = 0.2
eipxi = 0; % ξ= pi/6
alpha = deta*exp(1j*eipxi);
% define 3-PPM modulation
K = 3;
% p & q
q = deta*cos(eipxi);
p = deta*sin(eipxi);

% Variance matrix
V = zeros(2,2);
V(1,1) = (cosh(r))^2+(sinh(r))^2+cos(theta)*sinh(2*r);
V(1,2) = sin(theta)*sinh(2*r);
V(2,1) = sin(theta)*sinh(2*r);
V(2,2) = (cosh(r))^2+(sinh(r))^2-cos(theta)*sinh(2*r);

% Wigner Function
x = -6:0.01:6;
y = -6:0.01:6;
[ X, Y ] = meshgrid( x, y );
W = (1/2*pi)*exp(-0.5*(V(2,2).*(X-q).^2+V(1,1).*(Y-p).^2-2*V(1,2)*(X-q).*(Y-p)));
%figure(1)
subplot(2,1,1)
mesh(x,y,W) %mesh(x,y,W,'EdgeColor','k')

hold on
% mesh(x,y,0.1*ones(121,121))

% view([0,0,1])
%%
r = 0;
theta = 0;
z = r*exp(1j*theta);
% define displaced operator
deta = 0; % ∆ = 0.2
eipxi = 0; % ξ= pi/6
alpha = deta*exp(1j*eipxi);
% define 3-PPM modulation
K = 3;
% p & q
q = deta*cos(eipxi);
p = deta*sin(eipxi);

% Variance matrix
V = zeros(2,2);
V(1,1) = (cosh(r))^2+(sinh(r))^2+cos(theta)*sinh(2*r);
V(1,2) = sin(theta)*sinh(2*r);
V(2,1) = sin(theta)*sinh(2*r);
V(2,2) = (cosh(r))^2+(sinh(r))^2-cos(theta)*sinh(2*r);

% Wigner Function
x = -6:0.01:6;
y = -6:0.01:6;
[ X, Y ] = meshgrid( x, y );
W = (1/2*pi)*exp(-0.5*(V(2,2).*(X-q).^2+V(1,1).*(Y-p).^2-2*V(1,2)*(X-q).*(Y-p)));
mesh(x,y,W)
%figure(2)
subplot(2,1,2)
view([0,0,1])
xlabel('x')
ylabel('y')
zlabel('W(x,y)')







































































































function y = f(r,theta)
y = cosh(2*r)+tanh(r)*sinh(2*r)-(sinh(2*r)+tanh(r)*cosh(2*r)).*cos(theta);
end


