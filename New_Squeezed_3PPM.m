clc;
close all;
clear;
%% 两个挤压态，三维误码率
% define squeezed operator z0
r_0 = 0.1;
theta_0 = pi;
z_0 = r_0*exp(1j*theta_0);
% define squeezed operator z1
r_1 = 0.1;
theta_1 = pi;
z_1 = r_1*exp(1j*theta_1);
% define displaced operator α0
deta_0 = 0.1; % ∆ = 0.2,只影响对应的Ns
eipxi_0 = pi/6; % ξ= pi/6
alpha_0 = deta_0*exp(1j*eipxi_0);
Ns_0 = deta_0^2+(sinh(r_0))^2;
Ns_0 = 0.01:0.01:2.01;
% define displaced operator α1
deta_1 = 0.1; % ∆ = 0.2 ,只影响对应的Ns
eipxi_1 = pi/6; % ξ= pi/6
alpha_1 = deta_1*exp(1j*eipxi_1);
Ns_1 = deta_1^2+(sinh(r_1))^2;
Ns_1 = 0.01:0.01:2.01;
% define 3-PPM modulation
K = 3;
% 参数
mu_0 = cosh(r_0);
v_0 = sinh(r_0)*exp(1j*theta_0);
b0 = mu_0*alpha_0-v_0*alpha_0';
mu_1 = cosh(r_1);
v_1 = sinh(r_1)*exp(1j*theta_1);
b1 = mu_1*alpha_1-v_1*alpha_1';
A = mu_0*mu_1'-v_0*v_1';
B = v_0*mu_1-mu_0*v_1;
% calculate Pe
[ XX, YY ] = meshgrid( Ns_0, Ns_1 );
inner_pro = (A^(-1/2))*exp(-(A*(abs(b1)^2+abs(b0)^2)-2*b1*b0'+B*(b1')^2-B'*b0^2)/(2*A)).*sqrt((1/cosh(r_1)).*exp(-(YY-(sinh(r_1))^2).*f(r_1,theta_1))).*sqrt((1/cosh(r_0)).*exp(-(XX-(sinh(r_0))^2).*f(r_0,theta_0)));
P_e = 1-((sqrt(1+(K-1).*inner_pro)+(K-1).*sqrt(1-inner_pro)).^2)/K^2;
figure(1)
mesh(Ns_0, Ns_1, P_e); %'EdgeColor','b'
set(gca,'zscale','log')%将z轴上刻度单位设置为对数坐标型
xlabel(' Ns_0[photons/symbol]')
ylabel(' Ns_1[photons/symbol]')
zlabel('Error Probability Pe')
title('Two squeezed states in 3-PPM')


%% 固定Ns_1=0.5 1.0 1.5 2.0，与一个挤压态进行比较，的二维曲线图
XX = 0.5; % Ns_0 = 2
inner_pro = (A^(-1/2))*exp(-(A*(abs(b1)^2+abs(b0)^2)-2*b1*b0'+B*(b1')^2-B'*b0^2)/(2*A)).*sqrt((1/cosh(r_1)).*exp(-(YY-(sinh(r_1))^2).*f(r_1,theta_1))).*sqrt((1/cosh(r_0)).*exp(-(XX-(sinh(r_0))^2).*f(r_0,theta_0)));
P_e = 1-((sqrt(1+(K-1).*inner_pro)+(K-1).*sqrt(1-inner_pro)).^2)/K^2;
figure(2)
semilogy(Ns_1,P_e(:,1),'k:+')

hold on



XX = 1.0; % Ns_0 = 2
inner_pro = (A^(-1/2))*exp(-(A*(abs(b1)^2+abs(b0)^2)-2*b1*b0'+B*(b1')^2-B'*b0^2)/(2*A)).*sqrt((1/cosh(r_1)).*exp(-(YY-(sinh(r_1))^2).*f(r_1,theta_1))).*sqrt((1/cosh(r_0)).*exp(-(XX-(sinh(r_0))^2).*f(r_0,theta_0)));
P_e = 1-((sqrt(1+(K-1).*inner_pro)+(K-1).*sqrt(1-inner_pro)).^2)/K^2;

plot(Ns_1,P_e(:,1),'k--p')
hold on

XX = 1.5; % Ns_0 = 2
inner_pro = (A^(-1/2))*exp(-(A*(abs(b1)^2+abs(b0)^2)-2*b1*b0'+B*(b1')^2-B'*b0^2)/(2*A)).*sqrt((1/cosh(r_1)).*exp(-(YY-(sinh(r_1))^2).*f(r_1,theta_1))).*sqrt((1/cosh(r_0)).*exp(-(XX-(sinh(r_0))^2).*f(r_0,theta_0)));
P_e = 1-((sqrt(1+(K-1).*inner_pro)+(K-1).*sqrt(1-inner_pro)).^2)/K^2;

plot(Ns_1,P_e(:,1),'k-*')
hold on

XX = 2.0; % Ns_0 = 2
inner_pro = (A^(-1/2))*exp(-(A*(abs(b1)^2+abs(b0)^2)-2*b1*b0'+B*(b1')^2-B'*b0^2)/(2*A)).*sqrt((1/cosh(r_1)).*exp(-(YY-(sinh(r_1))^2).*f(r_1,theta_1))).*sqrt((1/cosh(r_0)).*exp(-(XX-(sinh(r_0))^2).*f(r_0,theta_0)));
P_e = 1-((sqrt(1+(K-1).*inner_pro)+(K-1).*sqrt(1-inner_pro)).^2)/K^2;

plot(Ns_1,P_e(:,1),'k:.')
hold on

%
Ns = 0.01:0.005:2;
Nr = Ns/log2(K);
theta = pi;
r = 0.1;
inner_pro = (1/cosh(r))*exp(-(Ns-(sinh(r))^2)*f(r,theta));
Pe = 1-((sqrt(1+(K-1).*inner_pro)+(K-1).*sqrt(1-inner_pro)).^2)/K^2;
plot(Ns,Pe,'k-')
xlabel('Received photons Ns_1[photons/symbol]')
ylabel('Error Probability Pe')
legend('Two squeezed-displaced states,Ns_0 =0.5','Two squeezed-displaced states,Ns_0 =1.0','Two squeezed-displaced states,Ns_0 =1.5','Two squeezed-displaced states,Ns_0 =2.0','One squeezed-displaced state,r=0.1,θ=π')
axis([0 2.01 0.001 1]) 
%

function y = f(r,theta)
y = cosh(2*r)+tanh(r)*sinh(2*r)-(sinh(2*r)+tanh(r)*cosh(2*r)).*cos(theta);
end


%%


