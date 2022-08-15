clc;
clear;
close all;
%% 可以验证测量结果 m 为整数 k 的概率是由几何分布给出的
% N = 0.3;
% for k = 1:30
%     p = (N^k)/((N+1)^k);
%     stem(k,p);
%     hold on
% end
%% 积分是关于复变量 \alpha 的，其中的被积函数具有二维高斯分布
% alpha_real = -9:0.2:9;
% alpha_Im = -9:0.2:9;
% N = 3;
% [ XX, YY ] = meshgrid( alpha_real, alpha_Im );
% alpha = XX+1i*YY;
% f = zeros(length(alpha_real),length(alpha_Im));
% gamma = 3+3*i;
% for i = 1:length(alpha_real)*length(alpha_Im) 
%   f(i) = exp(-    (      abs(   alpha(i)-gamma )      )     .^2 / N);
% end
% mesh(alpha_real,alpha_Im,f)

%% 3-PPM GUS SRM
n = 20;
h = 5;
N = 0.1;
v = N/(1+N);
deta = sqrt(5);
R_0 = zeros(n,n);   % 初始化R_0
R_deta = zeros(n,n);  % 初始化R_deta

% 定义 R_0
for i =1:1:n
    R_0(i,i) = (1-v)*v^i;
end
R_0 = R_0*(1/trace(R_0));
% 定义 R_deta 其为酋矩阵，迹为1
for i = 1:1:n       % i 行号
    for j = 1:1:n       % j 列号
        if j >= i
           R_deta(i,j) = (1-v)*v^j * sqrt(factorial(i)/factorial(j)) * (deta/N)^(j-i) * exp(-(1-v)*deta^2) * mlaguerre(i,j-i,-deta^2/(N*(N+1))) ;
        end
    end
end
for i = 1:1:n       % i 行号
    for j = 1:1:n       % j 列号
        if j >= i
           R_deta(j,i) = R_deta(i,j);
        end
    end
end
%%
[A,B] = eig(R_deta);
[C,D] = eig(R_0);
rou_up_1 = A(:,end-h+1:end)*B(end-h+1:end,end-h+1:end)*A(:,end-h+1:end)';
rou_up_0 = C(:,end-h+1:end)*D(end-h+1:end,end-h+1:end)*C(:,end-h+1:end)';
rou_down_0 = kron(rou_up_0,kron(rou_up_0,rou_up_1));

gamma_up_1 = A(:,end-h+1:end)*sqrt(B(end-h+1:end,end-h+1:end));
gamma_up_0 = C(:,end-h+1:end)*sqrt(D(end-h+1:end,end-h+1:end));
gamma_down_0 = kron(gamma_up_0,kron(gamma_up_0,gamma_up_1));
gamma_down_1 = kron(gamma_up_0,kron(gamma_up_1,gamma_up_0));
gamma_down_2 = kron(gamma_up_1,kron(gamma_up_0,gamma_up_1));
States_matrix = [gamma_down_0 gamma_down_1 gamma_down_2];
T = States_matrix*States_matrix';
G = States_matrix'*States_matrix;
S = (rou_down_0*T^(-1/2))^2;
Pc = trace(  (rou_down_0*     pinv(T^(0.5))    )^2           );





% 星座
% rou_0 = kron(kron(R_0,R_0),R_deta);
% rou_1 = kron(kron(R_0,R_deta),R_0);
% rou_2 = kron(R_deta,kron(R_0,R_0));
% StateM = [rou_0 rou_1 rou_2];
% T = StateM'*StateM;
% Pc = trace((rou_0*T^(-1/2))^2)


function L = mlaguerre(n,p,x)
 ret = 0;
 for i=0:n
     sum1 = ret + (power(-1,i)* (factorial(n+p)/(factorial(n-i) * factorial(p+i) * factorial(i)))*power(x,i));
     ret=sum1;
 end
 L=ret;
end

function y = L(m,n,gamma)
N = 0.2;
y = laguerreL(m,n-m,-gamma^2/(N*(N+1))); 
end