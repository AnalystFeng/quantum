clc;clear;
close all;
N = 0.001:0.001:0.1;

Pc_end = zeros(3,length(N));
alpha = 4; % 信号光子数 = 4
r = 0.5;
n = 10;
theta = pi;
ket = zeros(n,1);
noise = zeros(n,n);

%%  r = 0.6
mu = cosh (r);
v = sinh(r)*exp(1*i*theta);
B = mu*alpha-v*alpha';
for i = 1:1:10
    ket(i) = (sqrt(factorial(i))/mu) * (B/mu)^i * hermiteH(i,mu*v/(B^2)) * exp(  -0.5*(abs(B))^2 - B^2*(v'/2*mu));
end
rou = ket*ket';
% rou = rou+noise;
rou = rou/trace(rou);  % rou 带噪声挤压态的密度算子  ； 
R_deta = rou;
for k = 1:1:10
    k
% 定义 noise
    for i =1:1:n
        noise(i,i) = N(k)^i/(1+N(k))^(1+i);
    end
    noise = noise*(1/trace(noise));  % 噪声算子              
    R_0 = noise;
    h = 1;
    
    [A,E] = eig(R_deta);
    [C,D] = eig(R_0);
    rou_up_1 = A(:,end-h+1:end)*E(end-h+1:end,end-h+1:end)*A(:,end-h+1:end)';
    rou_up_0 = C(:,end-h+1:end)*D(end-h+1:end,end-h+1:end)*C(:,end-h+1:end)';
    rou_down_0 = kron(rou_up_0,kron(rou_up_0,rou_up_1));

    gamma_up_1 = A(:,end-h+1:end)*sqrt(E(end-h+1:end,end-h+1:end));
    gamma_up_0 = C(:,end-h+1:end)*sqrt(D(end-h+1:end,end-h+1:end));
    gamma_down_0 = kron(gamma_up_0,kron(gamma_up_0,gamma_up_1));
    gamma_down_1 = kron(gamma_up_0,kron(gamma_up_1,gamma_up_0));
    gamma_down_2 = kron(gamma_up_1,kron(gamma_up_0,gamma_up_1));
    States_matrix = [gamma_down_0 gamma_down_1 gamma_down_2];
    T = States_matrix*States_matrix';
    % G = States_matrix'*States_matrix;
    % S = (rou_down_0*T^(-1/2))^2;
    Pc = trace(  (rou_down_0*     pinv(T^(0.5))    )^2           );
    Pc_end(1,k) = Pc;
        
end
    
subplot(2,1,1)
plot(N',1-real(Pc_end(1,:)))

