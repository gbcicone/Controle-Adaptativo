clc;
clear;
%Tempo:
t = linspace(0,25,1000);


%Parâmetros reais:
Km_real = 0.838;
La_real = 0.0036;
b_real = 0.268;
Ra_real = 1.36;
J_real = 0.001;


%Cálculo dos K:
K1_real = Km_real/(J_real*La_real);
K2_real = (J_real*Ra_real+b_real*La_real)/(J_real*La_real);
K3_real = (Km_real^2)/(J_real*La_real);


%Entrada: 

Va = sin(pi/3*t);


%Função transferência:

ft = tf(K1_real,[1 K2_real K3_real]);

%Saída:
w = lsim(ft, Va, t);

%Derivadas:

dw = gradient(w,t(2)-t(1));
ddw = gradient(dw, t(2)-t(1));

phi = [w dw ddw];

theta = zeros(3,1); % [1/K1; K2/K1; K3/K1]
theta_hist = zeros(3,length(t));
e_hist = zeros(length(t),1);

%Fator de esquecimento:

beta = 0.1;

ms2 = 2.95;


%Matriz de Covariância:

P = [80, 0, 0, 
    0, 10, 0, 
    0, 0, 10];
P_hist = zeros(3,3,length(t));

dT = 25/length(t);


for i=1:length(t)

    e = Va(i)-phi(i,:)*theta;
    dP = beta*P-P*(phi(i,:)'*phi(i,:)/ms2)*P;
    P = P+dP*dT;
    dtheta = P*e*phi(i,:)';
    theta = theta + dtheta*dT;

    e_hist(i) = e;
    theta_hist(:,i) = theta;
    P_hist(:,:,i) = P;

end

figure;
plot(t, theta_hist(3,:), 'r', 'DisplayName','1/K1 estimado',"LineWidth", 2); hold on;
plot(t, theta_hist(2,:), 'g', 'DisplayName','K2/K1 estimado',"LineWidth", 2);
plot(t, theta_hist(1,:), 'b', 'DisplayName','K3/K1 estimado',"LineWidth", 2);
yline(1/K1_real, '--r', 'DisplayName', '1/K1 real',"LineWidth", 2);
yline(K2_real/K1_real, '--g', 'DisplayName', 'K2/K1 real',"LineWidth", 2);
yline(K3_real/K1_real, '--b', 'DisplayName', 'K3/K1 real',"LineWidth", 2);

legend;
title('Estimativa dos parâmetros (RLS contínuo)');
xlabel('Tempo [s]');
ylabel('Valor estimado');
grid on;

figure;
plot(t, e_hist);
title('Erro de predição');
xlim([0, 24]);
xlabel('Tempo [s]');
ylabel('Erro');
grid on;

