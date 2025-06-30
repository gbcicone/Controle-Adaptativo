%% Código GPT
clc;
clear;

% Parâmetros reais
M_real = 100;
f_real = 0.15;
k_real = 7;

% Tempo e entrada
t = linspace(0,25,1000);
u = 1 + cos(pi/3*t);

% Sistema original
ft = tf(1,[M_real f_real k_real]);
x = lsim(ft, u, t);
dx = gradient(x, t(2)-t(1));
ddx = gradient(dx, t(2)-t(1));

% Vetores para identificação
phi = [x(:), dx(:), ddx(:)]; % [x, dx, ddx] em colunas
Y = u';

% Inicialização
theta = zeros(3,1); % [k; f; M]
P = 1000*eye(3);
lambda = 0.99;

% Armazenar estimativas
theta_hist = zeros(length(t), 3);

% Loop RLS
for i = 1:length(t)
    phi_i = phi(i,:)';
    y_i = Y(i);
    K = (P * phi_i) / (lambda + phi_i' * P * phi_i);
    theta = theta + K * (Y(i) - phi_i' * theta);
    P = (1/lambda) * (P - K * phi_i' * P);
    theta_hist(i,:) = theta';
end

% Plot das estimativas
figure;
plot(t, theta_hist(:,1), 'r', t, theta_hist(:,2), 'g', t, theta_hist(:,3), 'b');
legend('k estimado', 'f estimado', 'M estimado');
title('Estimativas dos parâmetros via RLS');
grid on;



%% Código autoral:


clc;clear;


% Parâmetros reais
M_real = 100;
f_real = 0.15;
k_real = 7;

% Tempo:

t = linspace(0,25,1000);

% Entrada:

u = 1 + cos(pi/3*t);
ft = tf(1,[M_real f_real k_real]);

[x, t_return, x_pp] = lsim(ss(ft), u, t);

dx = gradient(x,t(2)-t(1));
ddx = gradient(dx,t(2)-t(1));

x_p = x_pp(:,1);
x_pp = x_pp(:,2);

%phi = [x x_p x_pp];
phi = [x dx ddx];


theta = zeros(3,1); % [k; f; M]
theta_hist = zeros(3,length(t));


beta = 1;

%ms2 = 1;

ms2 = zeros(1,length(t));

%Matriz de Covariância:

P = 100*eye(3);
P_hist = zeros(3,3,length(t));

dT = 25/length(t);

for i=1:length(t)
    ms2(i) = 1 + (phi(i,1))^2;
    e = (u(i)-phi(i,:)*theta)/ms2(i);
    dP = beta*P-P*(phi(i,:)'*phi(i,:)/ms2(i))*P;
    P = P+dP*dT;
    dtheta = P*e*phi(i,:)';
    theta = theta + dtheta*dT;

    
    theta_hist(:,i) = theta;
    P_hist(:,:,i) = P;
end



figure;
plot(t, theta_hist(1,:), 'r', 'DisplayName','k estimado'); hold on;
plot(t, theta_hist(2,:), 'g', 'DisplayName','f estimado');
plot(t, theta_hist(3,:), 'b', 'DisplayName','M estimado');

yline(k_real, 'r--', 'DisplayName','K');
yline(f_real, 'g--', 'DisplayName','f');
yline(M_real, 'b--', 'DisplayName','M');

legend;
title('Estimativa dos parâmetros (RLS contínuo)');
xlabel('Tempo [s]');
ylabel('Valor estimado');
grid on;



