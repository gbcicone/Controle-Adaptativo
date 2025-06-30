clc
clear 
close all

%% Parâmetros da Simulação
a = 1;
b = 5;
am = 6.7;
bm = 3.3;

k_star = (a + am) / b;
l_star = bm / b;

k0 = 1;
l0 = 1;

dt = 0.01;
t = 0:dt:7;
r = 15 * ones(1, length(t)); % Entrada constante r=15

%% Modelo de Referência
s = tf('s');
G = bm / (s + am);
xm = lsim(G, r, t);

%% MRAC Direto
gamma1 = 1.5; % Ganhos de adaptação para entrada constante
gamma2 = 1.5;

% Valores de x para o controlador adaptativo
x0 = 1;
x = zeros(1,length(t));
x(1) = x0;

x_hat_direto = zeros(1,length(t));
x0_hat = 0;
x_hat_direto(1) = x0_hat;

% Vetores para progressão dos ganhos e do erro
k_direto = zeros(1,length(t));
k_direto(1) = k0;
l_direto = zeros(1,length(t));
l_direto(1) = l0;
e_direto = zeros(1,length(t));
e_direto(1) = x_hat_direto(1)-xm(1);

u_direto = zeros(1,length(t));

e = zeros(1,length(t));
e(1) = x(1)-xm(1);

u_ideal = zeros(1,length(t));

for i=2:length(t)
    % X com controlador adaptativo
    e_direto(i-1) = x_hat_direto(i-1) - xm(i-1);
    
    k_direto_dot = gamma1*e_direto(i-1)*x_hat_direto(i-1)*sign(b);
    k_direto(i) = k_direto_dot*dt + k_direto(i-1);
    
    l_direto_dot = -gamma2*e_direto(i-1)*r(i-1)*sign(b);
    l_direto(i) = l_direto_dot*dt + l_direto(i-1);
    
    u_hat_direto = -k_direto(i)*x_hat_direto(i-1) + l_direto(i)*r(i-1);
    u_direto(i) = u_hat_direto;
    
    x_hat_direto_dot = a*x_hat_direto(i-1) + b*u_hat_direto;
    x_hat_direto(i) = x_hat_direto_dot * dt + x_hat_direto(i-1);
    
    % X como controlador ideal
    e(i-1) = x(i-1) - xm(i-1);
    u = -k_star*x(i-1) + l_star*r(i-1);
    u_ideal(i) = u;
    x_dot = a*x(i-1) + b*u;
    x_i = x(i-1) + x_dot*dt;
    x(i) = x_i;
end

%% MRAC Indireto

% Valores iniciais de adaptação
a0 = 1;
b0 = 1;
b_limite=1;

gamma1 = 2;
gamma2 = 2;

% Valores de x para o controlador adaptativo
x_hat_indireto = zeros(1,length(t));
x0_hat = 0;
x_hat_indireto(1) = x0_hat;

% Vetores para progressão dos ganhos e do erro
k_indireto = zeros(1,length(t));
k_indireto(1) = k0;
l_indireto = zeros(1,length(t));
l_indireto(1) = l0;
e_indireto = zeros(1,length(t));
e_indireto(1) = x_hat_indireto(1)-xm(1);

a_indireto = zeros(1,length(t));
a_indireto(1) = a0;
b_indireto = zeros(1,length(t));
b_indireto(1) = b0;

u_indireto = zeros(1,length(t));

for i=2:length(t)
    % X com controlador adaptativo
    e_indireto(i-1) = x_hat_indireto(i-1) - xm(i-1);
    
    k_indireto_i = (am+a_indireto(i-1))/b_indireto(i-1);
    k_indireto(i) = k_indireto_i;
    
    l_indireto_i = (bm)/b_indireto(i-1);
    l_indireto(i) = l_indireto_i;
    
    u_hat_indireto = -k_indireto(i)*x_hat_indireto(i-1) + l_indireto(i)*r(i-1);
    u_indireto(i) = u_hat_indireto;
    
    a_indireto_dot = gamma1*e_indireto(i-1)*x_hat_indireto(i-1);
    a_indireto(i) = a_indireto_dot*dt + a_indireto(i-1);
    
    b_indireto_dot = gamma2*e_indireto(i-1)*u_hat_indireto;
    b_indireto(i) = b_indireto_dot*dt + b_indireto(i-1);
    
    if ~ ( (abs(b_indireto(i)) > abs(b_limite)) || ( (abs(b_indireto(i)) > abs(b_limite)) && (e_indireto(i-1)*u_hat_indireto*sgn(b) > 0) ) )
        if b_indireto(i) < 0
            b_indireto(i) = - abs(b_limite);
        else
            b_indireto(i) = abs(b_limite);
        end
    end
    
    x_hat_indireto_dot = a*x_hat_indireto(i-1) + b*u_hat_indireto;
    x_hat_indireto(i) = x_hat_indireto_dot * dt + x_hat_indireto(i-1);
end

%% Plots

% Plot do controlador ideal

figure(1)
hold on
grid on
plot(t, x, 'b', "LineWidth", 3)
plot(t, xm, 'y--', "LineWidth", 3)
plot(t, e, 'r', "LineWidth", 2)
legend("Controlador Ideal", "Modelo de Referência", "Erro")
xlabel("Tempo [s]", "Fontsize", 15)
ylabel("Tensão", "Fontsize", 15)
title("Modelo de Referência ideal r = 15", "Fontsize", 15)


% Plot do controlador direto

figure(2)
hold on
grid on
plot(t, x_hat_direto,'b', "LineWidth", 3)
plot(t, xm, 'y--', "LineWidth", 3)
plot(t, e_direto, 'r', "LineWidth", 2)
legend("MRAC Direto", "Modelo de Referência", "Erro")
xlabel("Tempo [s]", "Fontsize", 15)
ylabel("Tensão", "Fontsize", 15)
title("MRAC Direto r = 15", "Fontsize", 15)


% Plot do controlador indireto
figure(3)
hold on
grid on
plot(t, x_hat_indireto, 'b', "LineWidth", 3)
plot(t, xm, 'y--', "LineWidth", 3)
plot(t, e_indireto, 'r', "LineWidth", 2)
legend("MRAC Indireto", "Modelo de Referência", "Erro")
xlabel("Tempo [s]", "Fontsize", 15)
ylabel("Tensão", "Fontsize", 15)
title("MRAC Indireto r = 15", "Fontsize", 15)

% Sinal de controle
figure(4)
hold on
grid on
plot(t, u_ideal, 'k', "LineWidth", 3)
plot(t, u_direto, 'c', "LineWidth", 2)
plot(t, u_indireto, 'm--', "LineWidth", 2) 
legend("Ideal","MRAC Direto","MRAC Indireto")
xlabel("Tempo [s]", "Fontsize", 15)
ylabel("Sinal u", "Fontsize", 15)
title("Sinal de Controle r = 15", "Fontsize", 15)

k_star_vector = k_star*ones(1, length(t));
figure(5)
hold on
grid on
plot(t, k_star_vector, 'k', "LineWidth", 3)
plot(t, k_direto, 'c', "LineWidth", 2)
plot(t, k_indireto, 'm', "LineWidth", 2)
legend("Ideal","MRAC Direto","MRAC Indireto")
xlabel("Tempo [s]", "Fontsize", 15)
ylabel("Parametro", "Fontsize", 15)
title("Parâmetro de controle k r = 15", "Fontsize", 15)

l_star_vector = l_star*ones(1, length(t));
figure(6)
hold on
grid on
plot(t, l_star_vector, 'k', "LineWidth", 3)
plot(t, l_direto, 'c', "LineWidth", 2)
plot(t, l_indireto, 'm', "LineWidth", 2)
legend("Ideal","MRAC Direto","MRAC Indireto")
xlabel("Tempo [s]", "Fontsize", 15)
ylabel("Parametro", "Fontsize", 15)
title("Parâmetro de controle l r = 15", "Fontsize", 15)

