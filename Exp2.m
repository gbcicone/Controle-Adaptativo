%% 2
clc
clear
close all

K1 = 1;
K2 = 12;
K3 = 100;

s = tf([1,0], 1);

t = linspace(0,25,100);

r = 2.7*cos(4*10e-7*t);

ft = tf(1,[K1 K2 K3]);

y = lsim(ft, r, t);

%% a

lambda = 1;
filtro = 1/(s+lambda)^2;

z = lsim(filtro, r, t);

% theta = [K1 K2 K3];

phi1 = lsim(s^2*filtro, y, t);
phi2 = lsim(s*filtro, y, t);
phi3 = lsim(1*filtro, y, t);

figure(1);
plot(t,z);
xlabel('Tempo [s]');
title('z(t) - exercício A');
grid on;

figure(2);
plot(t,phi1, t,phi2, t,phi3);
xlabel('Tempo [s]');
legend('\phi1', '\phi2', '\phi3');
title('Sinais \phi - exercício A');
grid on;

fig = 3;
%% a - diferentes lambdas

% lambdas = linspace(1,5,5);
% 
% legend_entries = {};
% fig = 7;
% for i=1:length(lambdas)
% 
%     lambda = lambdas(i);
%     filtro = 1/(s+lambda)^2;
% 
%     z = lsim(filtro, r, t);
% 
%     % theta = [K1 K2 K3];
% 
%     phi1 = lsim(s^2*filtro, y, t);
%     phi2 = lsim(s*filtro, y, t);
%     phi3 = lsim(1*filtro, y, t);
% 
%     legend_entries{length(legend_entries)+1} = sprintf('\\lambda = %.2f', lambda);
% 
%     figure(3);
%     plot(t,z);
%     xlabel('Tempo [s]');
%     title('z');
%     grid on;
%     hold on;
%     legend(legend_entries);
% 
%     figure(4);
%     plot(t,phi1);
%     xlabel('Tempo [s]');
%     title('\phi1');
%     grid on;
%     hold on;
%     legend(legend_entries);
% 
%     figure(5);
%     plot(t,phi2);
%     xlabel('Tempo [s]');
%     title('\phi2');
%     grid on;
%     hold on;
%     legend(legend_entries);
% 
%     figure(6);
%     plot(t,phi3);
%     xlabel('Tempo [s]');
%     title('\phi3');
%     grid on;
%     hold on;
%     legend(legend_entries);
% 
%     figure(fig);
%     plot(t,z);
%     xlabel('Tempo [s]');
%     title(sprintf('z - \\lambda = %.2f', lambda));
%     grid on;
%     hold on;
% 
%     figure(fig+1);
%     plot(t,phi1, t,phi2, t,phi3);
%     xlabel('Tempo [s]');
%     title(sprintf('\\phi - \\lambda = %.2f', lambda));
%     grid on;
%     hold on;
%     legend('\phi1','\phi2','\phi3');
% 
%     fig = fig+2;
% end

%% b
lambda = 1;
filtro = 1/(s+lambda)^2;

z = lsim(filtro, y, t);

% theta = [K1/K3 K2/K3 1/K3];

phi1 = lsim(-s^2*filtro, y, t);
phi2 = lsim(-s*filtro, y, t);
phi3 = lsim(1*filtro, r, t);

figure(fig);
plot(t,z);
xlabel('Tempo [s]');
title('z(t) - exercício B');
grid on;

figure(fig+1);
plot(t,phi1, t,phi2, t,phi3);
xlabel('Tempo [s]');
legend('\phi1', '\phi2', '\phi3');
title('Sinais \phi - exercício B');
grid on;

fig = fig+2;
%% c
lambda = 1;
filtro = 1/(s+lambda)^2;

z = lsim(K2*s*filtro, y, t) - lsim(filtro, r, t);

% theta = [K1 K3];

phi1 = lsim(-s^2*filtro, y, t);
phi2 = lsim(-1*filtro, y, t);

figure(fig);
plot(t,z);
xlabel('Tempo [s]');
title('z(t) - exercício C');
grid on;

figure(fig+1);
plot(t,phi1, t,phi2);
xlabel('Tempo [s]');
legend('\phi1', '\phi2');
title('Sinais \phi - exercício C');
grid on;