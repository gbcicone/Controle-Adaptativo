


t = linspace(0,25,100);

u = 1 + cos(pi/3*t);

lambda = 100;

filtro = tf(1,[1,2,1]);

M = 100;

f = 0.15;

k = 7;

ft = tf(1,[M f k]);

%%

xa = lsim(ft, u, t);

za = lsim(filtro, u, t);

phi1a = lsim(tf([1 0 0],[1 2 1]),xa,t);

phi2a = lsim(tf([1 0],[1 2 1]),xa,t);

phi3a = lsim(tf(1,[1 2 1]),xa,t);

figure
plot(t,za)
legend('Z(t)')
grid on 
title('Z(t) exercício A')

figure
plot(t,phi1a, t, phi2a, t, phi3a)
legend('\phi_1', '\phi_2', '\phi_3')
grid on
title('Sinais \phi Exercício A')


%%



xb = lsim(ft,u,t);

zb = lsim(tf([M 0 0],[1 2 1]),xb,t ) - lsim(tf(1,[1 2 1]),u,t );


phi1b = lsim(-tf([1 0],[1 2 1]),xb,t); 

phi2b = lsim(-tf(1,[1 2 1]),xb,t);


figure
plot(t,zb)
legend('Z(t)')
grid on 
title('Z(t) exercício B')

figure
plot(t,phi1b, t, phi2b)
legend('\phi_1', '\phi_2')
grid on
title('Sinais \phi Exercício B')

