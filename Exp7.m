%% Problema 1



clc; clear;





%Vetor de tempo

t = linspace(0,10,1000);





%Vetor de setpoint

r = 5*cos(t)+2*sin(pi*t);





% Parâmetros das equações reais:

a = 2;

b = 5;



% Parâmetros das equações estimadas:

am = 10;

bm = 1;



%Criação dos vetores do laço:

x = zeros(1, length(t));

k = zeros(1, length(t));

l = zeros(1, length(t));

e = zeros(1, length(t));

% Ganhos reais:

k_real = (am-a)/b;
l_real = bm/b;


% Xm estimado:

ft1 = tf(bm,[1 am]);



xm = lsim(ft1, r, t);



% Condições iniciais:

x(1) = 0;

xm(1) = 0;

k(1) = 0;

l(1) = 0;



% Gamma:

gamma = 40;



% Intervalo de tempo:

dt = t(2) - t(1);





for i=1:length(t)-1



    % Ação de controle (saída)

    u = -k(i)*r(i) + l(i)*x(i);



    % Encontrando o X:

    x(i+1) = x(i) + dt * (a*x(i) + b*u);



    % Cálculo do erro:

    e(i) = x(i) - xm(i);



    % Lei de adaptação:

    dk = gamma * e(i) * r(i);

    dl = -gamma * e(i) * x(i);



    k(i+1) = k(i) + dt * dk;

    l(i+1) = l(i) + dt * dl;





end



%Plote:



figure;

subplot(2,2,1);


plot(t, x, 'b', t, xm, 'r--');

legend('x(t)', 'x_m(t)');

title('Saídas da Planta e do Modelo');
grid on



subplot(2,2,2);

plot(t, e);

title('Erro e(t) = x(t) - x_m(t)');
grid on



subplot(2,2,3);

plot(t, k);
yline(k_real, '--b', 'k*');

title('Ganho adaptativo k(t)');
grid on



subplot(2,2,4);

plot(t, l);
yline(l_real,'--or', 'l*');

title('Ganho adaptativo l(t)');
grid on



%% Problema 2



clc;clear;



% Vetor de tempo:

t = linspace(0,5,1000);



% Seletor de casos de setpoint:

casos = 3;



if casos == 1

    r =5*ones(1,length(t));

elseif casos == 2

    r = sin(2*t);

elseif casos == 3

    r = 1./(1+t);

end



% Equação estimada:

ft = tf(2,[1 2]);



ym = lsim(ft, r, t);



% Parâmetros da equação real:

a = -1;



b = 12;


% Ganhos reais:

k_real = (2-a)/b;
l_real = 2/b;


% Iniciando vetores do laço:

yp = zeros(1, length(t));

k = zeros(1, length(t));

l = zeros(1, length(t));

e = zeros(1, length(t));



% Intervalo de tempo:

dt = t(2)-t(1);





% Gamma:

gamma = 10;





for i=1:length(t)-1



    % Ação de controle (sáida)

    u = -k(i)*r(i) + l(i)*yp(i);



    % Encontrando Yp:

    yp(i+1) = yp(i) + dt * (a*yp(i) + b*u);



    % Erro:

    e(i) = yp(i) - ym(i);



    % Lei de adaptação:

    dk = gamma * e(i) * r(i);

    dl = -gamma * e(i) * yp(i);



    k(i+1) = k(i) + dt * dk;

    l(i+1) = l(i) + dt * dl;







end



% Plote:



figure;

subplot(2,2,1);

plot(t, yp, 'b', t, ym, 'r--');

legend('yp(t)', 'y_m(t)');

title('Saídas da Planta e do Modelo');
grid on


subplot(2,2,2);

plot(t, e);

title('Erro e(t) = yp(t) - y_m(t)');

grid on

subplot(2,2,3);

plot(t, k);
yline(k_real, '--b', 'k*');
title('Ganho adaptativo k(t)');
grid on


subplot(2,2,4);

plot(t, l);
yline(l_real,'--or', 'l*');
title('Ganho adaptativo l(t)');



grid on





%% Problema 3

clc;clear;

% Vetor de tempo:


t = linspace(0,10,1000);


% Velocidade desejada:


Vs = 55*ones(1,length(t)); 


ft = tf(0.5,[1 0.5]);


Vm = lsim(ft,Vs,t);


% Iniciando vetores do laço:

V = zeros(1, length(t));

k = zeros(1, length(t));

l = zeros(1, length(t));

e = zeros(1, length(t));

delta = zeros(1, length(t));

% Casos dos parâmetros reais:

caso = 2;


if caso == 1

    a =0.5;

    b = 1.5;

    d= 10;

elseif caso == 2

    a = 0;

    b = 1.5;

    d = 0.2 + sin(0.02*t);

end


dt = t(2) - t(1);


gamma1 = 1;
gamma2 = 1;
gamma3 = 1;

%k(1) = 1;
%l(1) = 1;
%delta(1) = 1;



for i=1:length(t)-1

    % Ação de controle (sáida)

    theta = -k(i)*V(i) + l(i)*Vs(i)+delta(i);



    % Encontrando Yp:

    
    if caso == 2
        V(i+1) = V(i) + dt * (-(0.5 + 0.05/(1+V(i)))*V(i) + b*theta+d(i));
    elseif caso == 1
        V(i+1) = V(i) + dt * (-a*V(i) + b*theta+d);
    end
    


    % Erro:

    e(i) = V(i) - Vm(i);



    % Lei de adaptação:

    dk = gamma1 * e(i) * V(i)*sign(b);

    dl = -gamma2 * e(i) * Vs(i)*sign(b);

    ddelta = -gamma3 * e(i)*sign(b);


    k(i+1) = k(i) + dt * dk;

    l(i+1) = l(i) + dt * dl;

    delta(i+1) = delta(i)+dt*ddelta;


end


figure;

subplot(3,1,1); plot(t, V, t, Vm, '--'); legend('V(t)', 'Vm(t)');

title('Saída do sistema e modelo de referência'); grid on;



subplot(3,1,2); plot(t, k, t, l, t, delta);

legend('k(t)', 'l(t)', '\delta(t)'); title('Parâmetros adaptativos'); grid on;



subplot(3,1,3); plot(t, e); title('Erro e(t)'); grid on;




%% Problema 3 Método indireto

clc;clear;

% Vetor de tempo:


t = linspace(0,10,1000);


% Velocidade desejada:


Vs = 55*ones(1,length(t)); 


ft = tf(0.5,[1 0.5]);

am=0.5;
bm = 0.5;

Vm = lsim(ft,Vs,t);


% Iniciando vetores do laço:

V = zeros(1, length(t));


e = zeros(1, length(t));

delta = zeros(1, length(t));

a_hat = zeros(1,length(t));

b_hat =zeros(1,length(t));

d_hat = zeros(1,length(t));

% Casos dos parâmetros reais:

caso = 1;


if caso == 1

    a =0.5;

    b = 1.5;

    d= 10;

elseif caso == 2

    a = 0;

    b = 1.5;

    d = 0.2 + sin(0.02*t);

end


dt = t(2) - t(1);


gamma1 = 10;
gamma2 = 10;
gamma3 = 10;

%k(1) = 1;
%l(1) = 1;
%delta(1) = 1;



for i=1:length(t)-1

    % Ação de controle (sáida)

    theta = -((am+a_hat(i))/b_hat(i))*V(i) + (bm/b_hat(i))*Vs(i)+(a_hat(i)*Vs(i)-d_hat(i))/(b_hat(i));



    % Encontrando Yp:

    
    if caso == 2
        V(i+1) = V(i) + dt * (-(0.5 + 0.05/(1+V(i)))*V(i) + b*theta+d(i));
    elseif caso == 1
        V(i+1) = V(i) + dt * (-a*V(i) + b*theta+d);
    end
    


    % Erro:

    e(i) = V(i) - Vm(i);



    % Lei de adaptação:

    da_hat = gamma1 * e(i) * V(i)*sign(b);

    if abs(b_hat(i)) > b || (abs(b_hat(i)) == b && e(i)*theta*sign(b)>= 0)
        db_hat = -gamma2 * e(i) * theta*sign(b);
    else
        db_hat = 0;
    end
    dd_hat = -gamma3 * e(i)*sign(b);
    

    a_hat(i+1) = a_hat(i) + dt * da_hat;

    b_hat(i+1) = b_hat(i) + dt * db_hat;

    d_hat(i+1) = d_hat(i)+dt*dd_hat;


end


figure;

subplot(3,1,1); plot(t, V, t, Vm, '--'); legend('V(t)', 'Vm(t)');

title('Saída do sistema e modelo de referência'); grid on;



subplot(3,1,2); plot(t, a_hat, t, b_hat, t, d_hat);

legend('hat{a}(t)', 'hat{b}(t', 'hat{d}(t'); title('Parâmetros adaptativos'); grid on;



subplot(3,1,3); plot(t, e); title('Erro e(t)'); grid on;

%% Problema 8 a) Método LS


clc;clear;

% Vetor de tempo:


t = linspace(0,5,1000);
dt = t(2)-t(1);


% Entrada:

r = 1*ones(length(t),1);


% Parâmetros planta referência:

am = 10;
bm = 10;

% Parâmetros reais:

a = -5;
b =  3;

% Valores ideais de controle
k_real = (am-a)/b;
l_real = bm/b;

% Planta referência:
ft = tf(bm, [1 am]);

ym = lsim(ft, r, t);

gamma1 = 20;

gamma2 = 20;

% Iniciando vetores do laço:

u = zeros(1, length(t));

y = zeros(1, length(t));

k = zeros(1, length(t));

l = zeros(1, length(t));

e = zeros(1, length(t));

ms_2 = zeros(1, length(t));

for i=1:length(t)-1

    % Ação de controle (sáida)

    u = -k(i)*y(i) + l(i)*r(i);


    y(i+1) = y(i) + dt * (-a*y(i) + b*u);


    % Erro:

    ms_2(i) = 1 + y(i)^2;

    e(i) = (y(i) - ym(i))/ms_2(i);

    % Lei de adaptação:

    dk = gamma1 * e(i) * y(i)*sign(b);

    dl = -gamma2 * e(i) * r(i)*sign(b);

    


    k(i+1) = k(i) + dt * dk;

    l(i+1) = l(i) + dt * dl;

    

end


figure;

subplot(3,1,1); plot(t, y, t, ym, '--'); legend('y(t)', 'ym(t)');

title('Saída do sistema e modelo de referência'); grid on;



subplot(3,1,2); plot(t, k, t, l);
yline(l_real,'--or', 'l*');
yline(k_real, '--b', 'k*');
legend('k(t)', 'l(t)'); title('Parâmetros adaptativos'); grid on;



subplot(3,1,3); plot(t, e); title('Erro e(t)'); grid on;



%% Problema 8 a) Método Direto Gradient


clc;clear;

% Vetor de tempo:


t = linspace(0,5,1000);
dt = t(2)-t(1);


% Entrada:

r = 1*ones(length(t),1);


% Parâmetros planta referência:

am = 10;
bm = 10;

% Parâmetros reais:

a = -5;
b =  3;

% Valores ideais de controle
k_real = (am-a)/b;
l_real = bm/b;

% Planta referência:
ft = tf(bm, [1 am]);

ym = lsim(ft, r, t);

dym = gradient(ym,dt);

phi = [ym dym];

gamma1 = 20;

gamma2 = 20;

% Iniciando vetores do laço:

u = zeros(1, length(t));

uf = zeros(1, length(t));

theta = zeros(2,length(t));

epslon = zeros(1,length(t));

b_hat = zeros(1,length(t));

phi = zeros(2,length(t));

xi = zeros(1,length(t));

y = zeros(1, length(t));

k = zeros(1, length(t));

l = zeros(1, length(t));

e = zeros(1, length(t));

ms_2 = zeros(1, length(t));

P = 100*eye(2);


gamma = 1;

for i=1:length(t)-1

    % Ação de controle (sáida)

    u = -k(i)*y(i) + l(i)*r(i);


    uf_dot = -am*uf(i) + u;

    uf(i+1) = uf_dot*dt + uf(i);


    phi_dot_1 = -am*phi(1,i) + y(i);
    phi_dot_2 = -am*phi(2,i) - r(i);

    phi(1,i+1) = phi_dot_1*dt + phi(1,i);
    phi(2,i+1) = phi_dot_2*dt + phi(2,i);

    xi(i+1) = theta(:,i)'*phi(:,i+1) + uf(i+1);

    db_hat = gamma*epslon(i)*xi(i);
    b_hat(i+1) = db_hat*dt + b_hat(i);

    

    % Erro:

    e(i) = y(i) - ym(i);
    ms_2(i) = 1 + phi(:,i)'*phi(:,i) + uf(i)^2;
    epslon(i) = (e(i)-b_hat(i)*xi(i))/ms_2(i);

    % Lei de adaptação:


    dtheta = P*epslon(i)*phi(:,i);
    theta(:,i+1) = dtheta*dt + theta(:,i);


    k(i+1) = theta(1,i);
    l(i+1) = theta(2,i);
    

    y(i+1) = y(i) + dt * (-a*y(i) + b*u);

end


figure;

subplot(3,1,1); plot(t, y, t, ym, '--'); legend('y(t)', 'ym(t)');

title('Saída do sistema e modelo de referência'); grid on;



subplot(3,1,2); plot(t, k, t, l);
yline(l_real,'--or', 'l*');
yline(k_real, '--b', 'k*');
legend('k(t)', 'l(t)'); title('Parâmetros adaptativos'); grid on;



subplot(3,1,3); plot(t, e); title('Erro e(t)'); grid on;


