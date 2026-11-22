close all
clear 
% Definizione del vettore dei parametri
params = [    690;   % m_s
              1222;   % I_psi
              1.3;    % L_f
              1.5;    % L_r
              18000;  % k_sf
              1000;   % c_sf
              22000;  % k_sr
              1000;   % c_sr
              40;     % m_uf
              45;     % m_ur
              200000; % k_uf
              200000];% k_ur


x0 = zeros(8,1);  % Stato iniziale [z_s, z_s_dot, psi, psi_dot, z_uf, z_uf_dot, z_ur, z_ur_dot]
% x0= [0; pi/12; 0; 0; 0.04; 0; 0.04; 0];

t = 0:0.01:10;             % Tempo di simulazione
u1 = zeros(size(t));        % Segnale costante
u2 = zeros(size(t));        % Segnale costante
u = [t' u1' u2'];          % Crea una matrice [tempo, ingresso1, ingresso2]


[disturbances,z_rf,z_rr]=strada_dissestata1();


% Simulazione
result=sim("non_linear_model.slx");

figure()
% Plot della prima uscita
plot(result.y.Time, squeeze(result.y.Data(1,1,:)),'LineWidth',1.5,'Color','black');
xlabel('Tempo [s]');
ylabel('Spostamento [m]');
title('Z_s');
hold on
plot(result.y.Time,z_rf)
hold on
plot(result.y.Time,z_rr)
hold off
grid on;

figure()
% Plot della prima uscita
plot(result.y.Time, squeeze(result.y.Data(2,1,:)),'LineWidth',1.5);
xlabel('Tempo [s]');
ylabel('Spostamento [m]');
title('Z_s');
grid on;

figure()
% Plot della prima uscita
plot(result.y.Time, squeeze(result.y.Data(3,1,:)),'LineWidth',1.5);
xlabel('Tempo [s]');
ylabel('Spostamento [m]');
title('Z_s');
grid on;

figure()
% Plot della prima uscita
plot(result.y.Time, squeeze(result.y.Data(4,1,:)),'LineWidth',1.5);
xlabel('Tempo [s]');
ylabel('Spostamento [m]');
title('Z_s');
grid on;

figure()
% Plot della prima uscita
plot(result.y.Time, squeeze(result.y.Data(5,1,:)),'LineWidth',1.5);
xlabel('Tempo [s]');
ylabel('Spostamento [m]');
title('Z_s');
grid on;

figure()
% Plot della prima uscita
plot(result.y.Time, squeeze(result.y.Data(6,1,:)),'LineWidth',1.5);
xlabel('Tempo [s]');
ylabel('Spostamento [m]');
title('Z_s');
grid on;

figure()
% Plot della prima uscita
plot(result.y.Time, squeeze(result.y.Data(7,1,:)),'LineWidth',1.5);
xlabel('Tempo [s]');
ylabel('Spostamento [m]');
title('Z_s');
grid on;

figure()
% Plot della prima uscita
plot(result.y.Time, squeeze(result.y.Data(8,1,:)),'LineWidth',1.5);
xlabel('Tempo [s]');
ylabel('Spostamento [m]');
title('Z_s');
grid on;







