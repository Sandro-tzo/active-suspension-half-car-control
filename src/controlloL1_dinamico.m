clear; close all;

% Parametri del sistema
m_s = 690;    % Massa del corpo principale [kg]
I_psi = 1222;  % Momento d'inerzia rispetto all'asse [kg*m^2]
m_uf = 40;     % Massa della ruota anteriore [kg]
m_ur = 45;     % Massa della ruota posteriore [kg]
k_sf = 18000;  % Costante elastica sospensione anteriore [N/m]
k_sr = 22000;  % Costante elastica sospensione posteriore [N/m]
c_sf = 1000;   % Costante smorzamento sospensione anteriore [Ns/m]
c_sr = 1000;   % Costante smorzamento sospensione posteriore [Ns/m]
k_uf = 200000; % Costante elastica ruota anteriore [N/m]
k_ur = 200000; % Costante elastica ruota posteriore [N/m]
L_f = 1.3;     % Distanza tra CG e asse anteriore [m]
L_r = 1.5;     % Distanza tra CG e asse posteriore [m]

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

% Matrice A (dinamica del sistema)
A = [
    0, 1, 0, 0, 0, 0, 0, 0;
    -(k_sf + k_sr)/m_s, -(c_sf + c_sr)/m_s, (k_sf*L_f - k_sr*L_r)/m_s, (c_sf*L_f - c_sr*L_r)/m_s, k_sf/m_s, c_sf/m_s, k_sr/m_s, c_sr/m_s;
    0, 0, 0, 1, 0, 0, 0, 0;
    (k_sf*L_f - k_sr*L_r)/I_psi, (c_sf*L_f - c_sr*L_r)/I_psi, -(k_sf*L_f^2 + k_sr*L_r^2)/I_psi, -(c_sf*L_f^2 + c_sr*L_r^2)/I_psi, -k_sf*L_f/I_psi, -c_sf*L_f/I_psi, k_sr*L_r/I_psi, c_sr*L_r/I_psi;
    0, 0, 0, 0, 0, 1, 0, 0;
    k_sf/m_uf, c_sf/m_uf, -k_sf*L_f/m_uf, -c_sf*L_f/m_uf, -(k_sf + k_uf)/m_uf, -c_sf/m_uf, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1;
    k_sr/m_ur, c_sr/m_ur, k_sr*L_r/m_ur, c_sr*L_r/m_ur, 0, 0, -(k_sr + k_ur)/m_ur, -c_sr/m_ur
];

% Matrice B (ingressi)
B = [
    0, 0;
    1/m_s, 1/m_s;
    0, 0;
    -L_f/I_psi, L_r/I_psi;
    0, 0;
    -1/m_uf, 0;
    0, 0;
    0, -1/m_ur
];

% Matrice E (disturbi esterni, profili stradali)
E = [
    0, 0;
    0, 0;
    0, 0;
    0, 0;
    0, 0;
    k_uf/m_uf, 0;
    0, 0;
    0, k_ur/m_ur
];

%C=eye(8);
C=[
   1 0 0 0 0 0 0 0;   
   0 0 1 0 0 0 0 0; 
   0 0 0 0 1 0 0 0;
   0 0 0 0 0 0 1 0;];


C1=[1 0 0 0 0 0 0 0;
   0 1 0 0 0 0 0 0;
   0 0 1 0 0 0 0 0;
   0 0 0 1 0 0 0 0;
   0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0];

rho=0.00001;
rhoh2=0.000001;


D32=[0 0;
     0 0;
     0 0;
     0 0;
     rho 0;
     0 rho];

D31=zeros(6,2);

D=zeros(4,2);

x0 = zeros(8,1); % Tutti gli stati inizializzati a zero

t = 0:0.01:10;             % Tempo di simulazione
u1 = zeros(size(t));        % Segnale costante
u2 = zeros(size(t));        % Segnale costante
u = [t' u1' u2'];          % Crea una matrice [tempo, ingresso1, ingresso2]

% [disturbances,z_rf,z_rr]=strada_dissestata();
[disturbances,z_rf,z_rr]=dosso();

[ Ak, Bk, Ck, Dk ] = L1_dinamico(A, B, E, C1, D31,D32, C, D,2)


k_l1=L1_Regulation(A,B,E,C1,D31,D32,2);

% Simulazione
result_closed=sim("L1_dinamicoo.slx");

% Plot delle prime 4 uscite L1 vs pole placement
figure()
plot(t, z_rf,'--','LineWidth',1.5,'DisplayName','Bump');
hold on
plot(result_closed.y_l1.Time, squeeze(result_closed.y_l1.Data(1,1,:)),'LineWidth',1.5,'DisplayName','static L1 ');
xlabel('Time [s]');
ylabel('diplacement [m]');
title('Body displacement');
hold on
plot(result_closed.y_l1_dinamico.Time, squeeze(result_closed.y_l1_dinamico.Data(1,1,:)),'LineWidth',1.5,'DisplayName','dynamic L1');
hold off
legend show
grid on;

figure()
plot(result_closed.y_l1.Time, squeeze(result_closed.y_l1.Data(2,1,:)),'LineWidth',1.5,'DisplayName','static L1');
xlabel('Time [s]');
ylabel('velocity [m/s]');
title('Body velocity');
hold on
plot(result_closed.y_l1_dinamico.Time, squeeze(result_closed.y_l1_dinamico.Data(2,1,:)),'LineWidth',1.5,'DisplayName','dynamic L1');
hold off
legend show
grid on;

figure()
plot(result_closed.y_l1.Time, squeeze(result_closed.y_l1.Data(3,1,:)),'LineWidth',1.5,'DisplayName','static L1 ');
xlabel('Time [s]');
ylabel('angle [rad]');
title('pitch angle');
hold on
plot(result_closed.y_l1_dinamico.Time, squeeze(result_closed.y_l1_dinamico.Data(3,1,:)),'LineWidth',1.5,'DisplayName','dynamic L1');
hold off
legend show
grid on;

figure()
plot(result_closed.y_l1.Time, squeeze(result_closed.y_l1.Data(4,1,:)),'LineWidth',1.5,'DisplayName','static L1');
xlabel('Time [s]');
ylabel('velocity [rad/s]');
title('pitch velocity');
hold on
plot(result_closed.y_l1_dinamico.Time, squeeze(result_closed.y_l1_dinamico.Data(4,1,:)),'LineWidth',1.5,'DisplayName','dynamic L1');
hold off
legend show
grid on;

figure()
plot(result_closed.u_l1.Time, squeeze(result_closed.u_l1.Data(1,1,:)),'LineWidth',1.5,'DisplayName','static L1');
xlabel('Time [s]');
ylabel('force [N]');
title('first control input');
hold on
plot(result_closed.u_l1_dinamico.Time, squeeze(result_closed.u_l1_dinamico.Data(1,1,:)),'LineWidth',1.5,'DisplayName','dynamic L1');
hold off
legend show
grid on;

figure()
plot(result_closed.u_l1.Time, squeeze(result_closed.u_l1.Data(2,1,:)),'LineWidth',1.5,'DisplayName','static L1');
xlabel('Time [s]');
ylabel('force [N]');
title('second control input');
hold on
plot(result_closed.u_l1_dinamico.Time, squeeze(result_closed.u_l1_dinamico.Data(2,1,:)),'LineWidth',1.5,'DisplayName','dynamic L1');
hold off
legend show
grid on;