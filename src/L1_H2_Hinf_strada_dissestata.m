clear; close all;

% Parametri del sistema
m_s = 690;    % Massa del corpo principale [kg]
I_psi = 1222;  % Momento d'inerzia rispetto all'asse [kg*m^2]
m_uf = 40;     % Massa della ruota anteriore [kg]
m_ur = 45;     % Massa della ruota posteriore [kg]
k_sf = 18000;  % Costante elastica sospensione anteriore [N/m]
k_sr = 22000;  % Costante elastica sospensione posteriore [N/m]
k_uf = 200000; % Costante elastica ruota anteriore [N/m]
k_ur = 200000; % Costante elastica ruota posteriore [N/m]
c_sf = 1000;   % Costante smorzamento sospensione anteriore [Ns/m]
c_sr = 1000;   % Costante smorzamento sospensione posteriore [Ns/m]
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
C=[1 0 0 0 0 0 0 0;
   0 1 0 0 0 0 0 0;
   0 0 1 0 0 0 0 0;
   0 0 0 1 0 0 0 0;
   0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0];

rho=0.00001;
rhoh2=0.000001;


D12=[0 0;
     0 0;
     0 0;
     0 0;
     rho 0;
     0 rho];

D11=zeros(6,2);

D22=[0 0;
     0 0;
     0 0;
     0 0;
     rhoh2 0;
     0 rhoh2];


x0 = zeros(8,1); % Tutti gli stati inizializzati a zero

t = 0:0.01:10;             % Tempo di simulazione
u1 = zeros(size(t));        % Segnale costante
u2 = zeros(size(t));        % Segnale costante
u = [t' u1' u2'];          % Crea una matrice [tempo, ingresso1, ingresso2]

eig(A);

stabilitaLMI(A);

R = ctrb(A,B);

Rango = rank(R);

O = obsv(A,C);

Rango = rank(O);

[disturbances,z_rf,z_rr]=strada_dissestata();

k = place(A,B,[-90,-90,-80,-80,-70,-70,-60,-60;]);
k_pp=-k;

k_l1=L1_Regulation(A,B,E,C,D11,D12,2);

eig(A+B*k_l1);

k_h2=H2_Regulation(A,B,E,C,D22);

eig(A+B*k_h2);

k_hinf=Hinf_Regulation(A,B,E,C,D12,D11);

eig(A+B*k_hinf);


% Simulazione
result_open=sim("open_loop.slx");
% Simulazione
result_closed=sim("L1_H2_Hinf.slx");


% % Plot delle prime 4 uscite open vs pole placement
% figure()
% plot(t, z_rf,'--','LineWidth',1.5,'DisplayName','Uneven road');
% hold on
% plot(result_open.y.Time, squeeze(result_open.y.Data(1,1,:)),'LineWidth',1.5,'DisplayName','open loop');
% xlabel('Time [s]');
% ylabel('diplacement [m]');
% title('Body displacement');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(1,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_open.y.Time, squeeze(result_open.y.Data(2,1,:)),'LineWidth',1.5,'DisplayName','open loop');
% xlabel('Time [s]');
% ylabel('velocity [m/s]');
% title('Body velocity');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(2,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_open.y.Time, squeeze(result_open.y.Data(3,1,:)),'LineWidth',1.5,'DisplayName','open loop');
% xlabel('Time [s]');
% ylabel('angle [rad]');
% title('pitch angle');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(3,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_open.y.Time, squeeze(result_open.y.Data(4,1,:)),'LineWidth',1.5,'DisplayName','open loop');
% xlabel('Time [s]');
% ylabel('velocity [rad/s]');
% title('pitch velocity');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(4,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;


% % Plot delle prime 4 uscite L1 vs pole placement
% figure()
% plot(t, z_rf,'--','LineWidth',1.5,'DisplayName','Uneven road');
% hold on
% plot(result_closed.y_l1.Time, squeeze(result_closed.y_l1.Data(1,1,:)),'LineWidth',1.5,'DisplayName','L1 regulation');
% xlabel('Time [s]');
% ylabel('diplacement [m]');
% title('Body displacement');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(1,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_closed.y_l1.Time, squeeze(result_closed.y_l1.Data(2,1,:)),'LineWidth',1.5,'DisplayName','L1 regulation');
% xlabel('Time [s]');
% ylabel('velocity [m/s]');
% title('Body velocity');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(2,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_closed.y_l1.Time, squeeze(result_closed.y_l1.Data(3,1,:)),'LineWidth',1.5,'DisplayName','L1 regulation');
% xlabel('Time [s]');
% ylabel('angle [rad]');
% title('pitch angle');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(3,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_closed.y_l1.Time, squeeze(result_closed.y_l1.Data(4,1,:)),'LineWidth',1.5,'DisplayName','L1 regulation');
% xlabel('Time [s]');
% ylabel('velocity [rad/s]');
% title('pitch velocity');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(4,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_closed.u_l1.Time, squeeze(result_closed.u_l1.Data(1,1,:)),'LineWidth',1.5,'DisplayName','L1 regulation');
% xlabel('Time [s]');
% ylabel('force [N]');
% title('first control input');
% hold on
% plot(result_closed.u_pp.Time, squeeze(result_closed.u_pp.Data(1,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_closed.u_l1.Time, squeeze(result_closed.u_l1.Data(2,1,:)),'LineWidth',1.5,'DisplayName','L1 regulation');
% xlabel('Time [s]');
% ylabel('force [N]');
% title('second control input');
% hold on
% plot(result_closed.u_pp.Time, squeeze(result_closed.u_pp.Data(2,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;

% % Plot delle prime 4 uscite H2 vs pole placement
% figure()
% plot(t, z_rf,'--','LineWidth',1.5,'DisplayName','Uneven road');
% hold on
% plot(result_closed.y_h2.Time, squeeze(result_closed.y_h2.Data(1,1,:)),'LineWidth',1.5,'DisplayName','H2 regulation');
% xlabel('Time [s]');
% ylabel('diplacement [m]');
% title('Body displacement');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(1,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_closed.y_h2.Time, squeeze(result_closed.y_h2.Data(2,1,:)),'LineWidth',1.5,'DisplayName','H2 regulation');
% xlabel('Time [s]');
% ylabel('velocity [m/s]');
% title('Body velocity');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(2,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_closed.y_h2.Time, squeeze(result_closed.y_h2.Data(3,1,:)),'LineWidth',1.5,'DisplayName','H2 regulation');
% xlabel('Time [s]');
% ylabel('angle [rad]');
% title('pitch angle');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(3,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_closed.y_h2.Time, squeeze(result_closed.y_h2.Data(4,1,:)),'LineWidth',1.5,'DisplayName','H2 regulation');
% xlabel('Time [s]');
% ylabel('velocity [rad/s]');
% title('pitch velocity');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(4,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_closed.u_h2.Time, squeeze(result_closed.u_h2.Data(1,1,:)),'LineWidth',1.5,'DisplayName','H2 regulation');
% xlabel('Time [s]');
% ylabel('force [N]');
% title('first control input');
% hold on
% plot(result_closed.u_pp.Time, squeeze(result_closed.u_pp.Data(1,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_closed.u_h2.Time, squeeze(result_closed.u_h2.Data(2,1,:)),'LineWidth',1.5,'DisplayName','H2 regulation');
% xlabel('Time [s]');
% ylabel('force [N]');
% title('second control input');
% hold on
% plot(result_closed.u_pp.Time, squeeze(result_closed.u_pp.Data(2,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;

% % Plot delle prime 4 uscite Hinf vs pole placement
% figure()
% plot(t, z_rf,'--','LineWidth',1.5,'DisplayName','Uneven road');
% hold on
% plot(result_closed.y_hinf.Time, squeeze(result_closed.y_hinf.Data(1,1,:)),'LineWidth',1.5,'DisplayName','Hinf regulation');
% xlabel('Time [s]');
% ylabel('diplacement [m]');
% title('Body displacement');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(1,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_closed.y_hinf.Time, squeeze(result_closed.y_hinf.Data(2,1,:)),'LineWidth',1.5,'DisplayName','Hinf regulation');
% xlabel('Time [s]');
% ylabel('velocity [m/s]');
% title('Body velocity');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(2,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_closed.y_hinf.Time, squeeze(result_closed.y_hinf.Data(3,1,:)),'LineWidth',1.5,'DisplayName','Hinf regulation');
% xlabel('Time [s]');
% ylabel('angle [rad]');
% title('pitch angle');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(3,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_closed.y_hinf.Time, squeeze(result_closed.y_hinf.Data(4,1,:)),'LineWidth',1.5,'DisplayName','Hinf regulation');
% xlabel('Time [s]');
% ylabel('velocity [rad/s]');
% title('pitch velocity');
% hold on
% plot(result_closed.y_pp.Time, squeeze(result_closed.y_pp.Data(4,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_closed.u_hinf.Time, squeeze(result_closed.u_hinf.Data(1,1,:)),'LineWidth',1.5,'DisplayName','Hinf regulation');
% xlabel('Time [s]');
% ylabel('force [N]');
% title('first control input');
% hold on
% plot(result_closed.u_pp.Time, squeeze(result_closed.u_pp.Data(1,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;
% 
% figure()
% plot(result_closed.u_hinf.Time, squeeze(result_closed.u_hinf.Data(2,1,:)),'LineWidth',1.5,'DisplayName','Hinf regulation');
% xlabel('Time [s]');
% ylabel('force [N]');
% title('second control input');
% hold on
% plot(result_closed.u_pp.Time, squeeze(result_closed.u_pp.Data(2,1,:)),'LineWidth',1.5,'DisplayName','pole placement');
% hold off
% legend show
% grid on;


% Plot delle prime quattro uscite L1 vs H2 vs Hinf
figure()
plot(t, z_rf,'--','LineWidth',1.5,'DisplayName','Uneven road');
hold on
plot(result_closed.y_pp.Time, squeeze(result_closed.y_hinf.Data(1,1,:)),'LineWidth',1.5,'DisplayName','H inf');
xlabel('Time [s]');
ylabel('displacement [m]');
title('Body displacement');
hold on
plot(result_closed.y_pp.Time, squeeze(result_closed.y_h2.Data(1,1,:)),'LineWidth',1.5,'DisplayName','H2');
hold on
plot(result_closed.y_pp.Time, squeeze(result_closed.y_l1.Data(1,1,:)),'LineWidth',1.5,'DisplayName','L1');
hold off
legend show
grid on;

figure()
plot(result_closed.y_pp.Time, squeeze(result_closed.y_hinf.Data(2,1,:)),'LineWidth',1.5,'DisplayName','H inf');
xlabel('Time [s]');
ylabel('velocity [m/s]');
title('Body velocity');
hold on
plot(result_closed.y_pp.Time, squeeze(result_closed.y_h2.Data(2,1,:)),'LineWidth',1.5,'DisplayName','H2');
hold on
plot(result_closed.y_pp.Time, squeeze(result_closed.y_l1.Data(2,1,:)),'LineWidth',1.5,'DisplayName','L1');
hold off
legend show
grid on;

figure()
plot(result_closed.y_pp.Time, squeeze(result_closed.y_hinf.Data(3,1,:)),'LineWidth',1.5,'DisplayName','H inf');
xlabel('Time [s]');
ylabel('angle [rad]');
title('Pitch angle');
hold on
plot(result_closed.y_pp.Time, squeeze(result_closed.y_h2.Data(3,1,:)),'LineWidth',1.5,'DisplayName','H2');
hold on
plot(result_closed.y_pp.Time, squeeze(result_closed.y_l1.Data(3,1,:)),'LineWidth',1.5,'DisplayName','L1');
hold off
legend show
grid on;


figure()
plot(result_closed.y_pp.Time, squeeze(result_closed.y_hinf.Data(4,1,:)),'LineWidth',1.5,'DisplayName','H inf');
xlabel('Time [s]');
ylabel('velocity [rad/s]');
title('Pitch velocity');
hold on
plot(result_closed.y_pp.Time, squeeze(result_closed.y_h2.Data(4,1,:)),'LineWidth',1.5,'DisplayName','H2');
hold on
plot(result_closed.y_pp.Time, squeeze(result_closed.y_l1.Data(4,1,:)),'LineWidth',1.5,'DisplayName','L1');
hold off
legend show
grid on;

figure()
plot(result_closed.u_l1.Time, squeeze(result_closed.u_l1.Data(1,1,:)),'LineWidth',1.5,'DisplayName','L1 regulation');
xlabel('Time [s]');
ylabel('force [N]');
title('first control input');
hold on
plot(result_closed.u_h2.Time, squeeze(result_closed.u_h2.Data(1,1,:)),'LineWidth',1.5,'DisplayName','H2 regulation');
hold on
plot(result_closed.u_hinf.Time, squeeze(result_closed.u_hinf.Data(1,1,:)),'LineWidth',1.5,'DisplayName','Hinf regulation');
hold off
legend show
grid on;

figure()
plot(result_closed.u_l1.Time, squeeze(result_closed.u_l1.Data(2,1,:)),'LineWidth',1.5,'DisplayName','L1 regulation');
xlabel('Time [s]');
ylabel('force [N]');
title('second control input');
hold on
plot(result_closed.u_h2.Time, squeeze(result_closed.u_h2.Data(2,1,:)),'LineWidth',1.5,'DisplayName','H2 regulation');
hold on
plot(result_closed.u_hinf.Time, squeeze(result_closed.u_hinf.Data(2,1,:)),'LineWidth',1.5,'DisplayName','Hinf regulation');
hold off
legend show
grid on;
