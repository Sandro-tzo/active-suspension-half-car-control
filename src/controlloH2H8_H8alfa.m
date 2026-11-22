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
rho=0.000001;
D = [
    0, 0;
    0, 0;
    0, 0;
    0, 0;
    rho 0;
    0, rho];
Dw = zeros(6,2);

x0 = zeros(8,1); % Tutti gli stati inizializzati a zero

t = 0:0.01:10;             % Tempo di simulazione

 % [disturbances,z_rf,z_rr]=dosso();
[disturbances,z_rf,z_rr]=strada_dissestata();

k_h2hinf=H2H8_Regulation(A,B,E,C,D,Dw);

eig(A+B*k_h2hinf);

k_hinfalfa=HinfRstability_regulation(A,B,E,C,Dw,D,70,80,pi/4);

eig(A+B*k_hinfalfa);


% Simulazione
result_closed=sim("H2H8.slx");

% Plot delle prime 4 uscite
figure()
plot(t, z_rf,'--','LineWidth',1.5,'DisplayName','Uneven road');
hold on
plot(result_closed.y_hinfalfa.Time, squeeze(result_closed.y_hinfalfa.Data(1,1,:)),'LineWidth',1.5,'DisplayName',' Hinf-R ');
xlabel('Time [s]');
ylabel('diplacement [m]');
title('Body displacement');
hold on
plot(result_closed.y_h2hinf.Time, squeeze(result_closed.y_h2hinf.Data(1,1,:)),'LineWidth',1.5,'DisplayName','H2-Hinf ');
hold off
legend show
grid on;

figure()
plot(result_closed.y_hinfalfa.Time, squeeze(result_closed.y_hinfalfa.Data(2,1,:)),'LineWidth',1.5,'DisplayName',' Hinf-R');
xlabel('Time [s]');
ylabel('velocity [m/s]');
title('Body velocity');
hold on
plot(result_closed.y_h2hinf.Time, squeeze(result_closed.y_h2hinf.Data(2,1,:)),'LineWidth',1.5,'DisplayName','H2-Hinf ');
hold off
legend show
grid on;

figure()
plot(result_closed.y_hinfalfa.Time, squeeze(result_closed.y_hinfalfa.Data(3,1,:)),'LineWidth',1.5,'DisplayName',' Hinf-R');
xlabel('Time [s]');
ylabel('angle [rad]');
title('pitch angle');
hold on
plot(result_closed.y_h2hinf.Time, squeeze(result_closed.y_h2hinf.Data(3,1,:)),'LineWidth',1.5,'DisplayName','H2-Hinf ');
hold off
legend show
grid on;

figure()
plot(result_closed.y_hinfalfa.Time, squeeze(result_closed.y_hinfalfa.Data(4,1,:)),'LineWidth',1.5,'DisplayName',' Hinf-R');
xlabel('Time [s]');
ylabel('velocity [rad/s]');
title('pitch velocity');
hold on
plot(result_closed.y_h2hinf.Time, squeeze(result_closed.y_h2hinf.Data(4,1,:)),'LineWidth',1.5,'DisplayName','H2-Hinf ');
hold off
legend show
grid on;


figure()
plot(result_closed.u_hinfalfa.Time, squeeze(result_closed.u_hinfalfa.Data(1,1,:)),'LineWidth',1,'DisplayName',' Hinf-R');
xlabel('Time [s]');
ylabel('force [N]');
title('first control input');
hold on
plot(result_closed.u_h2hinf.Time, squeeze(result_closed.u_h2hinf.Data(1,1,:)),'LineWidth',1.5,'DisplayName','h2hinf ');
hold off
legend show
grid on;

figure()
plot(result_closed.u_hinfalfa.Time, squeeze(result_closed.u_hinfalfa.Data(2,1,:)),'LineWidth',1,'DisplayName',' Hinf-R');
xlabel('Time [s]');
ylabel('force [N]');
title('second control input');
hold on
plot(result_closed.u_h2hinf.Time, squeeze(result_closed.u_h2hinf.Data(2,1,:)),'LineWidth',1.5,'DisplayName','h2hinf');
hold off
legend show
grid on;
