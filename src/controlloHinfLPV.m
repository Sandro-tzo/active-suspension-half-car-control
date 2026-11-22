clear; close all;

% Parametri del sistema
m_s = 690;    % Massa del corpo principale [kg]
I_psi = 1222;  % Momento d'inerzia rispetto all'asse [kg*m^2]
m_uf = 40;     % Massa della ruota anteriore [kg]
m_ur = 45;     % Massa della ruota posteriore [kg]
k_sf_max = 20700; 
k_sf_min = 15300;  % Costante elastica sospensione anteriore [N/m]
k_sr_max = 25300;
k_sr_min = 18700;  % Costante elastica sospensione posteriore [N/m]
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
A1 = calcolaA(k_sf_min,k_sr_min);

A2 = calcolaA(k_sf_min,k_sr_max);

A3 = calcolaA(k_sf_max,k_sr_max);

A4 = calcolaA(k_sf_max,k_sr_min);

A = calcolaA(18000,22000);



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


C=eye(8);
D=zeros(8,4);

x0 = zeros(8,1);

C22=[1 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0;
    0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0];

rho=0.00001;

D22=[0 0;
     0 0;
     0 0;
     0 0;
     rho 0;
     0 rho];

D11=zeros(6,2);

% [disturbances,z_rf,z_rr]=dosso();

[disturbances,z_rf,z_rr]=strada_dissestata();

[k1,k2,k3,k4] = HinfLPV_Regulation(A1,A2,A3,A4,B,B,B,B,E,E,E,E,C22,D22,D11)

k_hinf=Hinf_Regulation(A,B,E,C22,D22,D11);

t = 0:0.01:10;        
% Simulazione
result_closed=sim("controllooHinfLPV.slx");

% Plot delle prime 4 uscite 
figure()
plot(t, z_rf,'--','LineWidth',1.5,'DisplayName','Bump');
hold on
plot(result_closed.y_hinf.Time, squeeze(result_closed.y_hinf.Data(1,1,:)),'LineWidth',1.5,'DisplayName','static Hinf');
xlabel('Time [s]');
ylabel('diplacement [m]');
title('Body displacement');
hold on
plot(result_closed.y_hinf_gs.Time, squeeze(result_closed.y_hinf_gs.Data(1,1,:)),'LineWidth',1.5,'DisplayName','Gain Scheduling Hinf ');
hold off
legend show
grid on;

figure()
plot(result_closed.y_hinf.Time, squeeze(result_closed.y_hinf.Data(2,1,:)),'LineWidth',1.5,'DisplayName','static Hinf');
xlabel('Time [s]');
ylabel('velocity [m/s]');
title('Body velocity');
hold on
plot(result_closed.y_hinf_gs.Time, squeeze(result_closed.y_hinf_gs.Data(2,1,:)),'LineWidth',1.5,'DisplayName','Gain Scheduling Hinf ');
hold off
legend show
grid on;

figure()
plot(result_closed.y_hinf.Time, squeeze(result_closed.y_hinf.Data(3,1,:)),'LineWidth',1.5,'DisplayName','static Hinf');
xlabel('Time [s]');
ylabel('angle [rad]');
title('pitch angle');
hold on
plot(result_closed.y_hinf_gs.Time, squeeze(result_closed.y_hinf_gs.Data(3,1,:)),'LineWidth',1.5,'DisplayName','Gain Scheduling Hinf ');
hold off
legend show
grid on;

figure()
plot(result_closed.y_hinf.Time, squeeze(result_closed.y_hinf.Data(4,1,:)),'LineWidth',1.5,'DisplayName','static Hinf');
xlabel('Time [s]');
ylabel('velocity [rad/s]');
title('pitch velocity');
hold on
plot(result_closed.y_hinf_gs.Time, squeeze(result_closed.y_hinf_gs.Data(4,1,:)),'LineWidth',1.5,'DisplayName','Gain Scheduling Hinf ');
hold off
legend show
grid on;

figure()
plot(result_closed.u_hinf.Time, squeeze(result_closed.u_hinf.Data(1,1,:)),'LineWidth',1.5,'DisplayName','static Hinf');
xlabel('Time [s]');
ylabel('force [N]');
title('first control input');
hold on
plot(result_closed.u_hinf_gs.Time, squeeze(result_closed.u_hinf_gs.Data(1,1,:)),'LineWidth',1.5,'DisplayName','Gain Scheduling Hinf');
hold off
legend show
grid on;

figure()
plot(result_closed.u_hinf.Time, squeeze(result_closed.u_hinf.Data(2,1,:)),'LineWidth',1.5,'DisplayName','static Hinf');
xlabel('Time [s]');
ylabel('force [N]');
title('second control input');
hold on
plot(result_closed.u_hinf_gs.Time, squeeze(result_closed.u_hinf_gs.Data(2,1,:)),'LineWidth',1.5,'DisplayName','Gain Scheduling Hinf');
hold off
legend show
grid on;



