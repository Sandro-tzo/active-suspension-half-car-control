clear; close all;

% Parametri del sistema
m_s = 690;    % Massa del corpo principale [kg]
I_psi = 1222;  % Momento d'inerzia rispetto all'asse [kg*m^2]
m_uf = 40;     % Massa della ruota anteriore [kg]
m_ur = 45;     % Massa della ruota posteriore [kg]
k_uf = 200000; % Costante elastica ruota anteriore [N/m]
k_ur = 200000; % Costante elastica ruota posteriore [N/m]
L_f = 1.3;     % Distanza tra CG e asse anteriore [m]
L_r = 1.5;     % Distanza tra CG e asse posteriore [m]
k_sf_confort=15000;
k_sr_confort=19000;
c_sf_confort=800;
c_sr_confort=800;
k_sf_normal=18000;
k_sr_normal=22000;
c_sf_normal=1000;
c_sr_normal=1000;
k_sf_sport=21000;
k_sr_sport=25000;
c_sf_sport=1300;
c_sr_sport=1300;

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

rhoh2=0.000001;

D22=[0 0;
     0 0;
     0 0;
     0 0;
     rhoh2 0;
     0 rhoh2];

D11=zeros(6,2);


[disturbances,z_rf,z_rr]=disturbo_stateflow();


A1=calcolaAssetto(k_sf_confort,k_sr_confort,c_sf_confort,c_sr_confort);


A2=calcolaAssetto(k_sf_normal,k_sr_normal,c_sf_normal,c_sr_normal);


A3=calcolaAssetto(k_sf_sport,k_sr_sport,c_sf_sport,c_sr_sport);


k1=H2H8_Regulation(A1,B,E,C22,D22,D11)

k2=H2H8_Regulation(A2,B,E,C22,D22,D11)

k3=H2H8_Regulation(A3,B,E,C22,D22,D11)

sigma=generateSigmaSignal();

t = 0:0.01:20;        
% Simulazione
result_closed=sim("State_flow_GS.slx");


% Plot delle prime 4 uscite 
figure()
plot(t, z_rf,'--','LineWidth',1.5,'DisplayName','road');
hold on
plot(result_closed.y_stateflow.Time, squeeze(result_closed.y_stateflow.Data(1,1,:)),'LineWidth',1.5,'DisplayName','H2Hinf');
xlabel('Time [s]');
ylabel('diplacement [m]');
title('Body displacement');
legend show
grid on;

figure()
plot(result_closed.y_stateflow.Time, squeeze(result_closed.y_stateflow.Data(2,1,:)),'LineWidth',1.5,'DisplayName','H2Hinf');
xlabel('Time [s]');
ylabel('velocity [m/s]');
title('Body velocity');
legend show
grid on;

figure()
plot(result_closed.y_stateflow.Time, squeeze(result_closed.y_stateflow.Data(3,1,:)),'LineWidth',1.5,'DisplayName','H2Hinf');
xlabel('Time [s]');
ylabel('angle [rad]');
title('pitch angle');
legend show
grid on;

figure()
plot(result_closed.y_stateflow.Time, squeeze(result_closed.y_stateflow.Data(4,1,:)),'LineWidth',1.5,'DisplayName','H2Hinf');
xlabel('Time [s]');
ylabel('velocity [rad/s]');
title('pitch velocity');
legend show
grid on;

figure()
plot(result_closed.u_stateflow.Time, squeeze(result_closed.u_stateflow.Data(1,1,:)),'LineWidth',1.5,'DisplayName','H2Hinf');
xlabel('Time [s]');
ylabel('force [N]');
title('first control input');
legend show
grid on;

figure()
plot(result_closed.u_stateflow.Time, squeeze(result_closed.u_stateflow.Data(2,1,:)),'LineWidth',1.5,'DisplayName','H2Hinf');
xlabel('Time [s]');
ylabel('force [N]');
title('second control input');
legend show
grid on;