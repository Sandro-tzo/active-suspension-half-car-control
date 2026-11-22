function [x_dot,y] = nonlinear_system(t, x, u, z_rf, z_rr, params)
% HALF_CAR_NONLINEAR Modello non lineare di un sistema "half-car" nello spazio di stato.
%
% Input:
%   t      - Tempo (non usato nel calcolo, ma utile per simulazioni)
%   x      - Stato del sistema [x1; x2; x3; x4; x5; x6; x7; x8]
%            x1: z_s    (posizione verticale della massa sospesa)
%            x2: z_s_dot (velocità verticale della massa sospesa)
%            x3: psi    (angolo di beccheggio)
%            x4: psi_dot (velocità angolare di beccheggio)
%            x5: z_uf   (posizione verticale ruota anteriore)
%            x6: z_uf_dot (velocità verticale ruota anteriore)
%            x7: z_ur   (posizione verticale ruota posteriore)
%            x8: z_ur_dot (velocità verticale ruota posteriore)
%   u      - Ingressi [u_f; u_r] (forze applicate alle ruote anteriore e posteriore)
%   z_rf   - Profilo stradale ruota anteriore
%   z_rr   - Profilo stradale ruota posteriore
%   params - Parametri del sistema (struttura con campi):
%            m_s, I_psi, L_f, L_r, k_sf, c_sf, k_sr, c_sr, m_uf, m_ur, k_uf, k_ur
%
% Output:
%   x_dot - Derivate dello stato [x1_dot; x2_dot; ...; x8_dot]

% Estrai parametri
m_s = params.m_s;       % Massa sospesa
I_psi = params.I_psi;   % Momento d'inerzia
L_f = params.L_f;       % Distanza asse anteriore - CG
L_r = params.L_r;       % Distanza asse posteriore - CG
k_sf = params.k_sf;     % Costante elastica sospensione anteriore
c_sf = params.c_sf;     % Smorzamento sospensione anteriore
k_sr = params.k_sr;     % Costante elastica sospensione posteriore
c_sr = params.c_sr;     % Smorzamento sospensione posteriore
m_uf = params.m_uf;     % Massa non sospesa anteriore
m_ur = params.m_ur;     % Massa non sospesa posteriore
k_uf = params.k_uf;     % Rigidezza pneumatici anteriori
k_ur = params.k_ur;     % Rigidezza pneumatici posteriori

% Estrai stati
z_s = x(1);      z_s_dot = x(2);
psi = x(3);      psi_dot = x(4);
z_uf = x(5);     z_uf_dot = x(6);
z_ur = x(7);     z_ur_dot = x(8);

% Estrai ingressi
u_f = u(1);  % Forza ruota anteriore
u_r = u(2);  % Forza ruota posteriore

% Calcola z_sf e z_sr (posizioni cinematiche)
z_sf = z_s - L_f * sin(psi);
z_sr = z_s + L_r * sin(psi);

% Calcola z_sf_dot e z_sr_dot (velocità cinematiche)
z_sf_dot = z_s_dot - L_f * cos(psi) * psi_dot;
z_sr_dot = z_s_dot + L_r * cos(psi) * psi_dot;

% Equazioni del sistema nello spazio di stato
z_s_ddot = (1/m_s) * ( ...
    -k_sf * (z_sf - z_uf) - c_sf * (z_sf_dot - z_uf_dot) ...
    -k_sr * (z_sr - z_ur) - c_sr * (z_sr_dot - z_ur_dot) ...
    + u_f + u_r);

psi_ddot = (1/I_psi) * ( ...
    L_f * cos(psi) * (k_sf * (z_sf - z_uf) + c_sf * (z_sf_dot - z_uf_dot)) ...
    - L_r * cos(psi) * (k_sr * (z_sr - z_ur) + c_sr * (z_sr_dot - z_ur_dot)) ...
    - L_f * cos(psi) * u_f + L_r * cos(psi) * u_r);

z_uf_ddot = (1/m_uf) * ( ...
    -k_uf * (z_uf - z_rf) - k_sf * (z_uf - z_sf) - c_sf * (z_uf_dot - z_sf_dot) ...
    - u_f);

z_ur_ddot = (1/m_ur) * ( ...
    -k_ur * (z_ur - z_rr) - k_sr * (z_ur - z_sr) - c_sr * (z_ur_dot - z_sr_dot) ...
    - u_r);

% Stato derivato
x_dot = [z_s_dot;           % x1_dot
         z_s_ddot;          % x2_dot
         psi_dot;           % x3_dot
         psi_ddot;          % x4_dot
         z_uf_dot;          % x5_dot
         z_uf_ddot;         % x6_dot
         z_ur_dot;          % x7_dot
         z_ur_ddot];        % x8_dot

y=x;
end
