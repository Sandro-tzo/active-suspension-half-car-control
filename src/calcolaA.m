function A=calcolaA(k_sf,k_sr)

m_s = 690;    % Massa del corpo principale [kg]
I_psi = 1222;  % Momento d'inerzia rispetto all'asse [kg*m^2]
m_uf = 40;     % Massa della ruota anteriore [kg]
m_ur = 45;     % Massa della ruota posteriore [kg]
k_uf = 200000; % Costante elastica ruota anteriore [N/m]
k_ur = 200000; % Costante elastica ruota posteriore [N/m]
c_sf = 1000;   % Costante smorzamento sospensione anteriore [Ns/m]
c_sr = 1000;   % Costante smorzamento sospensione posteriore [Ns/m]
L_f = 1.3;     % Distanza tra CG e asse anteriore [m]
L_r = 1.5;     % Distanza tra CG e asse posteriore [m]

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

end
