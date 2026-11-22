function [disturbances,z_rf,z_rr] = disturbo_stateflow()
% Parametri della strada
G_n0 = 25.6e-3;      % Valore G(n0) per amplificare il profilo stradale
n0 = 0.1;            % Frequenza di riferimento (cicli/m)
n_max = 512;         % Frequenza massima (cicli/m)
delta_n = 0.1;       % Intervallo frequenziale
v = 40 / 3.6;        % Velocità veicolo (m/s)

% Tempo di simulazione
t = (0:0.01:20)';    % Esteso fino a 20s per includere pianura e strada dissestata
x = v * t;           % Posizioni spaziali corrispondenti

% Calcolo dello spettro G_d
n = (delta_n:delta_n:n_max)';         % Vettore delle frequenze
G_d = G_n0 * (n / n0).^(-2);          % Densità spettrale di potenza

% Calcolo del profilo stradale dissestato
z_rf_rough = road_profile(x, n_max, delta_n, G_d, v); % Disturbo ruota anteriore

% Creazione della pianura (profilo piatto)
pianura_duration = 10; % Durata della pianura in secondi
pianura_steps = pianura_duration / 0.01; % Numero di passi temporali per la pianura
z_rf_flat = zeros(pianura_steps, 1);  % Profilo piatto per la pianura

% Combinazione di pianura e strada dissestata
z_rf = [z_rf_flat; z_rf_rough(1:length(t) - pianura_steps)];

% Introduzione dello sfasamento temporale per la ruota posteriore
L = 2.5;                              % Interasse (metri)
delay = L / v;                        % Ritardo temporale (s)

% Trova l'indice corrispondente al ritardo temporale
delay_steps = round(delay / 0.01);    % Converti il ritardo in numero di passi temporali

% Sfasamento per la ruota posteriore
z_rr = [zeros(delay_steps, 1); z_rf(1:end-delay_steps)];

% Creazione del segnale per Simulink
disturbances.time = t;
disturbances.signals.values = [z_rf, z_rr];  % Valori dei disturbi
disturbances.signals.dimensions = 2;         % Due segnali: anteriore e posteriore

% Plot dei profili stradali
figure;
subplot(2, 1, 1);
plot(t, z_rf, 'b', 'LineWidth', 1.5);
title('Road profile - front wheel');
xlabel('Time [s]');
ylabel('h_{front} [m]');
grid on;

subplot(2, 1, 2);
plot(t, z_rr, 'r', 'LineWidth', 1.5);
title('Road profile - rear wheel');
xlabel('Time [s]');
ylabel('h_{rear} [m]');
grid on;

end
