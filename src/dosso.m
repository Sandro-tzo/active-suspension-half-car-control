function [disturbances,z_rf,z_rr]= dosso()

% Tempo di simulazione
t = (0:0.01:10)';    % Vettore tempo: da 0 a 10 s con passo 0.01 
 
%dosso
h_dosso = 0.15;              
l_dosso = 4;                
v = 40/3.6;                
d_dosso = 10;             
t1_dosso = d_dosso / v;
t2_dosso = t1_dosso + l_dosso / v;
t_salita = (t2_dosso - t1_dosso) / 4;
t_discesa = t_salita;
t_pian = (t2_dosso - t1_dosso) - t_salita - t_discesa;
t_distdossi=(t2_dosso-t1_dosso)/0.01;
%segnale del dosso
z_rf = zeros(size(t));
z_rr = zeros(size(t));
for i = 1:length(t)
    if t(i) >= t1_dosso && t(i) < t1_dosso + t_salita
        % salita
        z_rr(i+t_distdossi/2) = h_dosso * (t(i) - t1_dosso) / t_salita;
        z_rf(i) = h_dosso * (t(i) - t1_dosso) / t_salita;
    elseif t(i) >= t1_dosso + t_salita && t(i) < t1_dosso + t_salita + t_pian
        % Parte piatta del dosso
        z_rr(i+t_distdossi/2) = h_dosso;
        z_rf(i) = h_dosso;
    elseif t(i) >= t1_dosso + t_salita + t_discesa && t(i) <= t2_dosso
        % discesa
        z_rr(i+t_distdossi/2) = h_dosso * (1 - (t(i) - (t1_dosso + t_salita + t_pian)) / t_discesa);
        z_rf(i) = h_dosso * (1 - (t(i) - (t1_dosso + t_salita + t_pian)) / t_discesa);
    end
end

% Creazione del vettore di disturbo per Simulink
disturbances.time = t;                     % Vettore del tempo
disturbances.signals.values = [z_rf, z_rr]; % Valori del profilo stradale [z_rf, z_rr]
disturbances.signals.dimensions = 2;       % Due segnali (anteriore e posteriore)

figure;
plot(t, z_rf, 'b', 'LineWidth', 1.5); hold on;
plot(t, z_rr, 'r--', 'LineWidth', 1.5);
title('Road Profile - Bump for Front and Rear Wheels');
xlabel('Time [s]');
ylabel('height [m]');
legend('front wheel', 'rear wheel');
grid on;

end