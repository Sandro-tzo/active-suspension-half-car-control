function h = road_profile(x, n_max, delta_n, G_d, v)
% ROAD_PROFILE Calcola il profilo stradale secondo la formula 2.0.3
%
% Input:
%   x       - Vettore delle posizioni spaziali (x = v * t se v è la velocità)
%   n_max   - Frequenza massima fino alla quale calcolare le armoniche
%   delta_n - Intervallo frequenziale minimo considerato
%   G_d     - Spettro di densità di potenza [G_d(i)] (vettore)
%   v       - Velocità del veicolo (costante)
%
% Output:
%   h       - Profilo stradale calcolato per ciascun valore di x

    % Numero di armoniche
    N = n_max / delta_n;
    
    % Inizializza il profilo stradale h come vettore di zeri
    h = zeros(size(x));
    
    % Ciclo per sommare tutte le armoniche
    for i = 1:N
        % Calcola la frequenza armonica corrente
        f_i = i * delta_n; 
        
        % Calcola l'ampiezza dello spettro per l'armonica corrente
        A_i = sqrt(2 * delta_n * G_d(i));
        
        % Genera una fase casuale tra 0 e 2*pi
        phi = 2 * pi * rand;
        
        % Aggiunge il contributo dell'armonica corrente
        h = h + A_i * cos(2 * pi * f_i * x / v + phi);
    end
end
