function C = randConvexCombination(M)
    % Numero di matrici nella cella
    n = numel(M);

    % Genera pesi casuali
    weights = rand(1, n);

    % Normalizza i pesi affinché la loro somma sia 1
    weights = weights / sum(weights);

    % Dimensione delle matrici (si assume che tutte siano uguali)
    C = zeros(size(M{1}));

    % Calcola la combinazione convessa
    for i = 1:n
        C = C + weights(i) * M{i};
    end
end
