 function [ Ak, Bk, Ck, Dk ] = L1_dinamico(A,B1, B2, C3, D31,D32, C, D, lambda)

    
    nx = size(A, 1);   
    nu = size(B1, 2);   
    nw = size(B2, 2);   
    nz = size(C3, 1);   
    np = size(C, 1);    

    % Variabili LMI
    X = sdpvar(nx, nx);
    Y = sdpvar(nx, nx);
    mu = sdpvar(1, 1);
    Ack = sdpvar(nx, nx);
    Bck = sdpvar(nx, np);
    Cck = sdpvar(nu, nx);
    Dck = sdpvar(nu, np);
    gamma = sdpvar(1, 1);

    %Vincolo LMI
    V21 = (Ack + (A + B1 * Dck*C)' + lambda * eye(nx));
    V31 = (B2 + B1 * Dck * D)';
    V32 = (X * B2 + Bck * D)';

    %
    V1 = [  (A * Y + B1 * Cck) + (A * Y + B1 * Cck)' + lambda * Y, V21', V31';
             V21, (X * A + Bck * C) + (X * A + Bck * C)' + lambda * X, V32';
              V31, V32, -mu * eye(nw)]<=0;

  
    V2 = [ lambda * Y, lambda * eye(nx), zeros(nx, nw), (C3 * Y + D32 * Cck)';
           lambda * eye(nx), lambda * X, zeros(nx, nw), (C3 + D32 * Dck * C)';
           zeros(nw, nx), zeros(nw, nx), (gamma - mu) * eye(nw), (D31 + D32 * Dck * D)';
           (C3 * Y + D32 * Cck), (C3 + D32 * Dck * C), (D31 + D32 * Dck * D), gamma * eye(nz) ] >= 0;

    V3 = gamma >= 0;

    V4=( X>=0);

    V5=(mu>=0);

    % Ottimizzazione
    V = [V1, V2, V3,V4,V5];
    opts = sdpsettings('solver', 'lmilab', 'verbose', 1);
    optimize(V, gamma, opts);

    % Risultati
    X1 = double(X);
    Y1 = double(Y);
    Ack1 = double(Ack);
    Bck1 = double(Bck);
    Cck1 = double(Cck);
    Dck1 = double(Dck);

    %
    [N, MT] = lu(eye(nx) - X1 * Y1);
    Dk = Dck1;
    Ck = (Cck1 - Dk * C * Y1) * inv(MT);
    Bk = inv(N) * (Bck1 - X1 * B1 * Dk);
    Ak = inv(N) * (Ack1 - N * Bk * C * Y1 - X1 * B1 * Ck * MT - X1 * (A + B1 * Dk * C) * Y1) * inv(MT);

end

