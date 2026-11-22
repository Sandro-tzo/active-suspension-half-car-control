function [K1,K2,K3,K4] = HinfLPV_Regulation(A1, A2, A3, A4, B1, B2, B3, B4, W1, W2, W3, W4, C, D,Dw)
    %% Dimensions
    n = size(A1, 1);  
    m = size(B1, 2);  
    mw = size(W1, 2); 
    p = size(C, 1);   

    yalmip('clear')

    %% Variables
    gamma = sdpvar(1, 1);
    X = sdpvar(n, n, 'symmetric'); 
    Y1 = sdpvar(m, n);
    Y2 = sdpvar(m, n);
    Y3 = sdpvar(m, n);
    Y4 = sdpvar(m, n);
    Q = sdpvar(p, p, 'symmetric');  

    A = {A1, A2, A3, A4};
    B = {B1, B2, B3, B4};
    W = {W1, W2, W3, W4};
    Y = {Y1, Y2, Y3, Y4};

    %% Vincoli LMI

     con = [];  

 
for i=1:4
    for j=1:4
    con=[con,[(A{i}*X+B{i}*Y{j})+(A{i}*X+B{i}*Y{j})'        W{i}           (C*X+D*Y{j})';
                W{i}'             -gamma*eye(mw)        Dw';
             (C*X+D*Y{j})               Dw          -gamma*eye(p)]<=0];
    end
end


    con=[con, X >= 0];         
    con=[con,gamma >= 0];      


    %% Risoluzione delle LMI
    %%risolvo le LMI
opts=sdpsettings;
    opts.solver='lmilab';
    optimize(con,gamma,opts)

    %% Calcolo del guadagno
    K1 = double(Y1) / (double(X));
    K2 = double(Y2) / (double(X));
    K3 = double(Y3) / (double(X));
    K4 = double(Y4) / (double(X));
end