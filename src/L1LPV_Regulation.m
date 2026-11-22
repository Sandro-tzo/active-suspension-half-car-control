function [K1,K2,K3,K4] = L1LPV_Regulation(A1, A2, A3, A4, B1, B2, B3, B4, W1, W2, W3, W4, C, D31,D32,lambda)
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
    mu=sdpvar(1,1);

    A = {A1, A2, A3, A4};
    B = {B1, B2, B3, B4};
    W = {W1, W2, W3, W4};
    Y = {Y1, Y2, Y3, Y4};

    %% Vincoli LMI

     con = [];  

   
        for j = 1:4       
            con= [con,[ lambda*X      zeros(n,mw)   (C*X+D32*Y{j})';
       zeros(mw,n)  (gamma-mu)*eye(mw)   D31';
       (C*X+D32*Y{j})    D31     gamma*eye(p)]>=0];
        end
    
% 
    
     for i = 1:4
        for j=1:4
            con = [con, [(A{i}*X+B{i}*Y{j})+(A{i}*X+B{i}*Y{j})'+lambda*X   W{i};
             W{i}'           -mu*eye(mw)]<=0];
        end
     end

    con=[con, X >= 0];         
    con=[con,gamma >= 0];      
    con=[con,mu >=0]; 
%%unisco i vincoli


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