function K = L1_Regulation(A,B1,B2,C,D31,D32,lambda)

% definiamo la LMI per la risoluzione del problema L1

%% Dimensions
    n = size(A)*[1;0];
    m = size(B1)*[0;1];
    mw = size(B2)*[0;1];
    p = size(C)*[1;0];

    yalmip('clear')

    gamma = sdpvar(1,1);
    mu=sdpvar(1,1);

    %% Variables
    X = sdpvar(n,n);
    Y = sdpvar(m,n);
    

%%imposto i vincoli LMI
F1=(X>=0);
F2=(gamma>=0);   
F3=([(A*X+B1*Y)+(A*X+B1*Y)'+lambda*X   B2;
             B2'           -mu*eye(mw)]<=0);
F4=([ lambda*X      zeros(n,mw)   (C*X+D32*Y)';
       zeros(mw,n)  (gamma-mu)*eye(mw)   D31';
    (C*X+D32*Y)    D31     gamma*eye(p)]>=0);
F5=(mu>=0);
%%unisco i vincoli
F=[F1,F2,F3,F4,F5];
%%risolvo le LMI
%% Solving LMI
    opts=sdpsettings;
    opts.solver='lmilab';
    optimize(F,gamma,opts)
    K = double(Y)/(double(X));

end