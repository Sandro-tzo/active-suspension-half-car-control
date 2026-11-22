function K = H2_Regulation(A,B,W,C,D)

% definiamo la LMI per la risoluzione del problema H-infinito variabili del nostro sistema

%% Dimensions
    n = size(A)*[1;0];
    m = size(B)*[0;1];
    mw = size(W)*[0;1];
    p = size(C)*[1;0];

    yalmip('clear')

    gamma = sdpvar(1,1);

    %% Variables
    X = sdpvar(n,n);
    Y = sdpvar(m,n);
    Q = sdpvar(p,p);


%%imposto i vincoli LMI
F1=(X>=0);
F2=(gamma>=0);   
F3=([(A*X+B*Y)+(A*X+B*Y)'   W;
             W'           -eye(mw)]<=0);
F4=([ Q      (C*X+D*Y);
    (C*X+D*Y)'    X]>=0);
F5=(trace(Q)<=gamma);
F6=(Q>=0);
%%unisco i vincoli
F=[F1,F2,F3,F4,F5,F6];
%%risolvo le LMI
%% Solving LMI
    opts=sdpsettings;
    opts.solver='lmilab';
    optimize(F,gamma,opts)
    K = double(Y)/(double(X));

end
